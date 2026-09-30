"""07a_prepare_kpmp: read KPMP open 10x sc/snRNA expression matrices (zip: barcodes/features/
matrix) for the study groups, per-cell QC, attach participant clinical fields and group.

Outputs data/processed/KPMP_{sn,sc}_qc_counts.h5ad and results/07a_prepare_kpmp/library_qc.tsv.
Non-10x (snDrop-seq, reference-only) and non-mtx formats are skipped and listed."""
from __future__ import annotations
import importlib, io, sys, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "07a_prepare_kpmp"
RAW = ROOT / "data/raw/kpmp"


def _member(names, stems):
    for n in names:
        base = n.split("/")[-1].lower()
        if any(base.endswith(f"{st}{v}{ext}") for st in stems for v in ("", "-v2")
               for ext in (".tsv.gz", ".tsv", ".mtx.gz", ".mtx")):
            return n
    raise KeyError(f"{stems} not in archive")


def _read(z, name, **kw):
    raw = z.read(name)
    if name.endswith(".gz"):
        raw = __import__("gzip").decompress(raw)
    return raw if kw.pop("binary", False) else pd.read_csv(io.BytesIO(raw), header=None, sep="\t", **kw)


def read_zip(path: Path):
    """10x triplet in a zip; handles features/genes naming, -v2 suffixes, gz or plain text."""
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        bc = _read(z, _member(names, ["barcodes"]))[0].values
        ft = _read(z, _member(names, ["features", "genes"]))
        M = scipy.io.mmread(io.BytesIO(_read(z, _member(names, ["matrix"]), binary=True))).tocsr().T
    if ft.shape[1] >= 3:
        keep = (ft[2] == "Gene Expression").values
        M, ft = M[:, keep], ft[keep]
    sym = pd.Index(ft[1].astype(str).str.replace(r"^hg38_+", "", regex=True).values)
    # collapse duplicated symbols by summing
    if sym.duplicated().any():
        codes, uniq = pd.factorize(sym)
        agg = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
        M, sym = (M @ agg).tocsr(), pd.Index(uniq)
    return M.astype(np.float32), bc, sym


def barcodes(path: Path) -> set:
    with zipfile.ZipFile(path) as z:
        return set(_read(z, _member(z.namelist(), ["barcodes"]))[0].astype(str))


def dedupe_reprocessed(files: pd.DataFrame, raw: Path, min_overlap: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Archives of one participant sharing > min_overlap of barcodes are the same library
    processed twice; keep the newer ('-v2') archive."""
    drop, log = set(), []
    for pid, g in files.groupby("redcap_id"):
        if len(g) < 2:
            continue
        bcs = {r.file_name: barcodes(raw / r.experimental_strategy.replace(" ", "_") / pid / r.file_name)
               for r in g.itertuples()}
        fs = list(bcs)
        for i in range(len(fs)):
            for j in range(i + 1, len(fs)):
                a, b = bcs[fs[i]], bcs[fs[j]]
                ov = len(a & b) / max(1, min(len(a), len(b)))
                if ov > min_overlap:
                    old = fs[i] if "-v2" not in fs[i] and "-v2" in fs[j] else fs[j] if "-v2" not in fs[j] and "-v2" in fs[i] else sorted([fs[i], fs[j]])[0]
                    drop.add(old)
                    log.append({"participant": pid, "file_a": fs[i], "file_b": fs[j], "overlap": ov, "dropped": old})
    return files[~files.file_name.isin(drop)], pd.DataFrame(log)


def qc_one(row, q):
    try:
        M, bc, sym = read_zip(RAW / row.experimental_strategy.replace(" ", "_") / row.redcap_id / row.file_name)
    except (OSError, KeyError, ValueError, zipfile.BadZipFile) as e:  # corrupted / unexpected archive -> reported, not silently dropped
        return None, {"file": row.file_name, "participant": row.redcap_id, "error": str(e)[:200]}
    ngenes = np.asarray((M > 0).sum(1)).ravel()
    mt = np.asarray(M[:, sym.str.startswith("MT-")].sum(1)).ravel() / np.maximum(np.asarray(M.sum(1)).ravel(), 1) * 100
    keep = (ngenes >= q["min_genes"]) & (ngenes <= q["max_genes"]) & (mt <= q["max_pct_mt"])
    a = ad.AnnData(X=M[keep], var=pd.DataFrame(index=sym))
    lib = row.file_name.split("_")[0][:8]
    a.obs_names = [f"{row.redcap_id}_{lib}_{b}" for b in bc[keep]]
    a.obs["participant"] = row.redcap_id; a.obs["library"] = f"{row.redcap_id}_{lib}"
    return a, {"file": row.file_name, "participant": row.redcap_id, "library": f"{row.redcap_id}_{lib}",
               "cells_raw": int(M.shape[0]), "cells_qc": int(keep.sum()),
               "median_genes_qc": float(np.median(ngenes[keep])) if keep.any() else 0.0}


def main():
    cfg = load_config(); set_global_seed(cfg["seed"])
    groups = importlib.import_module("00b_download_kpmp").participant_groups().set_index("Participant ID")
    man = pd.read_csv(RAW / "download_manifest_expr.tsv", sep="\t")
    man = man[man.experimental_strategy.isin(["Single-nucleus RNA-Seq", "Single-cell RNA-Seq"])
              & ~man.file_name.str.contains("METADATA")]
    use = man[(man.platform == "10x Genomics") & man.file_name.str.endswith(".zip")]
    exists = use.apply(lambda r: (RAW / r.experimental_strategy.replace(" ", "_") / r.redcap_id / r.file_name).exists(), axis=1)
    skipped = pd.concat([man.drop(use.index).assign(reason="non-10x or non-mtx"),
                         use[~exists].assign(reason="download failed")])
    use = use[exists]
    out = stage_dir(STAGE)
    use, dup = dedupe_reprocessed(use, RAW, cfg["kpmp_qc"]["duplicate_barcode_overlap"])
    dup.to_csv(out / "duplicate_archives.tsv", sep="\t", index=False)
    skipped[["redcap_id", "experimental_strategy", "platform", "file_name", "reason"]].to_csv(out / "skipped_files.tsv", sep="\t", index=False)
    qc_rows, dst = [], []
    for strat, tag in [("Single-nucleus RNA-Seq", "sn"), ("Single-cell RNA-Seq", "sc")]:
        q = cfg["kpmp_qc"][tag]
        sub = use[use.experimental_strategy == strat]
        res = Parallel(n_jobs=32)(delayed(qc_one)(r, q) for r in sub.itertuples())
        ads = [a for a, _ in res if a is not None and a.n_obs > 0]
        qc_rows += [dict(m, modality=tag) for _, m in res]
        # genes measured in every library (annotation versions differ between sites)
        genes = sorted(set.intersection(*[set(a.var_names) for a in ads]))
        assert len(genes) > 10000, f"only {len(genes)} shared genes in {tag}"
        A = ad.concat([a[:, genes] for a in ads], join="inner")
        g = groups.reindex(A.obs.participant.values)
        A.obs["group"] = g.group.values
        for col, new in [("Baseline eGFR (ml/min/1.73m2) (Binned)", "egfr_band"), ("Proteinuria (mg) (Binned)", "proteinuria_band"),
                         ("Albuminuria (mg) (Binned)", "albuminuria_band"), ("Sex", "sex"), ("Age (Years) (Binned)", "age_band"),
                         ("Diabetes History", "dm"), ("Hypertension History", "htn"), ("Tissue Source", "tissue_source")]:
            A.obs[new] = g[col].astype(str).values
        f = ROOT / cfg["paths"]["processed"] / f"KPMP_{tag}_qc_counts.h5ad"
        A.write_h5ad(f, compression="gzip"); dst.append(f)
        print(tag, A.shape, A.obs.groupby("group").participant.nunique().to_dict(), flush=True)
    pd.DataFrame(qc_rows).to_csv(out / "library_qc.tsv", sep="\t", index=False)
    write_provenance(STAGE, [RAW / "download_manifest_expr.tsv", RAW / cfg["kpmp"]["clinical_table"]],
                     dst + [out / "library_qc.tsv", out / "skipped_files.tsv", out / "duplicate_archives.tsv"], cfg["seed"])


if __name__ == "__main__":
    main()
