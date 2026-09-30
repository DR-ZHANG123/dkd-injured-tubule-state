"""03_pseudobulk_de: donor-level DESeq2 per lineage / fine cell type in the discovery atlas.

Model ~ batch + group (Control reference); batch = protocol batch from 02c (SC_A / SC_B / SN),
QC-failed libraries removed. Contrasts: dD = DKD-Control, dN = HKD-Control, dE = DKD-HKD.
Modes: all_libraries; one_lib_per_donor (library with most cells); SC_A_only and SN_only
(within-technology strata, used as an internal two-stratum replication).
Outputs results/03_pseudobulk_de/{level}/{celltype}.tsv and screen_summary.tsv
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "03_pseudobulk_de"
CONTRASTS = {"dD": ("DKD", "Control"), "dN": ("HKD", "Control"), "dE": ("DKD", "HKD")}


def eligible(obs: pd.DataFrame, cfg: dict) -> pd.Series:
    return obs.n_cells >= cfg["pseudobulk"]["min_cells_per_sample"]


def run_deseq(counts: pd.DataFrame, meta: pd.DataFrame, cfg: dict, n_cpus: int = 4) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    p = cfg["pseudobulk"]
    keep = (counts >= p["min_counts_gene"]).sum(0) >= p["min_samples_gene"]
    counts = counts.loc[:, keep]
    factors = ["batch", "group"] if meta.batch.nunique() > 1 else ["group"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(counts=counts.astype(int), metadata=meta[factors], design_factors=factors,
                           ref_level=["group", "Control"], refit_cooks=True, n_cpus=n_cpus, quiet=True)
        dds.deseq2()
        res = []
        for name, (a, b) in CONTRASTS.items():
            st = DeseqStats(dds, contrast=["group", a, b], quiet=True)
            st.summary()
            r = st.results_df[["baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]].copy()
            r.columns = ["baseMean"] + [f"{name}_{c}" for c in ["lfc", "se", "p", "padj"]]
            res.append(r if not res else r.drop(columns="baseMean"))
    return pd.concat(res, axis=1)


def one(level: str, ct: str, pb: ad.AnnData, cfg: dict, mode: str):
    o = pb.obs
    m = (o.celltype == ct) & eligible(o, cfg)
    sub = o[m]
    if mode == "one_lib_per_donor":
        sub = sub.sort_values("n_cells", ascending=False).drop_duplicates("donor")
    elif mode in ("SC_A_only", "SN_only"):
        sub = sub[sub.batch == mode.replace("_only", "")]
    nd = sub.groupby("group").donor.nunique()
    if any(nd.get(g, 0) < cfg["pseudobulk"]["min_donors_per_group"] for g in cfg["cohort"]["groups"]):
        return None
    idx = [pb.obs_names.get_loc(i) for i in sub.index]
    X = pb.X[idx]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    counts = pd.DataFrame(X, index=sub.index, columns=pb.var_names)
    res = run_deseq(counts, sub, cfg)
    res.insert(0, "gene", res.index)
    d = stage_dir(STAGE) / f"{level}_{mode}"
    d.mkdir(exist_ok=True)
    res.to_csv(d / f"{ct.replace('/', '_')}.tsv", sep="\t", index=False)
    row = {"level": level, "mode": mode, "celltype": ct, "n_genes": len(res),
           **{f"n_donors_{g}": int(nd.get(g, 0)) for g in cfg["cohort"]["groups"]},
           "n_samples": len(sub)}
    for c in CONTRASTS:
        for thr in [0.05, 0.10]:
            row[f"{c}_nDE_fdr{int(thr*100):02d}"] = int((res[f"{c}_padj"] < thr).sum())
    # shared-injury vs DKD-biased counts (plan rules 2 & 10)
    sig = lambda c: res[f"{c}_padj"] < 0.10
    same = np.sign(res.dD_lfc) == np.sign(res.dN_lfc)
    row["n_shared_injury"] = int((sig("dD") & sig("dN") & same & ~sig("dE")).sum())
    row["n_dkd_biased"] = int((sig("dE") & (np.sign(res.dE_lfc) == np.sign(res.dD_lfc))).sum())
    return row


def main():
    cfg = load_config(); set_global_seed(cfg["seed"])
    proc = ROOT / cfg["paths"]["processed"]
    jobs = []
    pbs = {}
    libqc = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    for level in ["lineage", "fine"]:
        pb = ad.read_h5ad(proc / f"pb_{level}.h5ad")
        pb.obs["batch"] = pb.obs.library.map(libqc.batch).astype(str)
        pbs[level] = pb[pb.obs.library.map(libqc.qc_pass).astype(bool).values].copy()
        for ct in sorted(pbs[level].obs.celltype.unique()):
            if ct == "Other":
                continue
            for mode in ["all_libraries", "one_lib_per_donor", "SC_A_only", "SN_only"]:
                jobs.append((level, ct, mode))
    rows = Parallel(n_jobs=24)(delayed(one)(lv, ct, pbs[lv], cfg, md) for lv, ct, md in jobs)
    summ = pd.DataFrame([r for r in rows if r])
    out = stage_dir(STAGE)
    summ.to_csv(out / "screen_summary.tsv", sep="\t", index=False)
    print(summ.to_string())
    write_provenance(STAGE, [proc / "pb_lineage.h5ad", proc / "pb_fine.h5ad",
                             ROOT / "results/02c_library_qc/library_qc.tsv"],
                     [out / "screen_summary.tsv"] + sorted(out.glob("*/*.tsv")), cfg["seed"])


if __name__ == "__main__":
    main()
