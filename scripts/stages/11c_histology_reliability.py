"""11c_histology_reliability: noise ceiling of the spot-level targets used in 11b (GSE211785, KPMP).

For each section and set, the measured genes are split at random into two halves (200 splits);
the Spearman-Brown corrected correlation between half-set spot scores estimates the reliability
of the full-set spot score. Genes and scoring follow 10_visium_common (per-section z of log1p CP10k)."""
from __future__ import annotations
import importlib, sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "11c_histology_reliability"
SETS = ["DKDiPT_up", "SharedInjury_up"]


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed)
    vc = importlib.import_module("10_visium_common")
    proc = ROOT / cfg["paths"]["processed"]
    X = scipy.io.mmread(proc / "GSE211785_ST_counts.mtx").tocsc().T.tocsr()     # spots x genes
    genes = pd.Index(open(proc / "GSE211785_ST_genes.txt").read().split())
    bcs = pd.Index(open(proc / "GSE211785_ST_barcodes.txt").read().split())
    meta = pd.read_csv(ROOT / "data/raw/GSE211785/GSE211785_ST_metadata.txt.gz", sep="\t", index_col=0)
    sets = vc.locked_sets()
    rng = np.random.default_rng(seed)
    rows = []
    mats = [("gse211785", sec, X[bcs.get_indexer(idx)], genes) for sec, idx in meta.groupby("orig.ident").groups.items()]
    import scanpy as sc
    st = pd.read_csv(ROOT / "results/10b_visium_kpmp/sample_table.tsv", sep="\t")
    done = set(pd.read_csv(ROOT / "results/11a_histology_embed/section_index.tsv", sep="\t").section)
    for r in st.itertuples():
        sec = f"{r.participant}__{r.sample}"
        if sec not in done:
            continue
        a = sc.read_10x_h5(r.h5); a.var_names_make_unique()
        pos = pd.read_csv(r.positions); pos = pos[pos.iloc[:, 1] == 1]
        a = a[a.obs_names.isin(pos.iloc[:, 0])]
        mats.append(("kpmp", sec, sp.csr_matrix(a.X), pd.Index(a.var_names)))
    for coh, sec, M, gn in mats:
        Z, _ = vc.lognorm_z(M, gn)
        for k in SETS:
            gg = [g for g in sets[k] if g in Z.columns]
            if len(gg) < 4:
                continue
            rs = []
            for _ in range(200):
                perm = rng.permutation(gg); a, b = perm[: len(gg) // 2], perm[len(gg) // 2:]
                r = np.corrcoef(Z[a].mean(1), Z[b].mean(1))[0, 1]
                rs.append(2 * r / (1 + r))
            rows.append({"cohort": coh, "section": sec, "set": k, "n_genes": len(gg), "reliability_sb": float(np.nanmedian(rs))})
    out = stage_dir(STAGE)
    R = pd.DataFrame(rows); R.to_csv(out / "spot_score_reliability.tsv", sep="\t", index=False)
    print(R.groupby(["cohort", "set"]).reliability_sb.describe().round(3))
    write_provenance(STAGE, [proc / "GSE211785_ST_counts.mtx"], [out / "spot_score_reliability.tsv"], seed)


if __name__ == "__main__":
    main()
