"""03b_technical_axis: does a DKD-HKD effect estimate track the nuclear/cytoplasmic axis?

Technical axis per lineage = mean log2 CPM of Control libraries in SN minus SC_A. For each
lineage and DE mode the Pearson r between dE log2FC and the axis is reported (all genes, and
genes with dE FDR < 0.10)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "03b_technical_axis"


def main():
    cfg = load_config(); set_global_seed(cfg["seed"])
    src = ROOT / cfg["paths"]["processed"] / "pb_lineage.h5ad"
    q = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    pb = ad.read_h5ad(src)
    pb.obs["batch"] = pb.obs.library.map(q.batch)
    ok = pb.obs.library.map(q.qc_pass).astype(bool) & (pb.obs.n_cells >= cfg["pseudobulk"]["min_cells_per_sample"])
    X = pb.X.toarray(); cpm = np.log2(X / X.sum(1, keepdims=True) * 1e6 + 1)
    rows, axes = [], {}
    de = ROOT / "results/03_pseudobulk_de"
    for lin in sorted(pb.obs.celltype.unique()):
        o = pb.obs
        m = ok & (o.celltype == lin) & (o.group == "Control")
        a, b = (m & (o.batch == "SN")).values, (m & (o.batch == "SC_A")).values
        if a.sum() < 2 or b.sum() < 2:
            continue
        tax = pd.Series(cpm[a].mean(0) - cpm[b].mean(0), index=pb.var_names)
        axes[lin] = tax
        for mode in ["all_libraries", "one_lib_per_donor", "SC_A_only", "SN_only"]:
            f = de / f"lineage_{mode}" / f"{lin}.tsv"
            if not f.exists():
                continue
            r = pd.read_csv(f, sep="\t").set_index("gene")
            g = r.index.intersection(tax.index); s = r.loc[g]; sig = (s.dE_padj < 0.10).values
            rows.append({"lineage": lin, "mode": mode, "n_genes": len(g), "n_sig": int(sig.sum()),
                         "r_all": np.corrcoef(s.dE_lfc, tax[g])[0, 1],
                         "r_sig": np.corrcoef(s.dE_lfc[sig], tax[g][sig])[0, 1] if sig.sum() > 5 else np.nan})
    out = stage_dir(STAGE)
    res = pd.DataFrame(rows); res.to_csv(out / "dE_vs_technical_axis.tsv", sep="\t", index=False)
    pd.DataFrame(axes).to_csv(out / "technical_axis_by_lineage.tsv", sep="\t")
    print(res.round(3).to_string())
    write_provenance(STAGE, [src, ROOT / "results/02c_library_qc/library_qc.tsv"],
                     [out / "dE_vs_technical_axis.tsv", out / "technical_axis_by_lineage.tsv"], cfg["seed"])


if __name__ == "__main__":
    main()
