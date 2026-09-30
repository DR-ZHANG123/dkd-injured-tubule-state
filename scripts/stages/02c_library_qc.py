"""02c_library_qc: per-library technical profile and protocol batch of the discovery atlas.

scRNA libraries fall into two protocol batches separable by the mitochondrial / ribosomal read
fraction of the library pseudobulk (SC_A: mt <= threshold; SC_B: mt > threshold). A scRNA
library with a ribosomal fraction below ribo_fail carries a nuclear-like profile and is excluded.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "02c_library_qc"


def main():
    cfg = load_config(); set_global_seed(cfg["seed"]); q = cfg["library_qc"]
    src = ROOT / cfg["paths"]["processed"] / "pb_lineage.h5ad"
    pb = ad.read_h5ad(src); g = pb.var_names
    X = pb.X.toarray() if hasattr(pb.X, "toarray") else np.asarray(pb.X)
    o = pb.obs.copy()
    o["tot"] = X.sum(1)
    o["mt"] = X[:, g.str.startswith("MT-")].sum(1)
    o["ribo"] = X[:, g.str.match(r"^RP[SL]\d")].sum(1)
    o["malat1"] = X[:, g.isin(["MALAT1"])].sum(1)
    lib = o.groupby("library", observed=True).agg(
        donor=("donor", "first"), tech=("tech", "first"), group=("group", "first"),
        sex=("sex", "first"), age=("age", "first"), n_cells=("n_cells", "sum"),
        tot=("tot", "sum"), mt=("mt", "sum"), ribo=("ribo", "sum"), malat1=("malat1", "sum"))
    for c in ["mt", "ribo", "malat1"]:
        lib[f"frac_{c}"] = lib[c] / lib.tot
    sc = lib.tech == "SC_RNA"
    lib["batch"] = np.where(~sc, "SN", np.where(lib.frac_mt > q["sc_mt_split"], "SC_B", "SC_A"))
    lib["qc_pass"] = ~(sc & (lib.frac_ribo < q["sc_ribo_fail"])) & (lib.n_cells >= q["min_cells_library"])
    lib = lib.drop(columns=["mt", "ribo", "malat1"]).reset_index()
    out = stage_dir(STAGE)
    lib.to_csv(out / "library_qc.tsv", sep="\t", index=False)
    tab = pd.crosstab([lib.batch, lib.qc_pass], lib.group)
    tab.to_csv(out / "batch_by_group.tsv", sep="\t")
    print(lib.round(4).to_string()); print(tab)
    write_provenance(STAGE, [src], [out / "library_qc.tsv", out / "batch_by_group.tsv"], cfg["seed"])


if __name__ == "__main__":
    main()
