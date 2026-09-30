"""02_discovery_pseudobulk: GSE211785 RNA cells -> raw-count object, library x cell-type
pseudobulk (fine types and lineages) and composition tables.

Outputs
  data/processed/GSE211785_rna_counts.h5ad   raw counts, RNA libraries of Control/HKD/DKD
  data/processed/pb_fine.h5ad, pb_lineage.h5ad   pseudobulk (obs = library x cell type)
  results/02_discovery_pseudobulk/library_table.tsv, composition_fine.tsv, pb_samples_*.tsv
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "02_discovery_pseudobulk"


def pseudobulk(X: sp.csr_matrix, obs: pd.DataFrame, key: str, var: pd.DataFrame) -> ad.AnnData:
    grp = obs["library"].astype(str) + "|" + obs[key].astype(str)
    codes, uniq = pd.factorize(grp)
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    M = (ind @ X).tocsr()
    meta = (obs.assign(_g=grp.values).groupby("_g", sort=False)
               .agg(library=("library", "first"), donor=("donor", "first"), tech=("tech", "first"),
                    group=("group", "first"), sex=("sex", "first"), age=("age", "first"),
                    celltype=(key, "first"), n_cells=("library", "size")).reindex(uniq))
    meta.index = uniq
    return ad.AnnData(X=M.astype(np.float32), obs=meta, var=var.copy())


def main():
    cfg = load_config(); set_global_seed(cfg["seed"])
    out = stage_dir(STAGE)
    proc = ROOT / cfg["paths"]["processed"]
    src = proc / "GSE211785_PreSCVI.h5ad"
    a = ad.read_h5ad(src)
    keep = a.obs.tech.isin(["SC_RNA", "SN_RNA"]) & a.obs.group.isin(cfg["cohort"]["groups"])
    a = a[keep.values].copy()
    X = a.layers["counts"].tocsr()
    assert np.allclose(X.data[:100000], np.round(X.data[:100000])), "counts layer is not integer"
    obs = pd.DataFrame({
        "library": a.obs.orig_ident.astype(str).values,
        "donor": a.obs["sample"].astype(str).str.extract(r"(HK\d+)")[0].values,
        "tech": a.obs.tech.astype(str).values, "group": a.obs.group.astype(str).values,
        "sex": a.obs.sex.astype(str).values, "age": a.obs.age.astype(float).values,
        "celltype": a.obs.Cluster_Idents.astype(str).values}, index=a.obs_names)
    lmap = {ct: lin for lin, cts in cfg["lineage_map"].items() for ct in cts}
    obs["lineage"] = obs.celltype.map(lmap).fillna("Other")
    var = pd.DataFrame(index=a.var_names.astype(str))
    rna = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var)
    rna.write_h5ad(proc / "GSE211785_rna_counts.h5ad", compression="gzip")

    pb_f = pseudobulk(X, obs, "celltype", var); pb_f.obs["lineage"] = pb_f.obs.celltype.map(lmap).fillna("Other")
    pb_l = pseudobulk(X, obs, "lineage", var)
    pb_f.write_h5ad(proc / "pb_fine.h5ad", compression="gzip")
    pb_l.write_h5ad(proc / "pb_lineage.h5ad", compression="gzip")

    lib = (obs.groupby("library").agg(donor=("donor", "first"), tech=("tech", "first"), group=("group", "first"),
                                     sex=("sex", "first"), age=("age", "first"), n_cells=("donor", "size"))
              .reset_index())
    lib.to_csv(out / "library_table.tsv", sep="\t", index=False)
    comp = pd.crosstab(obs.library, obs.celltype)
    comp.to_csv(out / "composition_fine.tsv", sep="\t")
    for name, pb in [("fine", pb_f), ("lineage", pb_l)]:
        pb.obs.to_csv(out / f"pb_samples_{name}.tsv", sep="\t")
    summ = (pb_l.obs[pb_l.obs.n_cells >= cfg["pseudobulk"]["min_cells_per_sample"]]
              .groupby(["celltype", "group"]).donor.nunique().unstack(fill_value=0))
    summ.to_csv(out / "donors_per_lineage_group.tsv", sep="\t")
    print(lib.groupby(["group", "tech"]).size()); print(summ)
    write_provenance(STAGE, [src], [proc / "GSE211785_rna_counts.h5ad", proc / "pb_fine.h5ad",
                     proc / "pb_lineage.h5ad", *sorted(out.glob("*.tsv"))], cfg["seed"],
                     extra={"n_cells": int(rna.n_obs), "n_genes": int(rna.n_vars)})


if __name__ == "__main__":
    main()
