"""02b_prepare_gse195460: read the 11 snRNA libraries of the replication cohort, apply QC,
attach donor clinical fields and save raw counts. No annotation here (stage 05 transfers labels
from the discovery reference)."""
from __future__ import annotations
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "02b_prepare_gse195460"


def main():
    cfg = load_config(); set_global_seed(cfg["seed"])
    q = cfg["replication_qc"]
    raw = ROOT / cfg["paths"]["raw"] / "GSE195460"
    out = stage_dir(STAGE)
    clin = pd.read_csv(ROOT / "results/01_cohort_audit/sample_table_GSE195460.tsv", sep="\t")
    files = sorted(raw.glob("*_filtered_feature_bc_matrix.h5"))
    adatas, rows = [], []
    for f in files:
        lib = re.search(r"_(Control\d|DN\d)_", f.name).group(1)
        a = sc.read_10x_h5(f); a.var_names_make_unique()
        a.obs_names = [f"{lib}_{b}" for b in a.obs_names]
        a.var["mt"] = a.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        n0 = a.n_obs
        keep = ((a.obs.n_genes_by_counts >= q["min_genes"]) & (a.obs.n_genes_by_counts <= q["max_genes"])
                & (a.obs.pct_counts_mt <= q["max_pct_mt"]))
        a = a[keep].copy()
        a.obs["library"] = lib
        adatas.append(a)
        rows.append({"library": lib, "cells_raw": n0, "cells_qc": a.n_obs,
                     "median_genes": float(np.median(a.obs.n_genes_by_counts))})
    genes = sorted(set.intersection(*[set(a.var_names) for a in adatas]))
    adata = ad.concat([a[:, genes] for a in adatas], join="inner")
    c = clin.set_index("library")
    adata.obs["diagnosis"] = adata.obs.library.map(c.diagnosis)
    for col in ["egfr", "ifta", "global glomerulosclerosis", "age", "sex"]:
        if col in c.columns:
            adata.obs[col.replace(" ", "_")] = adata.obs.library.map(c[col]).astype(str)
    adata.obs["donor"] = adata.obs.library
    dst = ROOT / cfg["paths"]["processed"] / "GSE195460_qc_counts.h5ad"
    adata.write_h5ad(dst, compression="gzip")
    qc = pd.DataFrame(rows).merge(c.reset_index()[["library", "diagnosis"] +
                                  [x for x in ["egfr", "ifta", "global glomerulosclerosis", "age", "sex"] if x in c.columns]],
                                  on="library", how="left")
    qc.to_csv(out / "library_qc.tsv", sep="\t", index=False)
    print(qc.to_string())
    write_provenance(STAGE, files, [dst, out / "library_qc.tsv"], cfg["seed"],
                     extra={"n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars)})


if __name__ == "__main__":
    main()
