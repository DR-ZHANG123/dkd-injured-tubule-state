"""04_scdisinfact: disentangled representation learning in the focus lineage (plan rule 4).

Cells of the focus lineage from the chosen libraries (primary: SN libraries; sensitivity: all
QC-passing libraries) are fitted with scDisInFact, condition = group (Control/HKD/DKD),
batch = library. Five seeds; per gene the condition-associated gene (CKG) score is averaged
over seeds and ranked. Outputs results/04_scdisinfact/<run>/ckg_scores.tsv, latent_*.tsv,
losses.tsv, and a per-seed rank-correlation table.
Usage: python 04_scdisinfact.py [primary|all_libraries]
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "04_scdisinfact"


def load_cells(cfg: dict, run: str) -> ad.AnnData:
    p = cfg["scdisinfact"]
    a = ad.read_h5ad(ROOT / cfg["paths"]["processed"] / "GSE211785_rna_counts.h5ad")
    q = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    a.obs["batch"] = a.obs.library.map(q.batch).astype(str)
    keep = a.obs.celltype.isin(p["focus_celltypes"]) & a.obs.library.map(q.qc_pass).astype(bool)
    if run == "primary":
        keep &= a.obs.batch == "SN"
    a = a[keep.values].copy()
    # libraries with too few focus cells cannot inform the batch-conditional MMD terms
    n = a.obs.library.value_counts()
    a = a[a.obs.library.isin(n[n >= p["min_cells_library"]].index).values].copy()
    sc.pp.filter_genes(a, min_cells=p["min_cells_gene"])
    # batch-stratified loess fits are singular for sparse strata; HVGs on pooled cells
    sc.pp.highly_variable_genes(a, n_top_genes=p["n_hvg"], flavor="seurat_v3", subset=False)
    return a[:, a.var.highly_variable.values].copy()


def fit_one(a: ad.AnnData, seed: int, cfg: dict):
    from scDisInFact import scdisinfact, create_scdisinfact_dataset
    p = cfg["scdisinfact"]
    set_global_seed(seed)
    counts = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
    meta = a.obs[["group", "library"]].copy()
    data_dict = create_scdisinfact_dataset(counts, meta, condition_key=["group"], batch_key="library")
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = scdisinfact(data_dict=data_dict, Ks=p["Ks"], batch_size=p["batch_size"], interval=p["interval"],
                        lr=p["lr"], reg_mmd_comm=p["reg_mmd_comm"], reg_mmd_diff=p["reg_mmd_diff"],
                        reg_gl=p["reg_gl"], reg_class=p["reg_class"], reg_kl_comm=p["reg_kl_comm"],
                        reg_kl_diff=p["reg_kl_diff"], seed=seed, device=dev)
    model.train()
    losses = model.train_model(nepochs=p["nepochs"], recon_loss="NB")
    model.eval()
    score = model.extract_gene_scores()[0]
    zc, zd, idx = [], [], []
    with torch.no_grad():
        for ds, mc in zip(data_dict["datasets"], data_dict["meta_cells"]):
            inf = model.inference(counts=ds.counts_norm.to(dev), batch_ids=ds.batch_id[:, None].to(dev), print_stat=False)
            zc.append(inf["mu_c"].cpu().numpy()); zd.append(inf["mu_d"][0].cpu().numpy()); idx += list(mc.index)
    return score, np.vstack(zc), np.vstack(zd), idx, losses


def main():
    cfg = load_config()
    run = sys.argv[1] if len(sys.argv) > 1 else "primary"
    out = stage_dir(STAGE) / run; out.mkdir(exist_ok=True)
    a = load_cells(cfg, run)
    print(run, a.shape, a.obs.groupby("group").library.nunique().to_dict(), flush=True)
    scores, loss_rows = {}, []
    for seed in cfg["scdisinfact"]["seeds"]:
        t = time.time()
        s, zc, zd, idx, losses = fit_one(a, seed, cfg)
        scores[f"seed{seed}"] = s
        lat = pd.DataFrame(np.hstack([zc, zd]), index=idx,
                           columns=[f"zc{i}" for i in range(zc.shape[1])] + [f"zd{i}" for i in range(zd.shape[1])])
        lat = lat.join(a.obs[["library", "donor", "group", "celltype"]])
        lat.to_csv(out / f"latent_seed{seed}.tsv.gz", sep="\t")
        loss_rows.append({"seed": seed, "minutes": (time.time() - t) / 60,
                          **{f"final_{k}": float(np.ravel(v)[-1]) for k, v in zip(
                              ["total", "recon", "kl_comm", "kl_diff", "mmd_comm", "mmd_diff", "class", "gl"], losses)
                             if hasattr(v, "__len__") and len(np.ravel(v))}})
        print(loss_rows[-1], flush=True)
    S = pd.DataFrame(scores, index=a.var_names)
    ranks = S.rank(ascending=False)
    S["mean_rank"] = ranks.mean(1); S["mean_score"] = S[list(scores)].mean(1)
    S["pct_rank"] = S.mean_rank.rank() / len(S)
    S = S.sort_values("mean_rank"); S.index.name = "gene"
    S.to_csv(out / "ckg_scores.tsv", sep="\t")
    ranks.corr(method="spearman").to_csv(out / "seed_rank_spearman.tsv", sep="\t")
    pd.DataFrame(loss_rows).to_csv(out / "losses.tsv", sep="\t", index=False)
    print(S.head(30)); print(ranks.corr(method="spearman").round(3))
    write_provenance(f"{STAGE}", [ROOT / cfg["paths"]["processed"] / "GSE211785_rna_counts.h5ad",
                                  ROOT / "results/02c_library_qc/library_qc.tsv"],
                     sorted(out.glob("*")), cfg["scdisinfact"]["seeds"][0],
                     extra={"run": run, "n_cells": int(a.n_obs), "n_genes": int(a.n_vars),
                            "params": cfg["scdisinfact"]})
    (out / "PROVENANCE.json").write_text((stage_dir(STAGE) / "PROVENANCE.json").read_text())


if __name__ == "__main__":
    main()
