"""06_validate_gse195460: locked programmes in an independent snRNA cohort (Control vs DKD).

1. Label transfer from the discovery SN reference (fine types); donor-grouped CV on reference.
2. iPT fraction of PT per donor (composition replication).
3. Donor pseudobulk of iPT cells -> programme scores (DKDiPT, DEonly, SharedInjury);
   DKD vs Control Hedges g, Welch / Mann-Whitney p, expression-matched random-set null.
Only the dD direction is testable here (no HKD arm)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from annotate import fit_reference, predict, donor_pseudobulk, zmatrix, programme_scores
from stats import compare, matched_random_null

STAGE = "06_validate_gse195460"


def load_reference(cfg):
    a = ad.read_h5ad(ROOT / cfg["paths"]["processed"] / "GSE211785_rna_counts.h5ad")
    q = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    keep = (a.obs.library.map(q.batch) == "SN") & a.obs.library.map(q.qc_pass).astype(bool) & (a.obs.lineage != "Other")
    return a[keep.values].copy()


def gene_sets(cfg):
    gs = pd.read_csv(ROOT / "results/05_lock_programme/gene_sets.tsv", sep="\t")
    return {k: v.gene.tolist() for k, v in gs.groupby("set_id")}


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed)
    out = stage_dir(STAGE)
    qf = ROOT / cfg["paths"]["processed"] / "GSE195460_qc_counts.h5ad"
    q = ad.read_h5ad(qf)
    ref = load_reference(cfg)
    model = fit_reference(ref, q.var_names, "celltype", seed=seed)
    cv = model["cv"]
    cv_acc = cv.assign(ok=cv.true == cv.pred).groupby("true").ok.mean()
    cv_acc.to_csv(out / "reference_cv_accuracy.tsv", sep="\t")
    lab = predict(model, q)
    q.obs = q.obs.join(lab)
    lab.join(q.obs[["donor", "diagnosis"]]).to_csv(out / "cell_labels.tsv.gz", sep="\t")
    # composition
    pt = q.obs[q.obs.label.isin(cfg["lineage_map"]["PT"])]
    frac = pt.groupby("donor").apply(lambda d: (d.label == "iPT").mean()).rename("iPT_frac")
    dx = q.obs.groupby("donor").diagnosis.first()
    comp = pd.concat([frac, dx, pt.groupby("donor").size().rename("n_PT")], axis=1)
    comp.to_csv(out / "iPT_fraction_by_donor.tsv", sep="\t")
    res = {"iPT_fraction_DKD_vs_Control": compare(comp.iPT_frac, comp.diagnosis, "DKD", "Control")}
    # programme scores in iPT pseudobulk (and PT-all as secondary)
    sets = gene_sets(cfg)
    rows, cov_all = [], {}
    for comp_name, labels in [("iPT", ["iPT"]), ("PT_all", cfg["lineage_map"]["PT"])]:
        m = q.obs.label.isin(labels).values
        n = q.obs[m].groupby("donor").size()
        pb = donor_pseudobulk(q, m)
        pb = pb.loc[n[n >= cfg["pseudobulk"]["min_cells_per_sample"]].index]
        sc_, cov = programme_scores(pb, sets); cov_all[comp_name] = cov
        z, mean_expr = zmatrix(pb)
        grp = dx.reindex(sc_.index)
        sc_.join(grp).to_csv(out / f"scores_{comp_name}.tsv", sep="\t")
        for k in sets:
            r = compare(sc_[k], grp, "DKD", "Control")
            r.update(matched_random_null(z, sets[k], grp, "DKD", "Control", mean_expr, seed=seed))
            rows.append({"compartment": comp_name, "set": k, **cov[k], **r})
    tab = pd.DataFrame(rows); tab.to_csv(out / "programme_tests.tsv", sep="\t", index=False)
    res["coverage"] = cov_all
    (out / "summary.json").write_text(json.dumps(res, indent=2, default=float))
    print(cv_acc.round(3).to_string()); print(comp); print(tab.round(4).to_string())
    write_provenance(STAGE, [qf, ROOT / "results/05_lock_programme/gene_sets.tsv"],
                     sorted(out.glob("*.tsv")) + [out / "summary.json"], seed)


if __name__ == "__main__":
    main()
