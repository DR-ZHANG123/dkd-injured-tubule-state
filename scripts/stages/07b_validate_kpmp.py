"""07b_validate_kpmp: primary KPMP test of the marker-gene programmes.

Per modality (sn, sc): label transfer from the matching discovery reference (SN or SC libraries),
participant-level iPT pseudobulk, programme scores, group contrasts, expression-matched random
null, and severity-adjusted models (eGFR band, cortical interstitial fibrosis %)."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from annotate import fit_reference, predict, donor_pseudobulk, zmatrix
from stats import compare, matched_random_null, hedges_g
from kpmp import read_descriptors

STAGE = "07b_validate_kpmp"
CONTRASTS = [("DKD", "HKD"), ("DKD", "HKD_withDM"), ("DKD", "DM_noCKD"), ("HKD_withDM", "HKD"),
             ("DKD", "Reference"), ("HKD", "Reference")]


def reference(cfg, tag):
    a = ad.read_h5ad(ROOT / cfg["paths"]["processed"] / "GSE211785_rna_counts.h5ad")
    q = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    b = a.obs.library.map(q.batch)
    keep = a.obs.library.map(q.qc_pass).astype(bool) & (a.obs.lineage != "Other")
    keep &= (b == "SN") if tag == "sn" else b.isin(["SC_A", "SC_B"])
    return a[keep.values].copy()


def band_lower(x: str) -> float:
    try:
        return float(str(x).replace("<", "").replace(">", "").split("-")[0].split()[0])
    except ValueError:
        return np.nan


def severity_models(scores: pd.DataFrame, meta: pd.DataFrame, sets) -> list[dict]:
    rows = []
    d = scores.join(meta)
    d = d[d.group.isin(["DKD", "HKD"])].copy()
    d["dkd"] = (d.group == "DKD").astype(int)
    for k in sets:
        d["y"] = d[k]
        for name, form, need in [("unadjusted", "y ~ dkd", []), ("eGFR_band", "y ~ dkd + egfr_lower", ["egfr_lower"]),
                                 ("IF_percent", "y ~ dkd + if_pct", ["if_pct"]),
                                 ("eGFR_and_IF", "y ~ dkd + egfr_lower + if_pct", ["egfr_lower", "if_pct"])]:
            dd = d.dropna(subset=["y", *need])
            if dd.dkd.nunique() < 2 or len(dd) < 8:
                continue
            m = smf.ols(form, dd).fit(cov_type="HC3")
            rows.append({"set": k, "model": name, "n": len(dd), "n_DKD": int(dd.dkd.sum()),
                         "beta_DKD": m.params["dkd"], "se": m.bse["dkd"], "p": m.pvalues["dkd"]})
    return rows


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed)
    out = stage_dir(STAGE)
    gs = pd.read_csv(ROOT / "results/05_lock_programme/gene_sets.tsv", sep="\t")
    sets = {k: v.gene.tolist() for k, v in gs.groupby("set_id")}
    desc = read_descriptors(ROOT / "data/raw/kpmp" / cfg["kpmp"]["descriptor_table"]).set_index("Participant ID")
    if_col = "Cortex: Interstitial Fibrosis: Percent"
    tests, sev, meta_all, cv_all = [], [], {}, {}
    for tag in ["sn", "sc"]:
        f = ROOT / cfg["paths"]["processed"] / f"KPMP_{tag}_qc_counts.h5ad"
        q = ad.read_h5ad(f)
        model = fit_reference(reference(cfg, tag), q.var_names, "celltype", seed=seed)
        cv = model["cv"]; cv_all[tag] = cv.assign(ok=cv.true == cv.pred).groupby("true").ok.mean()
        lab = predict(model, q); q.obs = q.obs.join(lab)
        q.obs[["participant", "group", "label", "max_prob"]].to_csv(out / f"cell_labels_{tag}.tsv.gz", sep="\t")
        meta = q.obs.groupby("participant").agg(group=("group", "first"), egfr_band=("egfr_band", "first"),
                                                 n_cells=("group", "size"))
        pt = q.obs[q.obs.label.isin(cfg["lineage_map"]["PT"])]
        meta["n_PT"] = pt.groupby("participant").size()
        meta["iPT_frac"] = pt.groupby("participant").apply(lambda d: (d.label == "iPT").mean())
        meta["egfr_lower"] = meta.egfr_band.map(band_lower)
        # one participant is listed twice with identical scores; collapse per participant
        meta["if_pct"] = pd.to_numeric(desc[if_col], errors="coerce").groupby(level=0).mean().reindex(meta.index)
        m = q.obs.label.eq("iPT").values
        n = q.obs[m].groupby("participant").size()
        meta["n_iPT"] = n
        pb = donor_pseudobulk(q, m, "participant")
        pb = pb.loc[n[n >= cfg["pseudobulk"]["min_cells_per_sample"]].index]
        z, mean_expr = zmatrix(pb)
        sc_ = pd.DataFrame({k: z[[g for g in v if g in z.columns]].mean(1) for k, v in sets.items()})
        cov = {k: len([g for g in v if g in z.columns]) for k, v in sets.items()}
        sc_.join(meta).to_csv(out / f"scores_{tag}.tsv", sep="\t")
        meta.to_csv(out / f"participants_{tag}.tsv", sep="\t"); meta_all[tag] = meta
        grp = meta.group.reindex(sc_.index)
        for a, b in CONTRASTS:
            tests.append({"modality": tag, "set": "iPT_fraction", "contrast": f"{a}-{b}",
                          **compare(meta.iPT_frac, meta.group, a, b)})
            for k in sets:
                r = compare(sc_[k], grp, a, b)
                if (a, b) in [("DKD", "HKD"), ("DKD", "HKD_withDM"), ("DKD", "DM_noCKD")] and r.get("n_a", 0) >= 3 and r.get("n_b", 0) >= 3:
                    r.update(matched_random_null(z, sets[k], grp, a, b, mean_expr, seed=seed))
                tests.append({"modality": tag, "set": k, "contrast": f"{a}-{b}", "n_measured": cov[k], **r})
        sev += [dict(r, modality=tag) for r in severity_models(sc_, meta, sets)]
        print(tag, "iPT pseudobulk participants", pb.shape[0], grp.value_counts().to_dict(), flush=True)
    T = pd.DataFrame(tests); T.to_csv(out / "programme_tests.tsv", sep="\t", index=False)
    S = pd.DataFrame(sev); S.to_csv(out / "severity_models.tsv", sep="\t", index=False)
    # random-effects (DerSimonian-Laird) pooling of Hedges g across modalities
    pooled = []
    for (k, c), g in T[T.hedges_g.notna()].groupby(["set", "contrast"]):
        if len(g) < 2:
            continue
        na, nb, gg = g.n_a.values, g.n_b.values, g.hedges_g.values
        v = (na + nb) / (na * nb) + gg ** 2 / (2 * (na + nb)); w = 1 / v
        fe = np.sum(w * gg) / w.sum(); Q = np.sum(w * (gg - fe) ** 2)
        tau2 = max(0, (Q - (len(g) - 1)) / (w.sum() - np.sum(w ** 2) / w.sum()))
        wr = 1 / (v + tau2); re = np.sum(wr * gg) / wr.sum(); se = np.sqrt(1 / wr.sum())
        from scipy.stats import norm
        pooled.append({"set": k, "contrast": c, "g_pooled": re, "se": se, "ci_lo": re - 1.96 * se,
                       "ci_hi": re + 1.96 * se, "p": 2 * norm.sf(abs(re / se)), "tau2": tau2})
    P = pd.DataFrame(pooled); P.to_csv(out / "pooled_effects.tsv", sep="\t", index=False)
    pd.DataFrame(cv_all).to_csv(out / "reference_cv_accuracy.tsv", sep="\t")
    ov = set(meta_all["sn"].index) & set(meta_all["sc"].index)
    (out / "summary.json").write_text(json.dumps({"participants_in_both_modalities": len(ov)}, indent=2))
    pd.set_option("display.width", 250)
    print(T[T.set.isin(["DKDiPT_up", "DEonly_up", "SharedInjury_up", "iPT_fraction"])][
        ["modality", "set", "contrast", "n_a", "n_b", "hedges_g", "welch_p", "p_emp_greater"]].round(4).to_string())
    print(P.round(4).to_string()); print(S.round(4).to_string())
    write_provenance(STAGE, [ROOT / cfg["paths"]["processed"] / "KPMP_sn_qc_counts.h5ad",
                             ROOT / cfg["paths"]["processed"] / "KPMP_sc_qc_counts.h5ad",
                             ROOT / "results/05_lock_programme/gene_sets.tsv"],
                     sorted(out.glob("*.tsv")) + [out / "summary.json"], seed)


if __name__ == "__main__":
    main()
