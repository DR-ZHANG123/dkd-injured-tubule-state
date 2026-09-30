"""11b_histology_association: same-section H&E morphology vs measured programme scores
(analysis_plan rule 9). Spots are not independent: inference is across sections.

A. Interpretable morphology: per section, partial Spearman correlation of each morphology feature
   with the measured DKDiPT_up spot score, adjusting for SharedInjury_up, prop_iPT and PT_share
   (ranks residualised on covariates). Across sections: mean Fisher z, one-sample t-test (df =
   sections - 1), sign count. Same for SharedInjury_up (adjusting DKDiPT_up) as comparison.
   A spot-level linear mixed model (random intercept per section) is reported alongside.
B. Frozen-encoder prediction (only if >= min_sections_prediction sections): ridge head on phikon-v2
   CLS embeddings predicting the within-section z of each measured score, leave-one-participant-out;
   per held-out participant Pearson r; baselines: morphology-feature ridge, and the same ridge on the
   SharedInjury_up target. Increment beyond injury/composition: partial r of prediction with the
   measured DKDiPT_up given measured SharedInjury_up + prop_iPT. Predictions are not used as
   evidence for the programme itself.
Usage: python 11b_histology_association.py [gse211785|kpmp|all]
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "11b_histology_association"
EMB = ROOT / "data/processed/histology_embeddings"
FEATS = ["hema_frac", "nuclei_per_1e4px_tissue", "eosin_frac", "eosin_mean", "hema_mean", "white_frac"]
COVS = ["prop_iPT", "PT_share"]


def load(cohort: str) -> pd.DataFrame:
    f = pd.read_csv(ROOT / f"results/11a_histology_embed/features_{cohort}.tsv.gz", sep="\t")
    if cohort == "gse211785":
        s = pd.read_csv(ROOT / "results/10a_visium_gse211785/spot_scores.tsv.gz", sep="\t")
        s = s.rename(columns={"donor": "participant"})
        m = f.merge(s, left_on=["section", "spot_id"], right_on=["section", "barcode"])
    else:
        s = pd.read_csv(ROOT / "results/10b_visium_kpmp/spot_scores.tsv.gz", sep="\t")
        s["section"] = s.participant.astype(str) + "__" + s["sample"].astype(str)
        m = f.merge(s, left_on=["section", "spot_id"], right_on=["section", "barcode"])
    return m[m.tissue_pass].copy()


def resid_rank(y, X):
    y = stats.rankdata(y); X = np.column_stack([np.ones(len(y))] + [stats.rankdata(x) for x in X.T])
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def morph_assoc(d: pd.DataFrame, cfg) -> tuple[pd.DataFrame, pd.DataFrame]:
    per, rows = [], []
    for target, adj in [("DKDiPT_up", ["SharedInjury_up"] + COVS), ("SharedInjury_up", ["DKDiPT_up"] + COVS)]:
        for sec, g in d.groupby("section"):
            g = g.dropna(subset=[target, *adj, *FEATS])
            if len(g) < 50:
                continue
            ry = resid_rank(g[target].values, g[adj].values)
            for ft in FEATS:
                rf = resid_rank(g[ft].values, g[adj].values)
                per.append({"target": target, "section": sec, "participant": g.participant.iloc[0], "feature": ft,
                            "n_spots": len(g), "partial_rho": stats.pearsonr(ry, rf)[0]})
    P = pd.DataFrame(per)
    for (t, ft), g in P.groupby(["target", "feature"]):
        z = np.arctanh(g.partial_rho.clip(-0.999, 0.999))
        tt = stats.ttest_1samp(z, 0)
        rows.append({"target": t, "feature": ft, "n_sections": len(g), "mean_partial_rho": float(np.tanh(z.mean())),
                     "t": tt.statistic, "p_sections": tt.pvalue, "n_positive": int((g.partial_rho > 0).sum())})
    S = pd.DataFrame(rows)
    S["fdr_sections"] = np.nan
    for t in S.target.unique():
        i = S.target == t
        S.loc[i, "fdr_sections"] = stats.false_discovery_control(S.loc[i, "p_sections"].values)
    return P, S


def mixed(d: pd.DataFrame) -> pd.DataFrame:
    import statsmodels.formula.api as smf
    rows = []
    dd = d.dropna(subset=["DKDiPT_up", "SharedInjury_up", *COVS, *FEATS]).copy()
    for c in ["DKDiPT_up", "SharedInjury_up", *COVS, *FEATS]:
        dd[c] = (dd[c] - dd[c].mean()) / dd[c].std()
    for ft in FEATS:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = smf.mixedlm(f"{ft} ~ DKDiPT_up + SharedInjury_up + prop_iPT + PT_share", dd, groups=dd.section).fit(reml=True)
        for term in ["DKDiPT_up", "SharedInjury_up"]:
            rows.append({"feature": ft, "term": term, "beta_sd": m.params[term], "se": m.bse[term],
                         "p_spot_level": m.pvalues[term], "n_spots": len(dd), "n_sections": dd.section.nunique()})
    return pd.DataFrame(rows)


def section_level(d: pd.DataFrame) -> pd.DataFrame:
    """Section means (tissue spots) of morphology and measured scores; OLS across sections:
    feature ~ DKDiPT_up + SharedInjury_up + prop_iPT (HC3 SE). One section per participant is used
    (the one with most spots) so that n = participants."""
    import statsmodels.formula.api as smf
    agg = d.groupby("section").agg(participant=("participant", "first"), n=("spot_id", "size"),
                                   **{c: (c, "mean") for c in FEATS + ["DKDiPT_up", "SharedInjury_up", "prop_iPT"]})
    agg = agg.sort_values("n", ascending=False).drop_duplicates("participant")
    z = agg.copy()
    for c in FEATS + ["DKDiPT_up", "SharedInjury_up", "prop_iPT"]:
        z[c] = (z[c] - z[c].mean()) / z[c].std()
    rows = []
    for ft in FEATS:
        m = smf.ols(f"{ft} ~ DKDiPT_up + SharedInjury_up + prop_iPT", z).fit(cov_type="HC3")
        for term in ["DKDiPT_up", "SharedInjury_up"]:
            rows.append({"feature": ft, "term": term, "beta_sd": m.params[term], "se": m.bse[term],
                         "p": m.pvalues[term], "n_participants": len(z)})
    return pd.DataFrame(rows)


def within_z(d, col):
    """Within-section z; a covariate constant within a section (e.g. no iPT signal) becomes 0."""
    z = d.groupby("section")[col].transform(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x * 0.0)
    return z.fillna(0.0) if col == "prop_iPT" else z


def predict(d: pd.DataFrame, cfg, seed) -> tuple[pd.DataFrame, pd.DataFrame]:
    from sklearn.linear_model import RidgeCV
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import make_pipeline
    emb = {}
    for sec in d.section.unique():
        z = np.load(EMB / f"{d.cohort_name.iloc[0]}__{sec}.npz")
        # Visium barcodes repeat across slides: key by (section, barcode)
        emb.update({(sec, sid): e for sid, e in zip(z["spot_id"], z["emb"].astype(np.float32))})
    d = d[[(a, b) in emb for a, b in zip(d.section, d.spot_id)]].copy()
    for c in ["DKDiPT_up", "SharedInjury_up", "prop_iPT"]:
        d[f"{c}_wz"] = within_z(d, c)
    d = d.dropna(subset=["DKDiPT_up_wz", "SharedInjury_up_wz", "prop_iPT_wz"])
    alphas = cfg["histology"]["ridge_alphas"]
    rows, preds = [], []
    parts = sorted(d.participant.unique())
    Eall = np.vstack([emb[(a, b)] for a, b in zip(d.section, d.spot_id)]); Fall = d[FEATS].fillna(d[FEATS].median()).values
    for p in parts:
        te = (d.participant == p).values; tr = ~te
        for target in ["DKDiPT_up_wz", "SharedInjury_up_wz"]:
            for name, X in [("encoder", Eall), ("morphology", Fall)]:
                mdl = make_pipeline(StandardScaler(), RidgeCV(alphas=alphas)).fit(X[tr], d[target].values[tr])
                yh = mdl.predict(X[te]); y = d[target].values[te]
                r = stats.pearsonr(yh, y)[0]
                row = {"held_out": p, "target": target.replace("_wz", ""), "features": name, "n_test_spots": int(te.sum()),
                       "pearson_r": r, "alpha": float(mdl[-1].alpha_)}
                if target == "DKDiPT_up_wz":
                    g = d[te]
                    A = np.c_[np.ones(te.sum()), g.SharedInjury_up_wz, g.prop_iPT_wz]
                    ry = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]; rp = yh - A @ np.linalg.lstsq(A, yh, rcond=None)[0]
                    row["partial_r_given_injury_iPT"] = stats.pearsonr(rp, ry)[0]
                rows.append(row)
                if name == "encoder":
                    preds.append(pd.DataFrame({"spot_id": d.spot_id.values[te], "section": d.section.values[te],
                                               "target": target, "pred": yh, "measured": y}))
    return pd.DataFrame(rows), pd.concat(preds)


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    out = stage_dir(STAGE)
    outs, ins = [], []
    for coh in (["gse211785", "kpmp"] if which == "all" else [which]):
        fi = ROOT / f"results/11a_histology_embed/features_{coh}.tsv.gz"
        if not fi.exists():
            continue
        d = load(coh).reset_index(drop=True); d["cohort_name"] = coh; ins.append(fi)
        P, S = morph_assoc(d, cfg); M = mixed(d); L = section_level(d)
        L.to_csv(out / f"morph_section_level_{coh}.tsv", sep="\t", index=False); outs.append(out / f"morph_section_level_{coh}.tsv")
        P.to_csv(out / f"morph_partial_by_section_{coh}.tsv", sep="\t", index=False)
        S.to_csv(out / f"morph_assoc_summary_{coh}.tsv", sep="\t", index=False)
        M.to_csv(out / f"morph_mixed_model_{coh}.tsv", sep="\t", index=False)
        outs += [out / f"morph_partial_by_section_{coh}.tsv", out / f"morph_assoc_summary_{coh}.tsv",
                 out / f"morph_mixed_model_{coh}.tsv"]
        pd.set_option("display.width", 250)
        print("==", coh, "sections", d.section.nunique(), "participants", d.participant.nunique(), "spots", len(d))
        print(S.round(4).to_string()); print(M.round(4).to_string()); print(L.round(4).to_string())
        if d.participant.nunique() >= cfg["histology"]["min_sections_prediction"]:
            R, Pr = predict(d, cfg, seed)
            R.to_csv(out / f"prediction_lopo_{coh}.tsv", sep="\t", index=False)
            Pr.to_csv(out / f"prediction_spots_{coh}.tsv.gz", sep="\t", index=False)
            summ = R.groupby(["target", "features"]).agg(n_folds=("pearson_r", "size"), mean_r=("pearson_r", "mean"),
                                                         median_r=("pearson_r", "median"),
                                                         min_r=("pearson_r", "min"), max_r=("pearson_r", "max"),
                                                         mean_partial_r=("partial_r_given_injury_iPT", "mean")).reset_index()
            summ["p_folds_r_gt0"] = [stats.ttest_1samp(np.arctanh(R[(R.target == t) & (R.features == f)].pearson_r), 0).pvalue
                                     for t, f in zip(summ.target, summ.features)]
            summ["p_folds_partial_gt0"] = [stats.ttest_1samp(np.arctanh(R[(R.target == t) & (R.features == f)]
                                           .partial_r_given_injury_iPT.dropna()), 0).pvalue if t == "DKDiPT_up" else np.nan
                                           for t, f in zip(summ.target, summ.features)]
            summ["n_folds_partial_gt0"] = [int((R[(R.target == t) & (R.features == f)].partial_r_given_injury_iPT > 0).sum())
                                           if t == "DKDiPT_up" else np.nan for t, f in zip(summ.target, summ.features)]
            summ.to_csv(out / f"prediction_summary_{coh}.tsv", sep="\t", index=False)
            outs += [out / f"prediction_lopo_{coh}.tsv", out / f"prediction_summary_{coh}.tsv"]
            print(summ.round(4).to_string())
    write_provenance(STAGE, ins + [ROOT / "results/10a_visium_gse211785/spot_scores.tsv.gz",
                                   ROOT / "results/10b_visium_kpmp/spot_scores.tsv.gz"], outs, seed,
                     extra={"cohorts": which})


if __name__ == "__main__":
    main()
