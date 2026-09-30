"""13b_wsi_association: participant-level links between whole-slide histology and programmes.

Unit = participant (one slide per stain). Slide representation = mean phikon-v2 embedding over
QC tiles ("all") or over the pre-specified tubule-rich subset ("tubule"); plus mean
interpretable features. Analyses (per stain x tile set):
  1. sanity: cross-validated ridge predicting TIV cortical IF% and tubular atrophy %
  2. programme scores (07b, per modality): CV ridge predicting DKDiPT_up and SharedInjury_up;
     partial target = DKDiPT_up residualised on SharedInjury_up and IF%
  3. context: CV logistic DKD vs HKD AUC (not a diagnostic model)
CV r / AUC get permutation p-values (y permuted, whole CV repeated). Groups used for 2-3:
DKD, HKD, HKD_withDM."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.linear_model import RidgeCV, LogisticRegressionCV
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from kpmp import read_descriptors

STAGE = "13b_wsi_association"
EMB = ROOT / "data/processed/wsi_embeddings"
GROUPS = ["DKD", "HKD", "HKD_withDM"]
IF, TA = "Cortex: Interstitial Fibrosis: Percent", "Cortex: Tubular Atrophy: Common Type: Percent"


def slide_table() -> tuple[dict, pd.DataFrame]:
    """{(stain, tileset): DataFrame participants x emb}, features table (participant, stain, tileset)."""
    embs, feats = {}, []
    for f in sorted(EMB.glob("*.npz")):
        pid, stain = f.name.split("__")[:2]
        z = np.load(f, allow_pickle=True)
        E, F, tub = z["emb"].astype(np.float32), pd.DataFrame(z["feats"], columns=z["feat_names"]), z["tubule"]
        for ts, m in [("all", np.ones(len(E), bool)), ("tubule", tub)]:
            if m.sum() < 20:
                continue
            embs.setdefault((stain, ts), {})[pid] = E[m].mean(0)
            feats.append({"participant": pid, "stain": stain, "tileset": ts, "n_tiles": int(m.sum()),
                          **F[m].mean().to_dict()})
    return {k: pd.DataFrame(v).T for k, v in embs.items()}, pd.DataFrame(feats)


def cv_predict(X, y, cv, seed):
    pred = np.zeros(len(y))
    for tr, te in cv.split(X, y):
        m = make_pipeline(StandardScaler(), PCA(min(20, len(tr) - 2), random_state=seed),
                          RidgeCV(alphas=CFG["ridge_alphas"]))
        m.fit(X[tr], y[tr]); pred[te] = m.predict(X[te])
    return pred


def _perm_null(fn, y, seed, n_perm):
    """Label permutations drawn up front from the fixed seed, evaluated in parallel workers."""
    from joblib import Parallel, delayed
    rng = np.random.default_rng(seed)
    perms = [rng.permutation(y) for _ in range(n_perm)]
    return np.array(Parallel(n_jobs=CFG.get("n_jobs", 48))(delayed(fn)(yp) for yp in perms))


def cv_r(X, y, seed, n_perm):
    cv = KFold(min(CFG["cv_folds"], len(y)), shuffle=True, random_state=seed)
    stat = lambda yy: stats.pearsonr(cv_predict(X, yy, cv, seed), yy)[0]
    r = stat(y)
    null = _perm_null(stat, y, seed, n_perm)
    return r, float((np.sum(null >= r) + 1) / (n_perm + 1))


def cv_auc(X, y, seed, n_perm):
    def run(yy):
        cv = StratifiedKFold(min(CFG["cv_folds"], int(min(np.bincount(yy)))), shuffle=True, random_state=seed)
        p = np.zeros(len(yy))
        for tr, te in cv.split(X, yy):
            m = make_pipeline(StandardScaler(), PCA(min(20, len(tr) - 2), random_state=seed),
                              LogisticRegressionCV(Cs=10, max_iter=2000))
            m.fit(X[tr], yy[tr]); p[te] = m.predict_proba(X[te])[:, 1]
        return roc_auc_score(yy, p)
    auc = run(y)
    null = _perm_null(run, y, seed, n_perm)
    return auc, float((np.sum(null >= auc) + 1) / (n_perm + 1))


def main():
    global CFG
    cfg_all = load_config(); CFG = cfg_all["wsi"]; seed = cfg_all["seed"]; set_global_seed(seed)
    n_perm = 20 if "--quick" in sys.argv else CFG["n_perm"]   # --quick: debugging only
    out = stage_dir(STAGE)
    embs, feats = slide_table()
    feats.to_csv(out / "slide_features.tsv", sep="\t", index=False)
    desc = read_descriptors(ROOT / "data/raw/kpmp" / cfg_all["kpmp"]["descriptor_table"]).set_index("Participant ID")
    path = pd.DataFrame({k: pd.to_numeric(desc[c], errors="coerce").groupby(level=0).mean() for k, c in [("IF_pct", IF), ("TA_pct", TA)]})
    scores = {t: pd.read_csv(ROOT / f"results/07b_validate_kpmp/scores_{t}.tsv", sep="\t", index_col=0) for t in ["sn", "sc"]}
    grp = pd.concat([s.group for s in scores.values()]).groupby(level=0).first()
    rows = []
    for (stain, ts), E in sorted(embs.items()):
        # 1. sanity: pathology descriptors (all participants with a slide)
        for tgt in ["IF_pct", "TA_pct"]:
            y = path[tgt].reindex(E.index).dropna()
            if len(y) >= 15:
                r, p = cv_r(E.loc[y.index].values, y.values.astype(float), seed, n_perm)
                rows.append({"analysis": "descriptor", "stain": stain, "tileset": ts, "modality": "-", "target": tgt,
                             "n": len(y), "metric": "cv_pearson_r", "value": r, "perm_p": p})
        # 2. programme scores
        for mod, s in scores.items():
            d = s[s.group.isin(GROUPS)].join(path).join(E, how="inner")
            if len(d) < 15:
                continue
            X = d[E.columns].values
            for tgt in ["DKDiPT_up", "SharedInjury_up"]:
                r, p = cv_r(X, d[tgt].values, seed, n_perm)
                rows.append({"analysis": "programme", "stain": stain, "tileset": ts, "modality": mod, "target": tgt,
                             "n": len(d), "metric": "cv_pearson_r", "value": r, "perm_p": p,
                             "n_DKD": int((d.group == "DKD").sum()), "n_HKD": int((d.group == "HKD").sum()),
                             "n_HKD_withDM": int((d.group == "HKD_withDM").sum())})
            dd = d.dropna(subset=["IF_pct"])
            if len(dd) >= 15:
                A = np.column_stack([np.ones(len(dd)), dd.SharedInjury_up, dd.IF_pct])
                res = dd.DKDiPT_up.values - A @ np.linalg.lstsq(A, dd.DKDiPT_up.values, rcond=None)[0]
                r, p = cv_r(dd[E.columns].values, res, seed, n_perm)
                rows.append({"analysis": "programme_partial", "stain": stain, "tileset": ts, "modality": mod,
                             "target": "DKDiPT_up | SharedInjury_up, IF%", "n": len(dd), "metric": "cv_pearson_r",
                             "value": r, "perm_p": p})
        # 3. context: DKD vs HKD
        g = grp.reindex(E.index)
        m = g.isin(["DKD", "HKD"])
        if m.sum() >= 15 and (g[m] == "HKD").sum() >= 5:
            auc, p = cv_auc(E[m.values].values, (g[m] == "DKD").astype(int).values, seed, n_perm)
            rows.append({"analysis": "classification", "stain": stain, "tileset": ts, "modality": "-",
                         "target": "DKD vs HKD", "n": int(m.sum()), "metric": "cv_auc", "value": auc, "perm_p": p,
                         "n_DKD": int((g[m] == "DKD").sum()), "n_HKD": int((g[m] == "HKD").sum())})
        print(stain, ts, "done", flush=True)
    R = pd.DataFrame(rows); R.to_csv(out / "wsi_cv_results.tsv", sep="\t", index=False)
    # interpretable features: partial Spearman with DKDiPT_up adjusting SharedInjury_up and IF%
    fr = []
    fcols = [c for c in feats.columns if c not in ("participant", "stain", "tileset", "n_tiles", "tissue_frac")]
    for mod, s in scores.items():
        for (stain, ts), F in feats.groupby(["stain", "tileset"]):
            d = s[s.group.isin(GROUPS)].join(path).join(F.set_index("participant")[fcols], how="inner").dropna(subset=["IF_pct"])
            if len(d) < 15:
                continue
            rk = d[["DKDiPT_up", "SharedInjury_up", "IF_pct", *fcols]].rank()
            A = np.column_stack([np.ones(len(rk)), rk.SharedInjury_up, rk.IF_pct])
            resid = lambda v: v - A @ np.linalg.lstsq(A, v, rcond=None)[0]
            ry = resid(rk.DKDiPT_up.values)
            for c in fcols:
                rho, p = stats.pearsonr(ry, resid(rk[c].values))
                fr.append({"modality": mod, "stain": stain, "tileset": ts, "feature": c, "n": len(d), "partial_rho": rho, "p": p})
    FR = pd.DataFrame(fr)
    if len(FR):
        FR["fdr"] = FR.groupby(["modality", "stain", "tileset"]).p.transform(lambda p: multipletests(p, method="fdr_bh")[1])
    FR.to_csv(out / "feature_partial_associations.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(R.round(4).to_string()); print(FR.sort_values("p").head(15).round(4).to_string() if len(FR) else "no features")
    write_provenance(STAGE, sorted(EMB.glob("*.npz")) + [ROOT / "results/07b_validate_kpmp/scores_sn.tsv",
                     ROOT / "results/07b_validate_kpmp/scores_sc.tsv"],
                     [out / "wsi_cv_results.tsv", out / "feature_partial_associations.tsv", out / "slide_features.tsv"], seed,
                     extra={"n_slides": len(list(EMB.glob("*.npz")))})


if __name__ == "__main__":
    main()
