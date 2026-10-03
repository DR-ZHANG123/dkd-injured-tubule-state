"""14b_albuminuria_models: covariate-adjusted albuminuria/proteinuria models for the iPT scores (KPMP).

Inputs are the participant tables written by 12c (DKDiPT_up = DKD-biased marker score, SharedInjury_up =
common injury score). Bands are ordinal: entered as percentile ranks (never midpoints) within each model's
complete-case sample, then z-scored together with the scores, so coefficients are standardized and
comparable between the two scores (raw score units differ in spread).
(1) z(rank outcome) ~ z(score) + C(group) [+ eGFR rank, IF% rank, RAAS], HC3 SE; reverse orientation
    z(score) ~ z(rank outcome) + ...; the 12c raw-unit estimate (score ~ C(group) + pct rank) is kept alongside.
(2) marker vs injury: joint model z(outcome) ~ z(marker) + z(injury) + C(group) with an HC3 Wald test of
    b_marker - b_injury, plus a group-stratified participant bootstrap of the difference in standardized
    coefficients from the two separate models.
Scopes: all diagnoses, DKD only, excluding HKD_withDM. sn is primary, sc secondary."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "14b_albuminuria_models"
IN = ROOT / "results/12c_kpmp_clinical"
MARK, INJ = "DKDiPT_up", "SharedInjury_up"
SCORES = [MARK, INJ]


def zs(v) -> np.ndarray:
    v = np.asarray(v, float)
    sd = v.std(ddof=1)
    return (v - v.mean()) / sd if sd > 0 else v * np.nan


def prep(d: pd.DataFrame, outcome: str, covs: dict) -> tuple[pd.DataFrame, list[str]]:
    """Complete-case frame: z_out, z_<score>, raw scores, pct-rank covariates, group. Drops constant covariates."""
    need = [outcome, *SCORES, "group", *covs.values()]
    x = d[need].dropna().copy()
    x["r_out"] = x[outcome].rank(pct=True)
    x["z_out"] = zs(x["r_out"])
    for s in SCORES:
        x[f"z_{s}"] = zs(x[s])
    terms = []
    for k, col in covs.items():
        if x[col].nunique() < 2:
            continue
        x[k] = x[col] if col == "raas" else x[col].rank(pct=True)
        terms.append(k)
    return x, terms


def rhs(x: pd.DataFrame, terms: list[str]) -> str:
    g = ["C(group)"] if x.group.nunique() > 1 else []
    return " + ".join(g + terms)


def ols(form: str, x: pd.DataFrame):
    return smf.ols(form, x).fit(cov_type="HC3")


def row(m, term: str) -> dict:
    ci = m.conf_int().loc[term]
    return {"estimate": m.params[term], "se_hc3": m.bse[term], "ci_lo": ci[0], "ci_hi": ci[1], "p": m.pvalues[term]}


def cov_sets(covs: dict) -> dict:
    s = {"diagnosis": {}}
    for k, col in covs.items():
        s[f"diagnosis+{k}"] = {k: col}
    s["diagnosis+all"] = dict(covs)
    return s


def model_rows(d, tag, outcome, scope, covs, min_n) -> list[dict]:
    rows = []
    for cname, cset in cov_sets(covs).items():
        x, terms = prep(d, outcome, cset)
        if len(x) < min_n or x[outcome].nunique() < 2:
            continue
        r = rhs(x, terms)
        base = {"modality": tag, "outcome": outcome, "scope": scope, "covariates": cname,
                "covariates_used": ",".join(terms) or "-", "n": len(x),
                "n_by_group": ";".join(f"{k}={v}" for k, v in x.group.value_counts().sort_index().items())}
        for s in SCORES:
            f1 = f"z_out ~ z_{s}" + (f" + {r}" if r else "")
            rows.append({**base, "score": s, "orientation": "outcome_on_score", **row(ols(f1, x), f"z_{s}")})
            f2 = f"z_{s} ~ z_out" + (f" + {r}" if r else "")
            rr = {**base, "score": s, "orientation": "score_on_outcome", **row(ols(f2, x), "z_out")}
            m3 = ols(f"{s} ~ r_out" + (f" + {r}" if r else ""), x)
            rr.update(estimate_raw=m3.params["r_out"], se_raw_hc3=m3.bse["r_out"], p_raw=m3.pvalues["r_out"])
            rows.append(rr)
    return rows


def design(x: pd.DataFrame, pred: list[str], groups: list[str]) -> np.ndarray:
    cols = [np.ones(len(x))] + [x[p].to_numpy() for p in pred]
    cols += [(x.group == g).to_numpy(float) for g in groups[1:]]
    return np.column_stack(cols)


def boot_once(x: pd.DataFrame, groups: list[str]) -> dict | None:
    """Standardized coefficients recomputed from scratch (ranks + z) in one resample."""
    x = x.copy()
    x["z_out"] = zs(x["r_raw"].rank(pct=True))
    for s in SCORES:
        x[f"z_{s}"] = zs(x[s])
    if x[["z_out", *[f"z_{s}" for s in SCORES]]].isna().any().any():
        return None
    out = {}
    for s in SCORES:
        X = design(x, [f"z_{s}"], groups)
        out[f"os_{s}"] = np.linalg.lstsq(X, x["z_out"].to_numpy(), rcond=None)[0][1]
        X = design(x, ["z_out"], groups)
        out[f"so_{s}"] = np.linalg.lstsq(X, x[f"z_{s}"].to_numpy(), rcond=None)[0][1]
    X = design(x, [f"z_{MARK}", f"z_{INJ}"], groups)
    b = np.linalg.lstsq(X, x["z_out"].to_numpy(), rcond=None)[0]
    out["joint_m"], out["joint_i"] = b[1], b[2]
    return out


def bootstrap(x: pd.DataFrame, n_boot: int, rng) -> pd.DataFrame:
    groups = sorted(x.group.unique())
    idx_by_g = [np.flatnonzero(x.group.to_numpy() == g) for g in groups]
    res = []
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(ix, len(ix), replace=True) for ix in idx_by_g])
        o = boot_once(x.iloc[idx], groups)
        if o is not None:
            res.append(o)
    return pd.DataFrame(res)


def boot_summary(v: np.ndarray) -> dict:
    v = v[np.isfinite(v)]
    p = 2 * min((v <= 0).mean(), (v >= 0).mean())
    return {"ci_lo": np.quantile(v, 0.025), "ci_hi": np.quantile(v, 0.975), "p_boot": min(1.0, p),
            "frac_positive": (v > 0).mean()}


def slope_rows(d, tag, outcome, scope, n_boot, min_n, seed) -> list[dict]:
    x, _ = prep(d, outcome, {})
    if len(x) < min_n or x[outcome].nunique() < 2:
        return []
    x["r_raw"] = x[outcome]
    r = rhs(x, [])
    base = {"modality": tag, "outcome": outcome, "scope": scope, "n": len(x)}
    rng = np.random.default_rng(seed)
    B = bootstrap(x, n_boot, rng)
    rows = []
    # (a) difference of standardized coefficients from separate models, both orientations
    for key, lab, lhs_fmt in [("os", "outcome_on_score", "z_out ~ z_{s}"), ("so", "score_on_outcome", "z_{s} ~ z_out")]:
        est = {}
        for s in SCORES:
            f = lhs_fmt.format(s=s) + (f" + {r}" if r else "")
            est[s] = ols(f, x).params["z_out" if key == "so" else f"z_{s}"]
        diff = B[f"{key}_{MARK}"] - B[f"{key}_{INJ}"]
        rows.append({**base, "test": f"separate_models_{lab}", "beta_marker": est[MARK], "beta_injury": est[INJ],
                     "difference": est[MARK] - est[INJ], "n_boot_valid": len(B), **boot_summary(diff.to_numpy())})
    # (b) joint model: both scores together, HC3 Wald test of b_marker - b_injury, plus bootstrap CI
    f = f"z_out ~ z_{MARK} + z_{INJ}" + (f" + {r}" if r else "")
    m = ols(f, x)
    wt = m.t_test(f"z_{MARK} - z_{INJ} = 0")
    jm, ji = row(m, f"z_{MARK}"), row(m, f"z_{INJ}")
    rows.append({**base, "test": "joint_model_outcome_on_both", "beta_marker": jm["estimate"], "beta_injury": ji["estimate"],
                 "difference": float(np.squeeze(wt.effect)), "p_marker": jm["p"], "p_injury": ji["p"],
                 "p_wald_diff_hc3": float(np.squeeze(wt.pvalue)), "n_boot_valid": len(B),
                 **boot_summary((B["joint_m"] - B["joint_i"]).to_numpy())})
    return rows


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); A = cfg["albuminuria_models"]
    out = stage_dir(STAGE)
    mrows, srows = [], []
    for j, tag in enumerate(["sn", "sc"]):
        d0 = pd.read_csv(IN / f"analysis_table_{tag}.tsv", sep="\t", index_col=0)
        for scope, groups in A["scopes"].items():
            d = d0[d0.group.isin(groups)]
            for k, outcome in enumerate(A["outcomes"]):
                mrows += model_rows(d, tag, outcome, scope, A["covariates"], A["min_n"])
                srows += slope_rows(d, tag, outcome, scope, A["n_boot"], A["min_n"], seed + 100 * j + 10 * k)
    M = pd.DataFrame(mrows)
    M["fdr_bh"] = np.nan  # BH within modality x outcome x orientation (all scopes, covariate sets, scores)
    for _, idx in M.groupby(["modality", "outcome", "orientation"]).groups.items():
        M.loc[idx, "fdr_bh"] = multipletests(M.loc[idx, "p"], method="fdr_bh")[1]
    S = pd.DataFrame(srows)
    M.to_csv(out / "albuminuria_models.tsv", sep="\t", index=False)
    S.to_csv(out / "slope_difference.tsv", sep="\t", index=False)
    (out / "README.md").write_text(
        "# 14b_albuminuria_models\n"
        "albuminuria_models.tsv: rank-OLS (HC3) of albuminuria/proteinuria band rank on each iPT score + diagnosis, "
        "then + eGFR band rank, + cortical IF% rank, + RAAS, + all; both orientations; standardized (z) coefficients, "
        "n per complete-case sample; estimate_raw = 12c raw-unit score ~ diagnosis + pct rank.\n"
        "slope_difference.tsv: marker (DKDiPT_up) minus common injury (SharedInjury_up) standardized coefficient; "
        f"separate models with {A['n_boot']} diagnosis-stratified participant bootstrap resamples (seeded), and joint model with HC3 Wald test.\n"
        "Scopes: all = DKD+HKD+HKD_withDM; DKD_only (no diagnosis term); no_HKD_withDM = DKD+HKD. sn primary, sc secondary.\n"
        "Bands are ordinal percentile ranks (no midpoints); p values are exploratory.\n")
    pd.set_option("display.width", 250, "display.max_rows", 500, "display.max_columns", 30)
    print(M.drop(columns=["n_by_group"]).round(4).to_string())
    print(S.round(4).to_string())
    write_provenance(STAGE, [IN / "analysis_table_sn.tsv", IN / "analysis_table_sc.tsv"],
                     [out / "albuminuria_models.tsv", out / "slope_difference.tsv", out / "README.md"], seed,
                     extra={"n_boot": A["n_boot"], "bootstrap": "stratified by diagnosis within scope"})


if __name__ == "__main__":
    main()
