"""14h_injury_intensity_adjustment: are the albuminuria association and the HKD -> HKD with diabetes
-> DKD trend of the marker-gene score explained by injury intensity within iPT cells?

Injury intensity = score of the discovery-atlas iPT signature (Abedini 2024, Supplementary Table 8;
results/14c scores). KPMP single-nucleus participants with DKD, HKD or HKD with diabetes.
(1) z(albuminuria rank) ~ z(marker) [+ z(iPT signature)] + diagnosis, OLS HC3; also within DKD.
(2) ordered trend: marker ~ ordinal group code + z(iPT signature) (OLS HC3; one-sided P for the
    group coefficient), and Kendall tau of the marker score residualised on the iPT signature,
    with label permutation."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "14h_injury_intensity_adjustment"
ORDER = {"HKD": 0, "HKD_withDM": 1, "DKD": 2}
z = lambda v: (v - v.mean()) / v.std(ddof=1)


def fit(form, d, term):
    m = smf.ols(form, d).fit(cov_type="HC3")
    ci = m.conf_int().loc[term]
    return {"n": int(m.nobs), "beta": m.params[term], "ci_lo": ci[0], "ci_hi": ci[1], "p": m.pvalues[term]}


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); c = cfg["injury_intensity"]
    f_sig = ROOT / "results/14c_published_pt_states/scores_sn.tsv"
    f_cl = ROOT / "results/12c_kpmp_clinical/analysis_table_sn.tsv"
    sig = pd.read_csv(f_sig, sep="\t", index_col=0)[[c["signature"]]].rename(columns={c["signature"]: "ipt_sig"})
    cl = pd.read_csv(f_cl, sep="\t", index_col=0)
    d = cl.join(sig, how="inner")
    d = d[d.group.isin(ORDER)].copy()
    rows = []
    # (1) albuminuria
    a = d.dropna(subset=["albuminuria", "DKDiPT_up", "ipt_sig"]).copy()
    a["z_out"] = z(a.albuminuria.rank(pct=True)); a["z_marker"] = z(a.DKDiPT_up); a["z_sig"] = z(a.ipt_sig)
    for scope, dd in [("all", a), ("DKD_only", a[a.group == "DKD"])]:
        g = " + C(group)" if scope == "all" else ""
        rows.append({"analysis": "albuminuria", "scope": scope, "model": "marker", **fit(f"z_out ~ z_marker{g}", dd, "z_marker")})
        rows.append({"analysis": "albuminuria", "scope": scope, "model": "marker + iPT signature",
                     **fit(f"z_out ~ z_marker + z_sig{g}", dd, "z_marker")})
        rows.append({"analysis": "albuminuria", "scope": scope, "model": "iPT signature alone", **fit(f"z_out ~ z_sig{g}", dd, "z_sig")})
    # (2) ordered trend
    t = d.dropna(subset=["DKDiPT_up", "ipt_sig"]).copy()
    t["code"] = t.group.map(ORDER); t["z_marker"] = z(t.DKDiPT_up); t["z_sig"] = z(t.ipt_sig)
    for model, form in [("marker", "z_marker ~ code"), ("marker + iPT signature", "z_marker ~ code + z_sig")]:
        r = fit(form, t, "code"); r["p_one_sided"] = r["p"] / 2 if r["beta"] > 0 else 1 - r["p"] / 2
        rows.append({"analysis": "exposure_trend", "scope": "all", "model": model, **r})
    resid = smf.ols("z_marker ~ z_sig", t).fit().resid.values
    tau = stats.kendalltau(t.code, resid).statistic
    rng = np.random.default_rng(seed)
    null = np.array([stats.kendalltau(rng.permutation(t.code.values), resid).statistic for _ in range(c["n_perm"])])
    rows.append({"analysis": "exposure_trend", "scope": "all", "model": "Kendall tau, marker residualised on iPT signature",
                 "n": len(t), "beta": tau, "p_one_sided": float((np.sum(null >= tau) + 1) / (len(null) + 1))})
    out = stage_dir(STAGE); T = pd.DataFrame(rows); T.to_csv(out / "injury_intensity_adjustment.tsv", sep="\t", index=False)
    print(T.round(4).to_string())
    write_provenance(STAGE, [f_sig, f_cl], [out / "injury_intensity_adjustment.tsv"], seed)


if __name__ == "__main__":
    main()
