"""12c_kpmp_clinical: clinical and pathology correlates of the iPT programme scores in KPMP.

Participants: DKD, HKD, HKD_withDM (scores from 07b, per modality). Public bands are used as
ordinal ranks (never midpoints); TIV descriptor grades 999 = not assessable -> missing.
(a) Spearman score vs variable (pooled); (b) score ~ group + rank(variable), HC3 SE.
Glycaemic variables (A1c band, diabetes duration) are analysed in diabetic participants only
(DKD + HKD_withDM). BH-FDR within each modality x analysis family."""
from __future__ import annotations
import importlib, sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from kpmp import read_descriptors

STAGE = "12c_kpmp_clinical"
RAW = ROOT / "data/raw/kpmp"
BANDS = {
    "proteinuria": ("Proteinuria (mg) (Binned)", ["<150 mg/g cr", "150 to <500 mg/g cr", "500 to <1000 mg/g cr", ">=1000 mg/g cr"]),
    "albuminuria": ("Albuminuria (mg) (Binned)", ["<30 mg/g cr", "30 to <300 mg/g cr", "300 to <500 mg/g cr",
                                                  "500 to <1000 mg/g cr", ">=1000 mg/g cr"]),
    "a1c": ("A1c (%) (Binned)", ["<6.5%", "6.5 to <7.5%", "7.5 to <8.5%", ">=8.5%"]),
}
PATH = {"if_pct": "Cortex: Interstitial Fibrosis: Percent", "ta_pct": "Cortex: Tubular Atrophy: Common Type: Percent",
        "tubular_injury_pct": "Cortex: Tubular Injury Other: Percent",
        "interstitial_wbc_pct": "Cortex: Interstitial Monuclear WBC: Percent",
        "arteriosclerosis_grade": "Cortex and Medulla: Arteriosclerosis: Arteriosclerosis: Presence",
        "arteriolar_hyalinosis_grade": "Cortex and Medulla: Arteriolar Hyalinosis: Arteriolar Hyalinosis: Presence"}
GLYC = {"a1c", "diabetes_duration"}
SCORES = ["DKDiPT_up", "SharedInjury_up"]
GROUPS = ["DKD", "HKD", "HKD_withDM"]


def lower_bound(x) -> float:
    try:
        return float(str(x).split("-")[0].split()[0])
    except ValueError:
        return np.nan


def clinical_table(cfg) -> pd.DataFrame:
    c = pd.read_csv(RAW / cfg["kpmp"]["clinical_table"]).set_index("Participant ID")
    t = pd.DataFrame(index=c.index)
    t["egfr_band"] = c["Baseline eGFR (ml/min/1.73m2) (Binned)"].map(lower_bound)
    for k, (col, order) in BANDS.items():
        t[k] = c[col].map({v: i for i, v in enumerate(order)})
    t["diabetes_duration"] = c["Diabetes Duration (Years)"].map(lower_bound)
    t["raas"] = c["On RAAS Blockade"].map({"Yes": 1, "No": 0})
    d = read_descriptors(RAW / cfg["kpmp"]["descriptor_table"]).set_index("Participant ID")
    for k, col in PATH.items():
        v = pd.to_numeric(d[col], errors="coerce")
        v = v.where(v != 999)
        t[k] = v.groupby(level=0).mean().reindex(t.index)
    return t


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); R = cfg["robustness"]
    out = stage_dir(STAGE)
    clin = clinical_table(cfg)
    rows = []
    for tag in ["sn", "sc"]:
        s = pd.read_csv(ROOT / f"results/07b_validate_kpmp/scores_{tag}.tsv", sep="\t", index_col=0)
        d = s[s.group.isin(GROUPS)][SCORES + ["group"]].join(clin)
        d.to_csv(out / f"analysis_table_{tag}.tsv", sep="\t")
        for var in clin.columns:
            dd = d[d.group.isin(["DKD", "HKD_withDM"])] if var in GLYC else d
            for sc_ in SCORES:
                x = dd[[sc_, var, "group"]].dropna()
                if len(x) < R["clinical_min_n"] or x[var].nunique() < 2:
                    continue
                rr = stats.spearmanr(x[sc_], x[var])
                rows.append({"modality": tag, "analysis": "spearman", "score": sc_, "variable": var, "n": len(x),
                             "estimate": rr.correlation, "p": rr.pvalue})
                x = x.assign(r=x[var].rank(pct=True), y=x[sc_])
                form = "y ~ C(group) + r" if x.group.nunique() > 1 else "y ~ r"
                m = smf.ols(form, x).fit(cov_type="HC3")
                rows.append({"modality": tag, "analysis": "group_adjusted_rank_ols", "score": sc_, "variable": var,
                             "n": len(x), "estimate": m.params["r"], "p": m.pvalues["r"]})
    T = pd.DataFrame(rows)
    T["fdr_bh"] = np.nan
    for _, idx in T.groupby(["modality", "analysis"]).groups.items():
        T.loc[idx, "fdr_bh"] = multipletests(T.loc[idx, "p"], method="fdr_bh")[1]
    T.to_csv(out / "clinical_associations.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_rows", 200)
    print(T.round(4).to_string())
    write_provenance(STAGE, [RAW / cfg["kpmp"]["clinical_table"], RAW / cfg["kpmp"]["descriptor_table"],
                             ROOT / "results/07b_validate_kpmp/scores_sn.tsv", ROOT / "results/07b_validate_kpmp/scores_sc.tsv"],
                     sorted(out.glob("*.tsv")), seed)


if __name__ == "__main__":
    main()
