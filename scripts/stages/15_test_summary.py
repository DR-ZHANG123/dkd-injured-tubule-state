"""15_test_summary: one table of the pre-specified tests and the main exploratory tests, with
effect size, 95% CI (primary KPMP arms: participant bootstrap; other Hedges g: large-sample SE), P and source file. No new statistics are
computed beyond the CI of already reported effect sizes."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "15_test_summary"
R = ROOT / "results"


def g_ci(g, na, nb):
    se = np.sqrt((na + nb) / (na * nb) + g ** 2 / (2 * (na + nb)))
    return g - 1.96 * se, g + 1.96 * se


def boot_g_ci(scores: pd.Series, groups: pd.Series, a: str, b: str, n_boot: int, seed: int):
    """Participant bootstrap (resampling within group) percentile CI for Hedges g."""
    from stats import hedges_g
    rng = np.random.default_rng(seed)
    xa, xb = scores[groups == a].dropna().values, scores[groups == b].dropna().values
    bs = [hedges_g(rng.choice(xa, len(xa)), rng.choice(xb, len(xb))) for _ in range(n_boot)]
    return tuple(np.nanpercentile(bs, [2.5, 97.5]))


def pick(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k].astype(str) == str(v)).values
    assert m.sum() == 1, (kw, m.sum())
    return df[m].iloc[0]


def main():
    cfg = load_config(); set_global_seed(cfg["seed"])
    rows, srcs = [], []

    def add(analysis, status, cohort, comparison, n, effect_name, eff, lo, hi, p, src):
        rows.append({"analysis": analysis, "status": status, "cohort": cohort, "comparison": comparison, "n": n,
                     "effect": effect_name, "estimate": eff, "ci_low": lo, "ci_high": hi, "p": p, "source": src})
        srcs.append(R / src)

    f = "07b_validate_kpmp/programme_tests.tsv"; t = pd.read_csv(R / f, sep="\t")
    for mod, lab in [("sn", "KPMP single-nucleus"), ("sc", "KPMP single-cell")]:
        r = pick(t, modality=mod, set="DKDiPT_up", contrast="DKD-HKD")
        sc_ = pd.read_csv(R / f"07b_validate_kpmp/scores_{mod}.tsv", sep="\t", index_col=0)
        srcs.append(R / f"07b_validate_kpmp/scores_{mod}.tsv")
        lo, hi = boot_g_ci(sc_.DKDiPT_up, sc_.group, "DKD", "HKD", cfg["test_summary"]["n_boot"], cfg["seed"])
        add("Marker-gene score, DKD vs HKD", "pre-specified primary", lab, "DKD vs HKD", f"{int(r.n_a)} vs {int(r.n_b)}",
            "Hedges g", r.hedges_g, lo, hi, r.welch_p, f)
    f = "07b_validate_kpmp/pooled_effects.tsv"; r = pick(pd.read_csv(R / f, sep="\t"), set="DKDiPT_up", contrast="DKD-HKD")
    add("Marker-gene score, DKD vs HKD", "pre-specified primary (pooled)", "KPMP both modalities", "DKD vs HKD", "random effects",
        "Hedges g", r.g_pooled, r.ci_lo, r.ci_hi, r.p, f)
    f = "06_validate_gse195460/programme_tests.tsv"; t = pd.read_csv(R / f, sep="\t")
    r = pick(t, compartment="iPT", set="DKDiPT_up"); lo, hi = g_ci(r.hedges_g, r.n_a, r.n_b)
    add("Marker-gene score", "pre-specified", "GSE195460 iPT", "DKD vs control", f"{int(r.n_a)} vs {int(r.n_b)}",
        "Hedges g", r.hedges_g, lo, hi, r.welch_p, f)
    f = "08_validate_ercb_arrays/programme_tests.tsv"; t = pd.read_csv(R / f, sep="\t")
    r = pick(t, platform="GPL24120", set="DKDiPT_up", contrast="DKD-HTN"); lo, hi = g_ci(r.hedges_g, r.n_a, r.n_b)
    add("Marker-gene score (5 of 12 genes)", "pre-specified", "ERCB tubulointerstitium", "DKD vs hypertensive nephropathy",
        f"{int(r.n_a)} vs {int(r.n_b)}", "Hedges g", r.hedges_g, lo, hi, r.welch_p, f)
    f = "10b_visium_kpmp/participant_tests.tsv"; t = pd.read_csv(R / f, sep="\t")
    r = pick(t, metric="DKDiPT_up", contrast="DKD-HKD"); lo, hi = g_ci(r.hedges_g, r.n_a, r.n_b)
    add("Marker-gene score in PT-rich spots", "pre-specified", "KPMP Visium", "DKD vs HKD", f"{int(r.n_a)} vs {int(r.n_b)}",
        "Hedges g", r.hedges_g, lo, hi, r.welch_p, f)
    f = "13b_wsi_association/wsi_cv_results.tsv"; t = pd.read_csv(R / f, sep="\t")
    r = pick(t, analysis="programme_partial", stain="HE", tileset="all", modality="sn")
    add("Whole-slide H&E prediction of marker score (adjusted)", "pre-specified primary histology", "KPMP whole slides",
        "cross-validated", int(r.n), "Pearson r", r.value, np.nan, np.nan, r.perm_p, f)
    f = "14a_concordance_permutation/permutation_concordance.tsv"; t = pd.read_csv(R / f, sep="\t")
    for mod, con, lab in [("sn", "DKD_vs_HKD", "KPMP single-nucleus"), ("sn", "DKD_vs_HKD_adj_egfr", "KPMP single-nucleus, eGFR-adjusted"),
                          ("sc", "DKD_vs_HKD", "KPMP single-cell")]:
        sub = t[(t.modality == mod) & (t.contrast.str.startswith(con)) & (t.gene_set == "disc_dE_fdr10") & (t.primary.astype(str) == "True")]
        if con == "DKD_vs_HKD":
            sub = sub[sub.contrast == "DKD_vs_HKD"]
        r = sub.iloc[0]
        add("Transcriptome-wide same-direction fraction", "exploratory", lab, "DKD vs HKD", int(r.n_genes),
            "fraction (permutation null mean)", r.observed_frac, r.perm_null_mean, r.perm_null_q95, r.perm_p, f)
    f = "14b_albuminuria_models/albuminuria_models.tsv"
    if (R / f).exists():
        srcs.append(R / f)
    out = stage_dir(STAGE)
    T = pd.DataFrame(rows)
    T.to_csv(out / "test_summary.tsv", sep="\t", index=False)
    print(T.round(3).to_string())
    write_provenance(STAGE, sorted(set(srcs)), [out / "test_summary.tsv"], cfg["seed"])


if __name__ == "__main__":
    main()
