"""14c_published_pt_states: is the DKD-associated iPT shift the known failed-repair / adaptive PT programme?

(1) Published PT injury-state signatures (Kirita 2020, Lake 2023, Muto 2021, Abedini 2024) parsed from
    downloaded supplementary tables / verbatim full-text sentences (scripts/lib/published_signatures.py).
(2) Overlap with the 12 DKDiPT_up markers and the discovery SN-stratum iPT dE genes (hypergeometric,
    background = genes with a discovery dE test), plus the direction of each signature in the discovery dE.
(3) KPMP sn/sc iPT participant pseudobulk: published-signature scores DKD vs HKD; marker score ~ DKD +
    published score (OLS, HC3); transcriptome sign concordance after excluding published injury genes."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from annotate import zmatrix
from stats import compare
from published_signatures import all_signatures

STAGE = "14c_published_pt_states"
MARKER = "DKDiPT_up"


def hyper(sig: set, hits: set, bg: set) -> dict:
    s, h = sig & bg, hits & bg
    k, M, n, N = len(s & h), len(bg), len(h), len(s)
    exp = n * N / M if M else np.nan
    p = stats.hypergeom.sf(k - 1, M, n, N) if N and n else np.nan
    a, b, c, d = k, N - k, n - k, M - N - n + k
    orr = (a + .5) * (d + .5) / ((b + .5) * (c + .5))
    p_dep = stats.hypergeom.cdf(k, M, n, N) if N and n else np.nan
    return {"k": k, "expected": exp, "odds_ratio": orr, "p_hyper": p, "p_depletion": p_dep,
            "genes": ",".join(sorted(s & h))}


def overlap_table(sigs: dict, disc: pd.DataFrame, markers: list, fdr: float) -> pd.DataFrame:
    bg = set(disc.index[disc.dE_p.notna()])
    de = disc[disc.dE_padj < fdr]
    hits = {"markers12": set(markers), "dE_all": set(de.index),
            "dE_up": set(de.index[de.dE_lfc > 0]), "dE_down": set(de.index[de.dE_lfc < 0])}
    stat = (disc.dE_lfc / disc.dE_se).dropna()
    rows = []
    for name, genes in sigs.items():
        s = set(genes)
        r = {"signature": name, "n_genes": len(s), "n_in_background": len(s & bg), "background": len(bg)}
        for hn, h in hits.items():
            o = hyper(s, h, bg)
            r.update({f"{hn}_n": len(h & bg), **{f"{hn}_{k}": v for k, v in o.items()}})
        ins = stat[stat.index.isin(s)]
        rest = stat[~stat.index.isin(s)]
        r.update({"disc_n_scored": len(ins), "disc_frac_dE_lfc_pos": float((ins > 0).mean()) if len(ins) else np.nan,
                  "disc_median_dE_wald": float(ins.median()) if len(ins) else np.nan,
                  "disc_median_dE_wald_rest": float(rest.median()),
                  "disc_mwu_p_vs_rest": stats.mannwhitneyu(ins, rest).pvalue if len(ins) > 2 else np.nan})
        rows.append(r)
    return pd.DataFrame(rows)


def score(z: pd.DataFrame, genes, drop=()) -> tuple[pd.Series, int]:
    g = [x for x in dict.fromkeys(genes) if x in z.columns and x not in set(drop)]
    return (z[g].mean(1) if g else pd.Series(np.nan, index=z.index)), len(g)


def ols_hc3(y: pd.Series, X: pd.DataFrame) -> sm.regression.linear_model.RegressionResultsWrapper:
    return sm.OLS(y, sm.add_constant(X), missing="drop").fit(cov_type="HC3")


def std(s: pd.Series) -> pd.Series:
    return (s - s.mean()) / s.std(ddof=1)


def model_rows(tag, df, pub_cols: list[str], label: str, n_genes: dict) -> dict:
    """marker (SD units within DKD+HKD) ~ DKD + published score(s) (SD units)."""
    y = std(df["marker"])
    base = ols_hc3(y, df[["DKD"]])
    X = pd.concat([df[["DKD"]], *[std(df[c]).rename(c) for c in pub_cols]], axis=1)
    fit = ols_hc3(y, X)
    r = {"modality": tag, "model": label, "published_terms": "+".join(pub_cols),
         "n_genes_published": ";".join(str(n_genes[c]) for c in pub_cols), "n": int(fit.nobs),
         "beta_DKD_unadj": base.params["DKD"], "p_DKD_unadj": base.pvalues["DKD"],
         "beta_DKD_adj": fit.params["DKD"], "se_DKD_adj": fit.bse["DKD"], "p_DKD_adj": fit.pvalues["DKD"],
         "frac_DKD_effect_retained": fit.params["DKD"] / base.params["DKD"], "r2_adj_model": fit.rsquared}
    for c in pub_cols:
        r[f"beta[{c}]"] = fit.params[c]; r[f"p[{c}]"] = fit.pvalues[c]
    if len(pub_cols) == 1:
        c = pub_cols[0]
        r["spearman_marker_published_all"] = stats.spearmanr(df.marker, df[c]).correlation
        res = lambda v: v - v.groupby(df.DKD).transform("mean")
        r["spearman_within_group"] = stats.spearmanr(res(df.marker), res(df[c])).correlation
    return r


def kpmp_models(sigs: dict, marker_genes: list, cfg: dict, out: Path):
    tests, models = [], []
    for tag in ["sn", "sc"]:
        pb = pd.read_csv(ROOT / f"results/07c_kpmp_concordance/iPT_pseudobulk_{tag}_all_calls.tsv.gz",
                         sep="\t", index_col=0)
        pb.index = pb.index.astype(str)
        part = pd.read_csv(ROOT / f"results/07b_validate_kpmp/participants_{tag}.tsv", sep="\t", index_col=0)
        part.index = part.index.astype(str)
        grp = part.group.reindex(pb.index)
        z, _ = zmatrix(pb)
        m, nm = score(z, marker_genes)
        tests.append({"modality": tag, "signature": MARKER, "variant": "locked", "n_genes_measured": nm,
                      **compare(m, grp, "DKD", "HKD")})
        full, excl, ng_full, ng_excl = {}, {}, {}, {}
        for name, genes in sigs.items():
            full[name], ng_full[name] = score(z, genes)
            excl[name], ng_excl[name] = score(z, genes, drop=marker_genes)
            for var, s, n in [("full", full[name], ng_full[name]), ("excl_markers", excl[name], ng_excl[name])]:
                tests.append({"modality": tag, "signature": name, "variant": var, "n_genes_measured": n,
                              **compare(s, grp, "DKD", "HKD"),
                              "spearman_with_marker": stats.spearmanr(s, m, nan_policy="omit").correlation})
        keep = grp.isin(["DKD", "HKD"])
        base = pd.DataFrame({"marker": m, "DKD": (grp == "DKD").astype(float)})[keep]
        for var, S, NG in [("full", full, ng_full), ("excl_markers", excl, ng_excl)]:
            df = base.join(pd.DataFrame(S)[keep])
            for name in sigs:
                if NG[name] >= 3:
                    models.append(model_rows(tag, df, [name], f"single_{var}", NG))
            joint = [c for c in cfg["conditioning"] if NG[c] >= 3]
            models.append(model_rows(tag, df, joint, f"joint_{var}", NG))
        pd.DataFrame(full).assign(group=grp, marker=m).to_csv(out / f"scores_{tag}.tsv", sep="\t")
    return pd.DataFrame(tests), pd.DataFrame(models)


def concordance(disc: pd.DataFrame, res: pd.DataFrame, masks: dict) -> list[dict]:
    j = disc.join(res[["log2FoldChange", "padj"]].rename(columns=lambda c: f"kpmp_{c}"), how="inner")
    j = j[j.kpmp_log2FoldChange.notna()]
    rows = []
    for name, mask in masks.items():
        jj = j[mask.reindex(j.index).fillna(False).astype(bool)]
        agree = int((np.sign(jj.dE_lfc) == np.sign(jj.kpmp_log2FoldChange)).sum())
        rows.append({"gene_set": name, "n_genes": len(jj), "n_sign_agree": agree,
                     "frac_agree": agree / max(1, len(jj)),
                     "binom_p": stats.binomtest(agree, len(jj), 0.5, alternative="greater").pvalue if len(jj) else np.nan,
                     "spearman_lfc": stats.spearmanr(jj.dE_lfc, jj.kpmp_log2FoldChange).correlation if len(jj) > 2 else np.nan})
    return rows


def concordance_excluding(sigs: dict, disc: pd.DataFrame, fdr: float) -> pd.DataFrame:
    union = set().union(*sigs.values())
    de = disc.dE_padj < fdr
    tested = disc.dE_padj.notna()
    idx = disc.index.to_series()
    masks = {"disc_dE_fdr10": de,
             "disc_dE_fdr10_excl_any_published": de & ~idx.isin(union),
             "disc_dE_fdr10_in_any_published": de & idx.isin(union),
             "all_tested": tested,
             "all_tested_excl_any_published": tested & ~idx.isin(union)}
    for name, genes in sigs.items():
        masks[f"disc_dE_fdr10_excl_{name}"] = de & ~idx.isin(set(genes))
    rows = []
    for tag in ["sn", "sc"]:
        res = pd.read_csv(ROOT / f"results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_{tag}.tsv", sep="\t", index_col=0)
        rows += [{"modality": tag, "n_published_union": len(union), **r} for r in concordance(disc, res, masks)]
    return pd.DataFrame(rows)


def main():
    cfg_all = load_config(); cfg = cfg_all["published_pt_states"]; seed = cfg_all["seed"]; set_global_seed(seed)
    out = stage_dir(STAGE)
    raw = ROOT / cfg["raw_dir"]
    P = all_signatures(raw, cfg)
    P.to_csv(out / "published_signatures.tsv", sep="\t", index=False)
    sigs = {k: v.gene.tolist() for k, v in P.groupby("signature", sort=False)}
    disc = pd.read_csv(ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    gs = pd.read_csv(ROOT / "results/05_lock_programme/gene_sets.tsv", sep="\t")
    markers = gs.loc[gs.set_id == MARKER, "gene"].drop_duplicates().tolist()
    all_sigs = {**sigs, "ANY_published_union": sorted(set(P.gene))}
    O = overlap_table(all_sigs, disc, markers, cfg["dE_fdr"])
    O.to_csv(out / "overlap.tsv", sep="\t", index=False)
    T, M = kpmp_models(sigs, markers, cfg, out)
    T.to_csv(out / "published_scores_DKD_vs_HKD.tsv", sep="\t", index=False)
    M.to_csv(out / "conditioned_models.tsv", sep="\t", index=False)
    C = concordance_excluding(sigs, disc, cfg["dE_fdr"])
    C.to_csv(out / "concordance_excluding_published.tsv", sep="\t", index=False)
    man = pd.read_csv(raw / "MANIFEST.tsv", sep="\t"); man.to_csv(out / "download_manifest.tsv", sep="\t", index=False)
    inputs = [raw / f for f in man.file] + [
        ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv", ROOT / "results/05_lock_programme/gene_sets.tsv",
        *[ROOT / f"results/07c_kpmp_concordance/iPT_pseudobulk_{t}_all_calls.tsv.gz" for t in ["sn", "sc"]],
        *[ROOT / f"results/07b_validate_kpmp/participants_{t}.tsv" for t in ["sn", "sc"]],
        *[ROOT / f"results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_{t}.tsv" for t in ["sn", "sc"]]]
    outputs = [out / f for f in ["published_signatures.tsv", "overlap.tsv", "published_scores_DKD_vs_HKD.tsv",
                                 "conditioned_models.tsv", "concordance_excluding_published.tsv",
                                 "scores_sn.tsv", "scores_sc.tsv", "download_manifest.tsv"]]
    write_provenance(STAGE, inputs, outputs, seed, extra={
        "n_signatures": len(sigs), "not_obtained": "Gerhardt et al. 2023 JASN doi:10.1681/ASN.0000000000000057 "
        "(not in PMC open-access subset; publisher/PMC pages refuse scripted access)"})
    print(O[["signature", "n_in_background", "markers12_k", "markers12_genes", "dE_all_k", "dE_all_expected",
             "dE_all_p_hyper", "disc_frac_dE_lfc_pos"]].to_string())
    print(T.to_string()); print(M.iloc[:, :12].to_string()); print(C.to_string())


if __name__ == "__main__":
    main()
