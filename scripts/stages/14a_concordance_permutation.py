"""14a_concordance_permutation: calibrated nulls for the KPMP iPT transcriptome concordance (07c).

The 07c binomial test compares the same-sign fraction of discovery dE genes against 50%, but
genes are correlated and ~59% of all tested genes already agree. Two nulls are computed here:
(1) participant-label permutation: DKD/HKD labels permuted among KPMP participants (within
    modality, group sizes kept), contrast refitted with a fast statistic each time;
(2) background-matched gene sets: random gene sets of the same size drawn from all tested genes,
    matched on discovery baseMean decile, scored with the observed KPMP DESeq2 log2FC.
Fast statistic (primary): the DESeq2 NB GLM (~ sex [+ eGFR_lower z] + group, offset log size
factor) refitted by vectorised IRLS with size factors and dispersions fixed at the DESeq2 estimates
on the observed labels (scripts/lib/perm_stats.py). The simpler log2 ratio of normalised group means
(+ pseudocount) is kept as a sensitivity row; it fails the >95% sign-agreement criterion in sn.
eGFR-adjusted contrast: labels permuted within eGFR_lower bands (config egfr_strata_edges).
The first n_deseq_check permutations are refitted with full DESeq2 (07c deseq_contrast)"""
from __future__ import annotations
import importlib.util, sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from perm_stats import deseq_fit, emp_p, mean_ratio_lfc, nb_glm_lfc, permute_labels, row_spearman

STAGE = "14a_concordance_permutation"
C07 = ROOT / "results/07c_kpmp_concordance"


def load_07c():
    spec = importlib.util.spec_from_file_location("s07c", ROOT / "scripts/stages/07c_kpmp_concordance.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def participant_meta(cfg, tag: str) -> pd.DataFrame:
    q = ad.read_h5ad(ROOT / cfg["paths"]["processed"] / f"KPMP_{tag}_qc_counts.h5ad", backed="r")
    meta = q.obs.groupby("participant", observed=True).agg(group=("group", "first"), sex=("sex", "first"))
    q.file.close()
    meta.index = meta.index.astype(str)
    meta["group"] = meta.group.astype(str); meta["sex"] = meta.sex.astype(str)
    desc = pd.read_csv(ROOT / f"results/07b_validate_kpmp/participants_{tag}.tsv", sep="\t", index_col=0)
    desc.index = desc.index.astype(str)
    return meta.join(desc[["egfr_lower"]])


def filtered_counts(pb: pd.DataFrame, md: pd.DataFrame, cfg) -> pd.DataFrame:
    """Same participant set and gene filter as 07c deseq_contrast."""
    p = cfg["pseudobulk"]
    c = pb.loc[md.index]
    return c.loc[:, (c >= p["min_counts_gene"]).sum(0) >= p["min_samples_gene"]].round().astype(int)


def background_null(ref: pd.Series, lfc: pd.Series, basemean: pd.Series, set_genes: list[str],
                    n_draw: int, n_bins: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """Random gene sets of len(set_genes) from all tested genes, matched on baseMean decile."""
    pool = ref.index.intersection(lfc.dropna().index).intersection(basemean.dropna().index)
    bins = pd.qcut(basemean[pool].rank(method="first"), n_bins, labels=False)
    need = bins[set_genes].value_counts()
    agree = (np.sign(ref[pool]) == np.sign(lfc[pool])).values
    pos = {k: np.where(bins.values == k)[0] for k in need.index}
    fr, rho = np.empty(n_draw), np.empty(n_draw)
    rv, lv = ref[pool].values, lfc[pool].values
    for i in range(n_draw):
        idx = np.concatenate([rng.choice(pos[k], int(m), replace=False) for k, m in need.items()])
        fr[i] = agree[idx].mean()
        rho[i] = stats.spearmanr(rv[idx], lv[idx]).correlation
    return fr, rho


def gene_sets(disc, lin, sets):
    shared = [g for g in sets["SharedInjury_up"] + sets["SharedInjury_down"] if g in lin.index]
    lin_ref = ((lin.dD_lfc + lin.dN_lfc) / 2).dropna()
    return {"disc_dE_fdr10": (disc.dE_lfc[disc.dE_padj.notna()], disc.baseMean,
                              disc.index[disc.dE_padj < 0.10].tolist()),
            "shared_injury_genes": (lin_ref, lin.baseMean, shared)}


def set_rows(tag, cname, stat_name, primary, S, genes, res, lab, obs_fast, null_fast, bnull):
    rows = []
    for sname, (ref, _, _) in S.items():
        sg = bnull[sname]["genes"]
        r = ref[sg].values; col = genes.get_indexer(sg); dlfc = res.log2FoldChange[sg].values
        f_obs = obs_fast[col]; nf = null_fast[:, col]
        of_frac = float((np.sign(r) == np.sign(f_obs)).mean())
        of_rho = float(stats.spearmanr(r, f_obs).correlation)
        pfrac = (np.sign(nf) == np.sign(r)[None]).mean(1); prho = row_spearman(r, nf)
        b = bnull[sname]
        rows.append({"modality": tag, "contrast": cname, "gene_set": sname, "n_genes": len(sg),
                     "n_dkd": int(lab.sum()), "n_hkd": int((~lab).sum()), "fast_stat": stat_name,
                     "primary": primary, "observed_frac": b["obs_frac"], "observed_frac_fast": of_frac,
                     "perm_null_mean": pfrac.mean(), "perm_null_q95": np.quantile(pfrac, 0.95),
                     "perm_p": emp_p(pfrac, of_frac),
                     "background_null_mean": b["frac"].mean(), "background_null_q95": np.quantile(b["frac"], 0.95),
                     "background_p": emp_p(b["frac"], b["obs_frac"]),
                     "observed_rho": b["obs_rho"], "observed_rho_fast": of_rho,
                     "rho_perm_null_mean": prho.mean(), "rho_perm_null_q95": np.quantile(prho, 0.95),
                     "rho_perm_p": emp_p(prho, of_rho),
                     "rho_background_null_mean": np.nanmean(b["rho"]), "rho_background_p": emp_p(b["rho"], b["obs_rho"]),
                     "sign_agree_fast_vs_deseq2": float((np.sign(f_obs) == np.sign(dlfc)).mean()),
                     "n_perm": len(pfrac), "n_background": len(b["frac"])})
    return rows


def run_contrast(tag, cname, pb, meta, cfg, S, rng, s07c, checks):
    pc_cfg = cfg["permutation_concordance"]
    adj = cname.endswith("egfr_lower"); covars = ["egfr_lower"] if adj else []
    md = meta.loc[meta.index.intersection(pb.index)]
    md = md[md.group.isin(["DKD", "HKD"])].dropna(subset=covars)
    counts = filtered_counts(pb, md, cfg)
    res = pd.read_csv(C07 / f"kpmp_iPT_{cname}_{tag}.tsv", sep="\t", index_col=0)
    sf, disp, lfc_refit = deseq_fit(counts, md, covars, "DKD", "HKD")
    refit_r = float(np.corrcoef(lfc_refit.reindex(res.index), res.log2FoldChange)[0, 1])
    lab = (md.group == "DKD").values
    X0 = np.column_stack([np.ones(len(md)), (md.sex == "Male").astype(float)])
    if adj:
        z = (md.egfr_lower - md.egfr_lower.mean()) / md.egfr_lower.std(ddof=1)
        X0 = np.column_stack([X0, z.values])
        strata = np.digitize(md.egfr_lower.values, pc_cfg["egfr_strata_edges"])
    else:
        strata = np.zeros(len(md), int)
    perms = permute_labels(lab, strata, pc_cfg["n_perm"], rng)
    bnull = {}
    for sname, (ref, bm, sgenes) in S.items():
        sg = [g for g in sgenes if g in ref.index and g in counts.columns and pd.notna(res.log2FoldChange.get(g))]
        r, d = ref[sg].values, res.log2FoldChange[sg].values
        fr, rho = background_null(ref, res.log2FoldChange, bm, sg, pc_cfg["n_background"], pc_cfg["basemean_bins"], rng)
        bnull[sname] = {"genes": sg, "frac": fr, "rho": rho, "obs_frac": float((np.sign(r) == np.sign(d)).mean()),
                        "obs_rho": float(stats.spearmanr(r, d).correlation)}
    genes = pd.Index(sorted(set().union(*[b["genes"] for b in bnull.values()])))
    y = counts[genes].values.astype(float); dsp = disp[genes].values
    nb = lambda L: nb_glm_lfc(y, sf, dsp, X0, L)
    obs_nb, null_nb = nb(lab[None])[0], nb(perms)
    rows = set_rows(tag, cname, "nb_glm_fixed_dispersion", True, S, genes, res, lab, obs_nb, null_nb, bnull)
    norm = y / sf[:, None]
    if not adj:   # task-specified simple statistic, kept as a sensitivity row
        mr = lambda L: mean_ratio_lfc(norm, L, pc_cfg["pseudocount"])
        rows += set_rows(tag, cname, "log2_mean_ratio_norm", False, S, genes, res, lab, mr(lab[None])[0], mr(perms), bnull)
    allg = counts.columns.intersection(res.log2FoldChange.dropna().index)
    obs_all = nb_glm_lfc(counts[allg].values.astype(float), sf, disp[allg].values, X0, lab[None], chunk=1)[0]
    val = {"modality": tag, "contrast": cname, "refit_vs_07c_lfc_pearson": refit_r,
           "sign_agree_nb_vs_deseq2_all_genes": float((np.sign(obs_all) == np.sign(res.log2FoldChange[allg].values)).mean()),
           "pearson_nb_vs_deseq2_all_genes": float(np.corrcoef(obs_all, res.log2FoldChange[allg].values)[0, 1])}
    for r in rows:
        val[f"sign_agree_{r['fast_stat']}_{r['gene_set']}"] = r["sign_agree_fast_vs_deseq2"]
    val["passes"] = all(r["sign_agree_fast_vs_deseq2"] > pc_cfg["min_sign_agreement"] for r in rows if r["primary"])
    if not adj and pc_cfg["n_deseq_check"] > 0:
        checks.extend(deseq_check(tag, cname, pb, md, cfg, bnull, perms[:pc_cfg["n_deseq_check"]],
                                  null_nb[:pc_cfg["n_deseq_check"]], genes, s07c, S))
    return rows, [val]


def deseq_check(tag, cname, pb, md, cfg, bnull, perms, null_fast, genes, s07c, S):
    """Refit the first permutations with full DESeq2 (07c deseq_contrast) and compare with the fast null."""
    ref = S["disc_dE_fdr10"][0]; sg0 = bnull["disc_dE_fdr10"]["genes"]
    out = []
    for i, L in enumerate(perms):
        m = md.copy(); m["group"] = np.where(L, "DKD", "HKD")
        r = s07c.deseq_contrast(pb.loc[m.index], m, cfg, "DKD", "HKD", [])
        sg = [g for g in sg0 if pd.notna(r.log2FoldChange.get(g))]
        d = r.log2FoldChange[sg].values; f = null_fast[i, genes.get_indexer(sg)]; rr = ref[sg].values
        out.append({"modality": tag, "contrast": cname, "perm": i, "n_genes": len(sg),
                    "frac_deseq2": float((np.sign(rr) == np.sign(d)).mean()),
                    "frac_fast": float((np.sign(rr) == np.sign(f)).mean()),
                    "rho_deseq2": float(stats.spearmanr(rr, d).correlation),
                    "rho_fast": float(stats.spearmanr(rr, f).correlation),
                    "sign_agree_fast_vs_deseq2": float((np.sign(d) == np.sign(f)).mean())})
        print(f"  deseq2 check {tag} perm {i}: frac_deseq2={out[-1]['frac_deseq2']:.3f} "
              f"frac_fast={out[-1]['frac_fast']:.3f} agree={out[-1]['sign_agree_fast_vs_deseq2']:.3f}", flush=True)
    return out


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed)
    rng = np.random.default_rng(seed)
    out = stage_dir(STAGE); s07c = load_07c()
    disc = pd.read_csv(ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    lin = pd.read_csv(ROOT / "results/03_pseudobulk_de/lineage_SN_only/PT.tsv", sep="\t").set_index("gene")
    gs = pd.read_csv(ROOT / "results/05_lock_programme/gene_sets.tsv", sep="\t")
    sets = {k: v.gene.tolist() for k, v in gs.groupby("set_id")}
    S = gene_sets(disc, lin, sets)
    rows, val, checks, inputs = [], [], [], []
    for tag in ["sn", "sc"]:
        pbf = C07 / f"iPT_pseudobulk_{tag}_all_calls.tsv.gz"; inputs.append(pbf)
        pb = pd.read_csv(pbf, sep="\t", index_col=0); pb.index = pb.index.astype(str)
        meta = participant_meta(cfg, tag)
        for cname in ["DKD_vs_HKD", "DKD_vs_HKD_adj_egfr_lower"]:
            inputs.append(C07 / f"kpmp_iPT_{cname}_{tag}.tsv")
            r, v = run_contrast(tag, cname, pb, meta, cfg, S, rng, s07c, checks)
            rows += r; val += v
            print(tag, cname, "done", flush=True)
    R = pd.DataFrame(rows); V = pd.DataFrame(val); K = pd.DataFrame(checks)
    R.to_csv(out / "permutation_concordance.tsv", sep="\t", index=False)
    V.to_csv(out / "fast_statistic_validation.tsv", sep="\t", index=False)
    outs = [out / "permutation_concordance.tsv", out / "fast_statistic_validation.tsv"]
    if len(K):
        K.to_csv(out / "deseq2_permutation_check.tsv", sep="\t", index=False)
        outs.append(out / "deseq2_permutation_check.tsv")
        print(K.groupby("modality")[["frac_deseq2", "frac_fast", "rho_deseq2", "rho_fast",
                                     "sign_agree_fast_vs_deseq2"]].mean().round(3).to_string())
    pd.set_option("display.width", 250)
    print(R.round(4).to_string()); print(V.round(4).T.to_string())
    write_provenance(STAGE, inputs + [ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv",
                                      ROOT / "results/03_pseudobulk_de/lineage_SN_only/PT.tsv",
                                      ROOT / "results/05_lock_programme/gene_sets.tsv",
                                      ROOT / "results/07b_validate_kpmp/participants_sn.tsv",
                                      ROOT / "results/07b_validate_kpmp/participants_sc.tsv"],
                     outs, seed, extra={"n_perm": cfg["permutation_concordance"]["n_perm"],
                                        "fast_stats": sorted(R.fast_stat.unique().tolist())})


if __name__ == "__main__":
    main()
