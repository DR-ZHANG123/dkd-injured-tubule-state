"""14d_contamination_controls: is the DKD-vs-HKD iPT shift leukocyte / ambient RNA carried into
iPT nuclei rather than a tubular state?

(1) Per sample (KPMP sn, KPMP sc participants; discovery SN libraries) contamination indices of the
    iPT pseudobulk: leukocyte-gene fraction (canonical and strict / no MHC-II), stromal-gene fraction,
    immunoglobulin fraction, non-PT marker fraction, ratio estimates of the immune / stromal share,
    and the correlation of the iPT profile with the same sample's non-PT (immune + stroma) pseudobulk.
    DKD vs HKD on every index.
(2) Programme scores (mean z, as in 07b/07c) refitted as score ~ DKD + index (OLS, HC3), and
    re-scored after residualising every gene on the index (conservative: removes any group signal
    collinear with the index).
(3) Expression specificity of the programme genes: participant-mean CPM in iPT / healthy PT / immune /
    stroma (KPMP sn; discovery SN for comparison), count share by compartment, and the CPM expected
    in iPT from contamination (ratio estimate x compartment CPM) versus the observed CPM.
(4) Transcriptome concordance (discovery SN iPT dE vs KPMP sn DKD vs HKD) after removing genes whose
    immune or stromal CPM exceeds a ratio x iPT CPM; HGNC locus_group and GO BP synaptic flags.
    Gene length is not assessed: no verified local gene-length annotation."""
from __future__ import annotations
import importlib, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from annotate import zmatrix
from stats import compare, hedges_g
from contamination import (compartment_map, grouped_pseudobulk, cpm, nonpt_markers, contamination_indices,
                           expected_contamination_share, residualize)

STAGE = "14d_contamination_controls"
m07c = importlib.import_module("07c_kpmp_concordance")
COMPS = ["iPT", "PT_healthy", "immune", "stroma"]
SETS = ["DKDiPT_up", "DKDiPT_down", "SharedInjury_up", "DEonly_up"]
INDICES = ["frac_leukocyte_genes", "frac_leukocyte_strict", "frac_stroma_strict", "frac_ig_genes",
           "frac_nonpt_markers", "c_immune", "c_stroma", "r_ambient_all", "r_ambient_markers"]
JOINT = ["frac_leukocyte_strict", "c_stroma", "frac_ig_genes", "r_ambient_markers"]


def split_compartments(pb: pd.DataFrame, n: pd.Series) -> tuple[dict, dict]:
    """'sample|compartment' rows -> {compartment: samples x genes}, {compartment: cells per sample}."""
    s = pb.index.str.rsplit("|", n=1).str[0]; c = pb.index.str.rsplit("|", n=1).str[1]
    out, nn = {}, {}
    for comp in COMPS:
        m = np.asarray(c == comp)
        out[comp] = pb[m].set_axis(s[m]); nn[comp] = n[m].set_axis(s[m])
    return out, nn


def load_kpmp(tag: str, cfg) -> tuple[dict, dict, pd.DataFrame]:
    q = ad.read_h5ad(ROOT / cfg["paths"]["processed"] / f"KPMP_{tag}_qc_counts.h5ad")
    lab = pd.read_csv(ROOT / f"results/07b_validate_kpmp/cell_labels_{tag}.tsv.gz", sep="\t", index_col=0)
    comp = lab.label.reindex(q.obs_names).map(compartment_map(cfg)).fillna("other")
    keys = (q.obs.participant.astype(str) + "|" + comp.values).values
    pb, n = grouped_pseudobulk(q.X, keys, q.var_names)
    meta = q.obs.groupby("participant", observed=True).agg(group=("group", "first"), sex=("sex", "first"))
    meta.index = meta.index.astype(str)
    del q
    comps, ncell = split_compartments(pb, n)
    return comps, ncell, meta


def load_discovery(cfg) -> tuple[dict, dict, pd.DataFrame]:
    """Discovery SN libraries from the stage-02 fine pseudobulk (raw counts summed per library x cell type)."""
    q = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    pbf = ad.read_h5ad(ROOT / cfg["paths"]["processed"] / "pb_fine.h5ad")
    o = pbf.obs
    sel = ((o.library.map(q.batch) == "SN") & o.library.map(q.qc_pass).astype(bool)).values
    comp = o.celltype.astype(str).map(compartment_map(cfg)).fillna("other")[sel]
    X = pbf.X[sel]
    keys = (o.library.astype(str)[sel] + "|" + comp).values
    pb, _ = grouped_pseudobulk(X, keys, pbf.var_names)
    ncell = (o[sel].assign(k=keys).groupby("k").n_cells.sum())
    meta = o[sel].groupby("library", observed=True).agg(group=("group", "first"), sex=("sex", "first"))
    meta.index = meta.index.astype(str)
    comps, nn = split_compartments(pb, ncell.reindex(pb.index))
    return comps, nn, meta


def cohort_block(name: str, comps: dict, ncell: dict, meta: pd.DataFrame, sets: dict, cfg) -> dict:
    C = cfg["contamination_controls"]; minc = cfg["pseudobulk"]["min_cells_per_sample"]
    pooled = pd.DataFrame({c: comps[c].sum(0) for c in COMPS}).T
    pooled_cpm = cpm(pooled)
    prog = set(sets["DKDiPT_up"] + sets["DKDiPT_down"])
    markers = nonpt_markers(pooled, prog, C["nonpt_marker_ratio"], C["nonpt_marker_min_cpm"])
    keep = ncell["iPT"][ncell["iPT"] >= minc].index
    ipt = comps["iPT"].loc[keep]
    nonpt = comps["immune"].add(comps["stroma"], fill_value=0)
    n_nonpt = ncell["immune"].add(ncell["stroma"], fill_value=0)
    idx = contamination_indices(ipt, nonpt, n_nonpt, pooled_cpm, C, markers, C["min_cells_compartment"])
    idx["group"] = meta.group.reindex(idx.index).values
    z, _ = zmatrix(ipt)
    for k in SETS:
        idx[k] = z[[g for g in sets[k] if g in z.columns]].mean(1)
    idx.insert(0, "cohort", name)
    return {"idx": idx, "z": z, "ipt": ipt, "pooled_cpm": pooled_cpm, "markers": markers, "comps": comps,
            "ncell": ncell, "meta": meta}


def index_tests(B: dict) -> list[dict]:
    d = B["idx"]
    return [{"cohort": d.cohort.iloc[0], "index": i, **compare(d[i], d.group, "DKD", "HKD")} for i in INDICES]


def adjusted_models(B: dict, sets: dict) -> tuple[list[dict], list[dict]]:
    d = B["idx"][B["idx"].group.isin(["DKD", "HKD"])].copy(); z = B["z"]
    d["dkd"] = (d.group == "DKD").astype(int); cohort = d.cohort.iloc[0]
    rows, genes = [], []
    for k in SETS:
        gg = [g for g in sets[k] if g in z.columns]
        for cov in ["none", *INDICES, "joint"]:
            cv = JOINT if cov == "joint" else ([] if cov == "none" else [cov])
            dd = d.dropna(subset=[k, *cv]).copy()
            if dd.dkd.sum() < 2 or (1 - dd.dkd).sum() < 2:
                continue
            if any(dd[c].std(ddof=1) == 0 for c in cv):   # constant index (e.g. no Ig reads) -> not estimable
                continue
            for c in cv:
                dd[c] = (dd[c] - dd[c].mean()) / dd[c].std(ddof=1)
            m = smf.ols(f"{k} ~ dkd" + "".join(f" + {c}" for c in cv), dd).fit(cov_type="HC3")
            r = {"cohort": cohort, "set": k, "covariate": cov, "n_DKD": int(dd.dkd.sum()), "n_HKD": int((1 - dd.dkd).sum()),
                 "beta_DKD": m.params["dkd"], "se_DKD": m.bse["dkd"], "p_DKD": m.pvalues["dkd"],
                 "beta_cov": m.params[cv[0]] if len(cv) == 1 else np.nan, "p_cov": m.pvalues[cv[0]] if len(cv) == 1 else np.nan}
            if len(cv) == 1:   # conservative: residualise every gene on the index, then re-score
                zr = residualize(z.loc[dd.index, gg], dd[cv[0]])
                s = zr.mean(1); grp = dd.group
                r["g_residualised"] = hedges_g(s[grp == "DKD"], s[grp == "HKD"])
                r["welch_p_residualised"] = compare(s, grp, "DKD", "HKD").get("welch_p", np.nan)
            else:
                r["g_residualised"] = hedges_g(dd.loc[dd.dkd == 1, k], dd.loc[dd.dkd == 0, k]) if cov == "none" else np.nan
            rows.append(r)
            if k.startswith("DKDiPT") and cov in ["none", "frac_leukocyte_strict", "c_stroma", "frac_nonpt_markers",
                                                  "r_ambient_markers", "joint"]:
                for g in gg:
                    dg = dd.assign(y=z.loc[dd.index, g])
                    mg = smf.ols("y ~ dkd" + "".join(f" + {c}" for c in cv), dg).fit(cov_type="HC3")
                    genes.append({"cohort": cohort, "set": k, "gene": g, "covariate": cov,
                                  "beta_DKD": mg.params["dkd"], "p_DKD": mg.pvalues["dkd"]})
    return rows, genes


def participant_mean_cpm(B: dict, samples, minc: int) -> pd.DataFrame:
    """compartment x gene: mean over samples (with >= minc cells in that compartment) of CPM."""
    out = {}
    for c in COMPS:
        s = [x for x in samples if B["ncell"][c].get(x, 0) >= minc]
        out[c] = cpm(B["comps"][c].loc[s]).mean(0)
        out[f"n_{c}"] = pd.Series(len(s), index=out[c].index)
    return pd.DataFrame(out).T


def specificity(B: dict, genes: list[str], minc: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    samples = B["meta"].index[B["meta"].group.isin(["DKD", "HKD"])]
    M = participant_mean_cpm(B, samples, minc)
    allg = M.loc[COMPS].T.astype(float)
    allg["ratio_immune_iPT"] = allg.immune / allg.iPT; allg["ratio_stroma_iPT"] = allg.stroma / allg.iPT
    allg["ratio_max_iPT"] = allg[["ratio_immune_iPT", "ratio_stroma_iPT"]].max(1)
    cnt = pd.DataFrame({c: B["comps"][c].reindex(samples).dropna(how="all").sum(0) for c in COMPS})
    share = cnt.div(cnt.sum(1), axis=0).add_prefix("count_share_")
    sp_ = allg.join(share)
    sp_.insert(0, "cohort", B["idx"].cohort.iloc[0])
    return sp_, sp_.loc[[g for g in genes if g in sp_.index]]


def expected_summary(B: dict, genes: list[str]) -> pd.DataFrame:
    ipt, idx = B["ipt"], B["idx"]
    E = expected_contamination_share(ipt, idx, B["pooled_cpm"], genes)
    E["group"] = E["sample"].map(idx.group)
    E = E[E.group.isin(["DKD", "HKD"])]
    rows = []
    for g, d in E.groupby("gene"):
        a, b = d[d.group == "DKD"], d[d.group == "HKD"]
        dobs = a.obs_cpm.mean() - b.obs_cpm.mean()
        rows.append({"cohort": idx.cohort.iloc[0], "gene": g, "obs_cpm_DKD": a.obs_cpm.mean(), "obs_cpm_HKD": b.obs_cpm.mean(),
                     "median_share_immune_DKD": a.share_immune.median(), "median_share_immune_HKD": b.share_immune.median(),
                     "median_share_stroma_DKD": a.share_stroma.median(), "median_share_stroma_HKD": b.share_stroma.median(),
                     "delta_obs_cpm": dobs,
                     "delta_explained_immune": (a.exp_immune_cpm.mean() - b.exp_immune_cpm.mean()) / dobs if dobs else np.nan,
                     "delta_explained_stroma": (a.exp_stroma_cpm.mean() - b.exp_stroma_cpm.mean()) / dobs if dobs else np.nan})
    return pd.DataFrame(rows)


def gene_flags(cfg) -> tuple[pd.Series, set]:
    h = pd.read_csv(ROOT / "data/raw/reference/hgnc_complete_set.txt", sep="\t", usecols=["symbol", "locus_group"],
                    low_memory=False).drop_duplicates("symbol").set_index("symbol").locus_group
    rx = re.compile(cfg["contamination_controls"]["synapse_term_regex"])
    syn = set()
    for line in open(ROOT / "data/raw/reference/c5.go.bp.v2024.1.Hs.symbols.gmt"):
        f = line.rstrip("\n").split("\t")
        if rx.search(f[0]):
            syn |= set(f[2:])
    return h, syn


def concordance_filtered(spec: dict, cfg, locus: pd.Series, syn: set) -> pd.DataFrame:
    C = cfg["contamination_controls"]
    disc = pd.read_csv(ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    kp = pd.read_csv(ROOT / "results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv", sep="\t", index_col=0)
    f10 = disc.dE_padj < 0.10
    lg = locus.reindex(disc.index)
    masks = {"disc_dE_fdr10": f10}
    for src, S in spec.items():
        r = S.ratio_max_iPT.reindex(disc.index)
        for thr_name, thr in [("", C["specificity_ratio"]), ("_strict", C["specificity_ratio_strict"])]:
            masks[f"fdr10_excl_nonPT_higher{thr_name}_{src}"] = f10 & (r <= thr)
        masks[f"fdr10_only_nonPT_higher_{src}"] = f10 & (r > C["specificity_ratio"])
    masks["fdr10_protein_coding"] = f10 & (lg == "protein-coding gene")
    masks["fdr10_non_coding_RNA"] = f10 & (lg == "non-coding RNA")
    masks["fdr10_non_synaptic"] = f10 & ~disc.index.isin(list(syn))
    masks["fdr10_protein_coding_non_synaptic_excl_nonPT_higher_KPMP_sn"] = (
        masks["fdr10_protein_coding"] & masks["fdr10_non_synaptic"] & masks["fdr10_excl_nonPT_higher_KPMP_sn"])
    rows = m07c.concordance(disc, kp, "dE_lfc", masks)
    base = pd.DataFrame(rows)
    comp = {"all_tested": disc.dE_padj.notna(), "fdr10": f10}
    lgc = [{"gene_set": k, "n": int(m.sum()), **{f"frac_{v.split()[0]}": float((lg[m] == v).mean())
            for v in ["protein-coding gene", "non-coding RNA"]}, "frac_not_in_HGNC": float(lg[m].isna().mean()),
            "frac_synaptic": float(disc.index[m].isin(list(syn)).mean())}
           for k, m in comp.items()]
    return base, pd.DataFrame(lgc)


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); C = cfg["contamination_controls"]
    out = stage_dir(STAGE); minc = C["min_cells_compartment"]
    gs = pd.read_csv(ROOT / "results/05_lock_programme/gene_sets.tsv", sep="\t")
    sets = {k: v.gene.tolist() for k, v in gs.groupby("set_id")}
    prog = sets["DKDiPT_up"] + sets["DKDiPT_down"]
    report_genes = prog + [g for g in C["leukocyte_genes"] + C["stroma_strict"] if g not in prog]
    blocks = {}
    for name, loader in [("discovery_SN", lambda: load_discovery(cfg)), ("KPMP_sn", lambda: load_kpmp("sn", cfg)),
                         ("KPMP_sc", lambda: load_kpmp("sc", cfg))]:
        comps, ncell, meta = loader()
        blocks[name] = cohort_block(name, comps, ncell, meta, sets, cfg)
        print(name, "iPT samples", len(blocks[name]["idx"]), "non-PT markers", len(blocks[name]["markers"]), flush=True)
    I = pd.concat([b["idx"] for b in blocks.values()]); I.to_csv(out / "contamination_indices_by_sample.tsv", sep="\t")
    T = pd.DataFrame([r for b in blocks.values() for r in index_tests(b)])
    T.to_csv(out / "index_DKD_vs_HKD.tsv", sep="\t", index=False)
    rows, genes = [], []
    for b in blocks.values():
        r, g = adjusted_models(b, sets); rows += r; genes += g
    A = pd.DataFrame(rows); A.to_csv(out / "score_models_adjusted.tsv", sep="\t", index=False)
    G = pd.DataFrame(genes); G.to_csv(out / "per_gene_models_adjusted.tsv", sep="\t", index=False)
    locus, syn = gene_flags(cfg)
    spec_all, spec_prog, exp_rows = {}, [], []
    for name in ["KPMP_sn", "discovery_SN", "KPMP_sc"]:
        sa, sp_ = specificity(blocks[name], report_genes, minc)
        spec_all[name] = sa; spec_prog.append(sp_)
        exp_rows.append(expected_summary(blocks[name], prog))
    S = pd.concat(spec_prog); S["set"] = S.index.map(lambda g: next((k for k in ["DKDiPT_up", "DKDiPT_down"] if g in sets[k]), "index_gene"))
    S["locus_group"] = S.index.map(locus); S["synaptic_GO"] = S.index.isin(list(syn))
    S.index.name = "gene"; S.to_csv(out / "marker_specificity.tsv", sep="\t")
    spec_all["KPMP_sn"].to_csv(out / "all_gene_specificity_KPMP_sn.tsv.gz", sep="\t")
    E = pd.concat(exp_rows); E.to_csv(out / "expected_contamination_programme_genes.tsv", sep="\t", index=False)
    Cc, L = concordance_filtered({"KPMP_sn": spec_all["KPMP_sn"], "discovery_SN": spec_all["discovery_SN"]}, cfg, locus, syn)
    Cc.to_csv(out / "concordance_specificity_filtered.tsv", sep="\t", index=False)
    L.to_csv(out / "gene_class_composition.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(T[["cohort", "index", "n_a", "n_b", "mean_a", "mean_b", "hedges_g", "welch_p"]].round(4).to_string())
    print(A[A.set.isin(["DKDiPT_up", "DKDiPT_down"])][["cohort", "set", "covariate", "n_DKD", "n_HKD", "beta_DKD", "p_DKD",
                                                       "p_cov", "g_residualised", "welch_p_residualised"]].round(4).to_string())
    print(S[S.cohort == "KPMP_sn"][["iPT", "PT_healthy", "immune", "stroma", "ratio_immune_iPT", "ratio_stroma_iPT",
                                    "count_share_iPT", "count_share_immune", "locus_group"]].round(3).to_string())
    print(E.round(3).to_string()); print(Cc.round(4).to_string()); print(L.round(3).to_string())
    write_provenance(STAGE, [ROOT / "data/processed/KPMP_sn_qc_counts.h5ad", ROOT / "data/processed/KPMP_sc_qc_counts.h5ad",
                             ROOT / "data/processed/pb_fine.h5ad", ROOT / "results/02c_library_qc/library_qc.tsv",
                             ROOT / "results/07b_validate_kpmp/cell_labels_sn.tsv.gz",
                             ROOT / "results/07b_validate_kpmp/cell_labels_sc.tsv.gz",
                             ROOT / "results/05_lock_programme/gene_sets.tsv",
                             ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv",
                             ROOT / "results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv",
                             ROOT / "data/raw/reference/hgnc_complete_set.txt",
                             ROOT / "data/raw/reference/c5.go.bp.v2024.1.Hs.symbols.gmt"],
                     sorted(out.glob("*.tsv")) + sorted(out.glob("*.tsv.gz")), seed,
                     extra={"n_nonpt_markers": {k: len(b["markers"]) for k, b in blocks.items()},
                            "gene_length": "not assessed: no verified local gene-length annotation"})


if __name__ == "__main__":
    main()
