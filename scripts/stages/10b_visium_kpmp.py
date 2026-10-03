"""10b_visium_kpmp: independent spatial test of the locked programmes on KPMP Visium
(fresh-frozen whole-transcriptome sections; participant is the unit of analysis).

Rules (fixed in 10_visium_common / config `visium` before any KPMP spot was scored):
- spots under tissue; PT-rich spot = top pt_rich_top_frac of a section by the PT-lineage marker
  score (10_visium_common.pt_lineage_markers; NNLS proportions are kept as descriptive columns only,
  being mis-calibrated across platforms). Spot resolution cannot separate iPT from healthy PT (10a
  positive control), so injury indices enter as covariates;
- participant pseudobulk = summed counts of PT-rich spots over all of the participant's sections
  (participants with < min_spots_participant PT-rich spots dropped); per-gene z across participants;
  set score = mean z of measured genes (>= min_genes);
- contrasts DKD vs HKD / HKD_withDM / DM_noCKD / Reference (Hedges g, Welch, expression-matched
  random-set null) and OLS adjusted for the shared-injury score and the injury index;
- within-section spot association (DKDiPT_up vs injury index across PT-rich spots), averaged per
  participant, compared across groups (independent test of the 10a observation).
Outputs results/10b_visium_kpmp/ (spot_scores.tsv.gz for the morphology stage)."""
from __future__ import annotations
import glob, importlib, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from stats import compare, matched_random_null
from annotate import zmatrix
vc = importlib.import_module("10_visium_common")

STAGE = "10b_visium_kpmp"
BASE = ROOT / "data/processed/kpmp_visium"
CONTRASTS = [("DKD", "HKD"), ("DKD", "HKD_withDM"), ("DKD", "DM_noCKD"), ("HKD_withDM", "HKD"), ("DKD", "Reference")]


def samples() -> pd.DataFrame:
    rows = []
    tifs = glob.glob(str(ROOT / "data/raw/kpmp/Spatial_Transcriptomics/*/*.tif"))
    for h5 in sorted(glob.glob(str(BASE / "*/*/**/filtered_feature_bc_matrix.h5"), recursive=True)):
        p = Path(h5)
        part = p.relative_to(BASE).parts[0]; samp = p.relative_to(BASE).parts[1]
        pos = next(iter(glob.glob(str(p.parent / "spatial/tissue_positions*.csv"))), None)
        sf = next(iter(glob.glob(str(p.parent / "spatial/scalefactors_json.json"))), None)
        tif = [t for t in tifs if samp in Path(t).name]
        rows.append({"participant": part, "sample": samp, "h5": h5, "positions": pos, "scalefactors": sf,
                     "fullres_tif": tif[0] if tif else None,
                     "hires_png": str(p.parent / "spatial/tissue_hires_image.png")})
    return pd.DataFrame(rows)


def read_sample(r) -> sc.AnnData:
    a = sc.read_10x_h5(r.h5); a.var_names_make_unique()
    pos = pd.read_csv(r.positions, header=0 if r.positions.endswith("tissue_positions.csv") else None)
    pos.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_fullres", "pxl_col_fullres"]
    pos = pos.set_index("barcode")
    a = a[a.obs_names.isin(pos.index[pos.in_tissue == 1])].copy()
    a.obs = a.obs.join(pos[["array_row", "array_col", "pxl_row_fullres", "pxl_col_fullres"]])
    a.obs["participant"] = r.participant; a.obs["sample"] = r.sample
    return a


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); V = cfg["visium"]
    out = stage_dir(STAGE)
    groups = importlib.import_module("00b_download_kpmp").participant_groups().set_index("Participant ID")
    S = samples()
    S["group"] = S.participant.map(groups.group)
    S.to_csv(out / "sample_table.tsv", sep="\t", index=False)
    man = pd.read_csv(ROOT / "data/raw/kpmp/download_manifest_visium_tif.tsv", sep="\t")
    missing = man[~man.file_name.apply(lambda f: any(Path(t).name == f for t in S.fullres_tif.dropna()))]
    missing[["redcap_id", "file_name", "file_size"]].to_csv(out / "fullres_tif_missing.tsv", sep="\t", index=False)

    sets = {k: v for k, v in vc.locked_sets().items() if k in vc.SETS_USED}
    mk = pd.read_csv(ROOT / "results/10a_visium_gse211785/ipt_markers.tsv", sep="\t")
    sets["iPT_marker"] = mk[mk.set == "iPT_marker"].gene.tolist()
    sets["healthyPT_marker"] = mk[mk.set == "healthyPT_marker"].gene.tolist()
    a0 = read_sample(S.iloc[0])
    P, _ = vc.reference_signatures(cfg, V["deconv_markers"], V["deconv_min_cp10k"])
    P = P.loc[P.index.intersection(a0.var_names)]
    P, dmk = vc.reference_signatures_from(P, V["deconv_markers"], V["deconv_min_cp10k"])
    ptm = vc.pt_lineage_markers(P, V["deconv_markers"], V["deconv_min_cp10k"])
    pd.Series(ptm, name="gene").to_csv(out / "pt_lineage_markers.tsv", sep="\t", index=False)
    dmk.to_csv(out / "deconvolution_markers.tsv", sep="\t", index=False)

    spots, pbs, cov_rows = [], {}, []
    for r in S.itertuples():
        a = read_sample(r)
        W, _ = vc.nnls_deconvolve(a.X, a.var_names, P, dmk); W.index = a.obs_names
        Z, _ = vc.lognorm_z(a.X, a.var_names); Z.index = a.obs_names
        sc_, cov = vc.set_scores(Z, sets, V["min_genes"])
        cov_rows.append(cov.assign(sample=r.sample))
        d = a.obs[["participant", "sample", "pxl_row_fullres", "pxl_col_fullres"]].join(sc_)
        d["PT_share"] = (W.PT_healthy + W.iPT).values
        for c in W.columns:
            d[f"prop_{c}"] = W[c].values
        d["PT_lineage_score"] = Z[[g for g in ptm if g in Z.columns]].mean(1).values
        d["PT_rich"] = d.PT_lineage_score >= d.PT_lineage_score.quantile(1 - V["pt_rich_top_frac"])
        d["injury_index"] = d.iPT_marker - d.healthyPT_marker
        d["group"] = r.group; d["fullres_tif_available"] = r.fullres_tif is not None
        spots.append(d)
        m = d.PT_rich.values
        if m.any():
            v = np.asarray(a.X[m].sum(0)).ravel()
            pbs[r.participant] = pbs.get(r.participant, 0) + pd.Series(v, index=a.var_names)
    SP = pd.concat(spots); SP.index.name = "barcode"
    SP.to_csv(out / "spot_scores.tsv.gz", sep="\t")
    C = pd.concat(cov_rows); C.groupby("set")[["n_set", "n_measured"]].median().to_csv(out / "set_coverage.tsv", sep="\t")

    # --- participant level ---
    npt = SP.groupby("participant").PT_rich.sum()
    keep = npt[npt >= V["min_spots_participant"]].index
    PB = pd.DataFrame(pbs).T.loc[keep].fillna(0)
    z, me = zmatrix(PB)
    scores = pd.DataFrame({k: z[[g for g in v if g in z.columns]].mean(1) for k, v in sets.items()})
    meas = {k: len([g for g in v if g in z.columns]) for k, v in sets.items()}
    part = pd.DataFrame({"group": groups.group.reindex(scores.index), "n_PT_rich_spots": npt.reindex(scores.index),
                         "n_sections": S.groupby("participant").size().reindex(scores.index),
                         "PT_rich_frac": SP.groupby("participant").PT_rich.mean().reindex(scores.index)})
    # within-section association, averaged per participant
    assoc = []
    for (p, s_), d in SP[SP.PT_rich].groupby(["participant", "sample"]):
        if len(d) >= V["min_spots_section"]:
            assoc.append({"participant": p, "sample": s_, "n": len(d),
                          "r_injury": stats.spearmanr(d.DKDiPT_up, d.injury_index).correlation,
                          "partial_r_injury_given_SharedInjury_up": vc.partial_r(d.DKDiPT_up, d.injury_index, d.SharedInjury_up)})
    AS = pd.DataFrame(assoc); AS.to_csv(out / "spot_association_by_section.tsv", sep="\t", index=False)
    part = part.join(AS.groupby("participant")[["r_injury", "partial_r_injury_given_SharedInjury_up"]].mean())
    part = part.join(scores)
    part.to_csv(out / "participant_scores.tsv", sep="\t")
    grp = part.group
    rows, adj = [], []
    for a_, b_ in CONTRASTS:
        for k in list(sets) + ["r_injury", "partial_r_injury_given_SharedInjury_up", "PT_rich_frac"]:
            r = compare(part[k], grp, a_, b_)
            if k in sets and k.startswith(("DKDiPT", "DEonly")) and (a_, b_) in CONTRASTS[:3] and r.get("n_a", 0) >= 3 and r.get("n_b", 0) >= 3:
                r.update(matched_random_null(z, sets[k], grp, a_, b_, me, seed=seed))
            rows.append({"metric": k, "contrast": f"{a_}-{b_}", "n_measured": meas.get(k), **r})
        d = part[grp.isin([a_, b_])].copy(); d["case"] = (d.group == a_).astype(int)
        if d.case.nunique() == 2 and len(d) >= 8:
            for k in ["DKDiPT_up", "DKDiPT_down", "DEonly_up"]:
                if meas.get(k, 0) < V["min_genes"]:
                    continue
                for form in [f"{k} ~ case + SharedInjury_up", f"{k} ~ case + SharedInjury_up + injury_index"]:
                    dd = d.assign(injury_index=d.iPT_marker - d.healthyPT_marker)
                    m = smf.ols(form, dd).fit(cov_type="HC3")
                    adj.append({"set": k, "contrast": f"{a_}-{b_}", "model": form, "n": len(dd),
                                "beta_case": m.params["case"], "se": m.bse["case"], "p": m.pvalues["case"]})
    T = pd.DataFrame(rows); T.to_csv(out / "participant_tests.tsv", sep="\t", index=False)
    # transcriptome-wide sign concordance with the discovery SN iPT dE (as 07c for snRNA);
    # effect = difference of group means of log2 CPM of the PT-rich pseudobulk (participant unit)
    disc = pd.read_csv(ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    cpm = np.log2(PB.div(PB.sum(1), axis=0) * 1e6 + 1)
    cc = []
    for a_, b_ in [("DKD", "HKD"), ("DKD", "HKD_withDM"), ("HKD_withDM", "HKD")]:
        diff = cpm[grp == a_].mean() - cpm[grp == b_].mean()
        tt = stats.ttest_ind(cpm[grp == a_], cpm[grp == b_], equal_var=False)
        j = disc.join(pd.DataFrame({"vis_diff": diff, "vis_t": tt.statistic}, index=cpm.columns), how="inner")
        j = j[(cpm[j.index] > 1).mean() >= 0.5]
        for name, m in [("disc_dE_fdr10", j.dE_padj < 0.10), ("all_tested", j.dE_padj.notna())]:
            jj = j[m & j.vis_t.notna()]
            ag = int((np.sign(jj.dE_lfc) == np.sign(jj.vis_diff)).sum())
            cc.append({"contrast": f"{a_}-{b_}", "gene_set": name, "n_genes": len(jj), "n_sign_agree": ag,
                       "frac_agree": ag / max(1, len(jj)),
                       "binom_p": stats.binomtest(ag, len(jj), 0.5, alternative="greater").pvalue if len(jj) else np.nan,
                       "spearman": stats.spearmanr(jj.dE_lfc, jj.vis_t).correlation if len(jj) > 2 else np.nan})
    CC = pd.DataFrame(cc); CC.to_csv(out / "transcriptome_concordance.tsv", sep="\t", index=False)
    print(CC.round(4).to_string())
    A = pd.DataFrame(adj); A.to_csv(out / "injury_adjusted.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(S.groupby("group").participant.nunique(), part.group.value_counts())
    print(pd.Series(meas))
    print(T[T.metric.isin(["DKDiPT_up", "DKDiPT_down", "DEonly_up", "SharedInjury_up", "r_injury",
                           "partial_r_injury_given_SharedInjury_up"])][
        ["metric", "contrast", "n_a", "n_b", "hedges_g", "welch_p", "p_emp_greater"]].round(3).to_string())
    print(A.round(3).to_string())
    write_provenance(STAGE, [*S.h5.map(Path), ROOT / "results/05_lock_programme/gene_sets.tsv",
                             ROOT / "results/10a_visium_gse211785/ipt_markers.tsv",
                             ROOT / "data/raw/kpmp" / cfg["kpmp"]["clinical_table"]],
                     sorted(out.glob("*.tsv")) + [out / "spot_scores.tsv.gz"], seed, extra={"visium": V})


if __name__ == "__main__":
    main()
