"""14f_site_sensitivity: KPMP recruitment site (participant-ID prefix) as a confounder of DKD vs HKD.

Restricted to sites that contributed both DKD and HKD participants in a modality, the iPT
pseudobulk DESeq2 contrast is refitted with ~ sex + site + group and compared with discovery
(same-sign agreement, label-permutation within site), and the marker-gene score is modelled
~ group + site (OLS, HC3)."""
from __future__ import annotations
import importlib, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "14f_site_sensitivity"
R = ROOT / "results"


def deseq_site(pb: pd.DataFrame, md: pd.DataFrame, cfg) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    p = cfg["pseudobulk"]
    counts = pb.loc[md.index]
    counts = counts.loc[:, (counts >= p["min_counts_gene"]).sum(0) >= p["min_samples_gene"]].round().astype(int)
    md = md.copy(); md["site"] = "s" + md.site.astype(str)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(counts=counts, metadata=md[["sex", "site", "group"]], design_factors=["sex", "site", "group"],
                           ref_level=["group", "HKD"], n_cpus=16, quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["group", "DKD", "HKD"], quiet=True); st.summary()
    return st.results_df


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); c = cfg["site_sensitivity"]
    out = stage_dir(STAGE)
    disc = pd.read_csv(R / "03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    genes = disc.index[disc.dE_padj < 0.10]
    rows, score_rows, inputs = [], [], []
    for mod in ["sn", "sc"]:
        f_pb = R / f"07c_kpmp_concordance/iPT_pseudobulk_{mod}_all_calls.tsv.gz"
        f_sc = R / f"07b_validate_kpmp/scores_{mod}.tsv"
        inputs += [f_pb, f_sc]
        pb = pd.read_csv(f_pb, sep="\t", index_col=0)
        sc = pd.read_csv(f_sc, sep="\t", index_col=0)
        meta = pd.read_csv(R / f"07b_validate_kpmp/participants_{mod}.tsv", sep="\t", index_col=0)
        q = pd.read_csv(ROOT / "data/raw/kpmp" / cfg["kpmp"]["clinical_table"]).set_index("Participant ID")
        md = meta.loc[meta.index.intersection(pb.index)]
        md = md[md.group.isin(["DKD", "HKD"])].copy()
        md["sex"] = q.Sex.reindex(md.index).where(lambda s: s.isin(["Male", "Female"]), "Unknown").values
        md["site"] = md.index.str.split("-").str[0]
        both = md.groupby("site").group.nunique()
        keep_sites = both[both == 2].index
        mdk = md[md.site.isin(keep_sites)]
        res = deseq_site(pb, mdk, cfg)
        j = disc.loc[genes].join(res[["log2FoldChange"]], how="inner").dropna()
        agree = float((np.sign(j.dE_lfc) == np.sign(j.log2FoldChange)).mean())
        # label permutation within site (fast statistic: sign of the site-adjusted mean difference of log CPM)
        cpm = np.log2(pb.loc[mdk.index, j.index].div(pb.loc[mdk.index].sum(1), axis=0) * 1e6 + 1)
        resid = cpm - cpm.groupby(mdk.site).transform("mean")
        rng = np.random.default_rng(seed)

        def frac(lab):
            d = resid[lab == "DKD"].mean() - resid[lab == "HKD"].mean()
            return float((np.sign(d) == np.sign(j.dE_lfc)).mean())
        obs_fast = frac(mdk.group.values)
        null = []
        for _ in range(c["n_perm"]):
            lab = mdk.group.copy()
            for s in keep_sites:
                m = (mdk.site == s).values
                lab.values[m] = rng.permutation(lab.values[m])
            null.append(frac(lab.values))
        null = np.array(null)
        rows.append({"modality": mod, "sites_kept": ",".join(keep_sites), "n_DKD": int((mdk.group == "DKD").sum()),
                     "n_HKD": int((mdk.group == "HKD").sum()), "n_genes": len(j), "frac_agree_deseq_site": agree,
                     "frac_agree_fast": obs_fast, "perm_null_mean": null.mean(), "perm_null_q95": np.quantile(null, 0.95),
                     "perm_p": (np.sum(null >= obs_fast) + 1) / (len(null) + 1)})
        d = sc.loc[sc.index.intersection(mdk.index), ["DKDiPT_up"]].join(mdk[["group", "site"]])
        d["dkd"] = (d.group == "DKD").astype(int)
        for name, form, dd in [("all_sites_unadjusted", "DKDiPT_up ~ dkd", sc[sc.group.isin(["DKD", "HKD"])].assign(dkd=lambda x: (x.group == "DKD").astype(int))),
                               ("shared_sites_site_adjusted", "DKDiPT_up ~ dkd + C(site)", d)]:
            m = smf.ols(form, dd).fit(cov_type="HC3")
            score_rows.append({"modality": mod, "model": name, "n": len(dd), "beta_DKD": m.params["dkd"],
                               "se": m.bse["dkd"], "p": m.pvalues["dkd"]})
        res.to_csv(out / f"kpmp_iPT_DKD_vs_HKD_site_adjusted_{mod}.tsv", sep="\t")
    A = pd.DataFrame(rows); A.to_csv(out / "site_adjusted_agreement.tsv", sep="\t", index=False)
    B = pd.DataFrame(score_rows); B.to_csv(out / "site_adjusted_scores.tsv", sep="\t", index=False)
    print(A.round(3).to_string()); print(B.round(4).to_string())
    write_provenance(STAGE, inputs + [R / "03_pseudobulk_de/fine_SN_only/iPT.tsv"],
                     sorted(out.glob("*.tsv")), seed)


if __name__ == "__main__":
    main()
