"""07c_kpmp_concordance: transcriptome-wide concordance complementing the KPMP primary test.

(a) Participant-level iPT pseudobulk DESeq2 in KPMP (DKD vs HKD, ~ sex + group) per modality,
    and transcriptome-wide concordance with the discovery SN-stratum iPT dE.
(b) Per-gene replication of the locked DKDiPT programme.
(c) Programme scores recomputed from high-confidence iPT calls (posterior >= threshold)."""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from annotate import donor_pseudobulk, zmatrix
from stats import compare, matched_random_null

STAGE = "07c_kpmp_concordance"


def deseq_contrast(pb: pd.DataFrame, meta: pd.DataFrame, cfg, a: str, b: str, covars: list[str]) -> pd.DataFrame:
    """DESeq2 on participant pseudobulk: ~ sex + covars + group, contrast a vs b."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    p = cfg["pseudobulk"]
    md = meta.loc[pb.index, ["group", "sex", *covars]].copy()
    md = md[md.group.isin([a, b])].dropna(subset=covars)
    md["sex"] = md.sex.where(md.sex.isin(["Male", "Female"]), "Unknown")
    n = {g: int((md.group == g).sum()) for g in [a, b]}
    md["group"] = md.group.str.replace("_", "")      # pydeseq2 rejects '_' in factor levels
    a, b = a.replace("_", ""), b.replace("_", "")
    for c in covars:   # scale continuous covariates for a stable fit
        md[c] = (md[c] - md[c].mean()) / md[c].std(ddof=1)
    counts = pb.loc[md.index]
    counts = counts.loc[:, (counts >= p["min_counts_gene"]).sum(0) >= p["min_samples_gene"]].round().astype(int)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(counts=counts, metadata=md, design_factors=["sex", *covars, "group"],
                           continuous_factors=covars or None, ref_level=["group", b], n_cpus=16, quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["group", a, b], quiet=True); st.summary()
    r = st.results_df; r.attrs["n"] = n
    return r


def concordance(disc: pd.DataFrame, res: pd.DataFrame, lfc_col: str, masks: dict) -> list[dict]:
    j = disc.join(res[["log2FoldChange", "pvalue", "padj"]].rename(columns=lambda c: f"kpmp_{c}"), how="inner")
    rows = []
    for name, mask in masks.items():
        jj = j[mask.reindex(j.index).fillna(False).astype(bool) & j.kpmp_log2FoldChange.notna()]
        agree = int((np.sign(jj[lfc_col]) == np.sign(jj.kpmp_log2FoldChange)).sum())
        rows.append({"gene_set": name, "n_genes": len(jj), "n_sign_agree": agree, "frac_agree": agree / max(1, len(jj)),
                     "binom_p": stats.binomtest(agree, len(jj), 0.5, alternative="greater").pvalue if len(jj) else np.nan,
                     "spearman_lfc": stats.spearmanr(jj[lfc_col], jj.kpmp_log2FoldChange).correlation if len(jj) > 2 else np.nan,
                     "n_kpmp_padj10_same_sign": int(((jj.kpmp_padj < 0.10) & (np.sign(jj[lfc_col]) == np.sign(jj.kpmp_log2FoldChange))).sum())})
    return rows


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed)
    out = stage_dir(STAGE)
    disc = pd.read_csv(ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    gs = pd.read_csv(ROOT / "results/05_lock_programme/gene_sets.tsv", sep="\t")
    sets = {k: v.gene.tolist() for k, v in gs.groupby("set_id")}
    conc, pergene, hc = [], [], []
    for tag in ["sn", "sc"]:
        q = ad.read_h5ad(ROOT / cfg["paths"]["processed"] / f"KPMP_{tag}_qc_counts.h5ad")
        lab = pd.read_csv(ROOT / f"results/07b_validate_kpmp/cell_labels_{tag}.tsv.gz", sep="\t", index_col=0)
        q.obs = q.obs.join(lab[["label", "max_prob"]])
        meta = q.obs.groupby("participant").agg(group=("group", "first"), sex=("sex", "first"))
        for conf_name, thr in [("all_calls", 0.0), ("high_confidence", cfg["kpmp_concordance"]["min_posterior"])]:
            m = (q.obs.label.eq("iPT") & (q.obs.max_prob >= thr)).values
            n = q.obs[m].groupby("participant").size()
            pb = donor_pseudobulk(q, m, "participant")
            pb = pb.loc[n[n >= cfg["pseudobulk"]["min_cells_per_sample"]].index]
            pb.to_csv(out / f"iPT_pseudobulk_{tag}_{conf_name}.tsv.gz", sep="\t")
            grp = meta.group.reindex(pb.index)
            z, me = zmatrix(pb)
            for k in ["DKDiPT_up", "DEonly_up", "SharedInjury_up"]:
                s = z[[g for g in sets[k] if g in z.columns]].mean(1)
                r = compare(s, grp, "DKD", "HKD")
                r.update(matched_random_null(z, sets[k], grp, "DKD", "HKD", me, seed=seed))
                hc.append({"modality": tag, "cells": conf_name, "set": k, **r})
            if conf_name != "all_calls":
                continue
            desc = pd.read_csv(ROOT / f"results/07b_validate_kpmp/participants_{tag}.tsv", sep="\t", index_col=0)
            meta2 = meta.join(desc[["egfr_lower", "if_pct"]])
            lin = pd.read_csv(ROOT / "results/03_pseudobulk_de/lineage_SN_only/PT.tsv", sep="\t").set_index("gene")
            shared = set(sets["SharedInjury_up"] + sets["SharedInjury_down"])
            masks = {"disc_dE_fdr10": disc.dE_padj < 0.10,
                     "disc_dE_fdr10_and_dD_concordant": (disc.dE_padj < 0.10) & (np.sign(disc.dE_lfc) == np.sign(disc.dD_lfc)),
                     "all_tested": disc.dE_padj.notna()}
            runs = [("DKD", "HKD", []), ("DKD", "HKD", ["egfr_lower"]), ("DKD", "HKD", ["if_pct"]),
                    ("DKD", "HKD_withDM", []), ("HKD_withDM", "HKD", [])]
            for a, b, cv in runs:
                res = deseq_contrast(pb, meta2, cfg, a, b, cv)
                tagname = f"{a}_vs_{b}" + ("_adj_" + "_".join(cv) if cv else "")
                res.to_csv(out / f"kpmp_iPT_{tagname}_{tag}.tsv", sep="\t")
                for r in concordance(disc, res, "dE_lfc", masks):
                    conc.append({"modality": tag, "contrast": tagname, "n_a": res.attrs["n"][a],
                                 "n_b": res.attrs["n"][b], **r})
                if (a, b, cv) == ("DKD", "HKD", []):
                    # negative control: shared-injury genes should not separate DKD from HKD
                    sh = lin.loc[lin.index.intersection(list(shared))].assign(dE_lfc=lambda d: (d.dD_lfc + d.dN_lfc) / 2)
                    for r in concordance(sh, res, "dE_lfc", {"shared_injury_genes": pd.Series(True, index=sh.index)}):
                        conc.append({"modality": tag, "contrast": tagname, "n_a": res.attrs["n"][a],
                                     "n_b": res.attrs["n"][b], **r})
                    main_res = res
            res = main_res
            for d in ["up", "down"]:
                for g in sets[f"DKDiPT_{d}"]:
                    pergene.append({"modality": tag, "set": f"DKDiPT_{d}", "gene": g,
                                    "disc_dE_lfc": disc.dE_lfc.get(g, np.nan),
                                    "kpmp_lfc": res.log2FoldChange.get(g, np.nan), "kpmp_p": res.pvalue.get(g, np.nan)})
        print(tag, "done", flush=True)
    C = pd.DataFrame(conc); C.to_csv(out / "transcriptome_concordance.tsv", sep="\t", index=False)
    G = pd.DataFrame(pergene); G["same_sign"] = np.sign(G.disc_dE_lfc) == np.sign(G.kpmp_lfc)
    G.to_csv(out / "per_gene_replication.tsv", sep="\t", index=False)
    H = pd.DataFrame(hc); H.to_csv(out / "scores_by_call_confidence.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(C.round(4).to_string()); print(G.round(3).to_string())
    print(H[["modality", "cells", "set", "n_a", "n_b", "hedges_g", "welch_p", "p_emp_greater"]].round(4).to_string())
    write_provenance(STAGE, [ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv",
                             ROOT / "results/05_lock_programme/gene_sets.tsv"],
                     sorted(out.glob("*.tsv")), seed)


if __name__ == "__main__":
    main()
