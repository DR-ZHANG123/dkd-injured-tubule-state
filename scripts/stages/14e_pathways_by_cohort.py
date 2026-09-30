"""14e_pathways_by_cohort: preranked GSEA of the DKD-vs-HKD iPT shift in each cohort separately.

Rankings (Wald statistics, iPT pseudobulk DKD vs HKD):
  discovery_SN          results/03_pseudobulk_de/fine_SN_only/iPT.tsv (dE_lfc / dE_se)
  discovery_SN_ribo_adj results/12a_technical_robustness/discovery_iPT_DKD_vs_HKD_adj_frac_ribo.tsv
  KPMP_sn               results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv
  KPMP_sn_ribo_adj      results/12a_technical_robustness/kpmp_sn_iPT_DKD_vs_HKD_adj_frac_ribo.tsv
Gene filters: 'excl_technical' (ribosomal / mitochondrial / translation genes removed, as in 12a;
Hallmark, Reactome, GO BP) and 'full' (all genes; Hallmark and Reactome, so translation terms stay
testable). Each ranking uses only its own cohort's tested genes. A pathway is called consistent
for a discovery x KPMP pair when NES has the same sign and FDR < consistency_fdr in both."""
from __future__ import annotations
import importlib, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "14e_pathways_by_cohort"
m09 = importlib.import_module("09_signature_biology")
m12a = importlib.import_module("12a_technical_robustness")
REF = ROOT / "data/raw/reference"
SOURCES = {
    "discovery_SN": ("results/03_pseudobulk_de/fine_SN_only/iPT.tsv", "gene", "dE_lfc", "dE_se"),
    "discovery_SN_ribo_adj": ("results/12a_technical_robustness/discovery_iPT_DKD_vs_HKD_adj_frac_ribo.tsv", 0, "dE_lfc", "dE_se"),
    "KPMP_sn": ("results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv", 0, "log2FoldChange", "lfcSE"),
    "KPMP_sn_ribo_adj": ("results/12a_technical_robustness/kpmp_sn_iPT_DKD_vs_HKD_adj_frac_ribo.tsv", 0, "log2FoldChange", "lfcSE"),
}
PAIRS = [("discovery_SN", "KPMP_sn"), ("discovery_SN_ribo_adj", "KPMP_sn_ribo_adj"),
         ("discovery_SN", "KPMP_sn_ribo_adj"), ("discovery_SN_ribo_adj", "KPMP_sn")]


def rankings(ex: set) -> dict[tuple[str, str], pd.Series]:
    out = {}
    for name, (path, idx, lfc, se) in SOURCES.items():
        df = pd.read_csv(ROOT / path, sep="\t", index_col=idx)
        z = m09.wald(df, lfc, se).dropna()
        z = z[~z.index.duplicated()]
        out[(name, "full")] = z
        out[(name, "excl_technical")] = z[~m12a.is_technical(z.index, ex)]
    return out


def run_gsea(R: dict, cfg, seed) -> pd.DataFrame:
    full_libs = cfg["pathways_by_cohort"]["full_gene_libraries"]
    rows = []
    for (name, filt), z in R.items():
        for lib, f in m09.GMTS.items():
            if filt == "full" and lib not in full_libs:
                continue
            r = m09.prerank(z, REF / f, seed, cfg)
            r.insert(0, "library", lib); r.insert(0, "gene_filter", filt); r.insert(0, "ranking", name)
            r["n_ranked_genes"] = len(z)
            rows.append(r)
        print("gsea", name, filt, flush=True)
    return pd.concat(rows, ignore_index=True)


def consistency(G: pd.DataFrame, thr: float) -> pd.DataFrame:
    rows = []
    for filt, d in G.groupby("gene_filter"):
        w = d.pivot_table(index=["library", "term"], columns="ranking", values=["NES", "fdr"])
        for a, b in PAIRS:
            if ("NES", a) not in w.columns or ("NES", b) not in w.columns:
                continue
            x = pd.DataFrame({"NES_disc": w[("NES", a)], "fdr_disc": w[("fdr", a)],
                              "NES_kpmp": w[("NES", b)], "fdr_kpmp": w[("fdr", b)]}).dropna()
            x["same_sign"] = np.sign(x.NES_disc) == np.sign(x.NES_kpmp)
            x["consistent"] = x.same_sign & (x.fdr_disc < thr) & (x.fdr_kpmp < thr)
            x["direction"] = np.where(x.NES_disc > 0, "up_in_DKD", "down_in_DKD")
            x = x.reset_index().assign(gene_filter=filt, discovery=a, kpmp=b)
            rows.append(x)
    return pd.concat(rows, ignore_index=True)


def nes_agreement(Cn: pd.DataFrame) -> pd.DataFrame:
    """Transcriptome-wide pathway-level agreement: sign agreement and Spearman r of NES."""
    from scipy import stats
    rows = []
    for (filt, a, b, lib), d in Cn.groupby(["gene_filter", "discovery", "kpmp", "library"]):
        sig = d[(d.fdr_disc < 0.25)]
        rows.append({"gene_filter": filt, "discovery": a, "kpmp": b, "library": lib, "n_terms": len(d),
                     "spearman_NES": stats.spearmanr(d.NES_disc, d.NES_kpmp).correlation,
                     "n_disc_fdr25": len(sig), "frac_sign_agree_disc_fdr25": float(sig.same_sign.mean()) if len(sig) else np.nan,
                     "n_consistent": int(d.consistent.sum()),
                     "n_consistent_up": int((d.consistent & (d.NES_disc > 0)).sum()),
                     "n_consistent_down": int((d.consistent & (d.NES_disc < 0)).sum())})
    return pd.DataFrame(rows)


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); P = cfg["pathways_by_cohort"]
    out = stage_dir(STAGE)
    ex = m12a.exclusion_genes(cfg)
    R = rankings(ex)
    pd.DataFrame([{"ranking": n, "gene_filter": f, "n_genes": len(z)} for (n, f), z in R.items()]).to_csv(
        out / "ranking_sizes.tsv", sep="\t", index=False)
    G = run_gsea(R, cfg, seed); G.to_csv(out / "gsea_by_cohort.tsv", sep="\t", index=False)
    Cn = consistency(G, P["consistency_fdr"]); Cn.to_csv(out / "pathway_consistency_all_terms.tsv", sep="\t", index=False)
    Cn[Cn.consistent].sort_values(["gene_filter", "discovery", "kpmp", "NES_disc"]).to_csv(
        out / "consistent_pathways.tsv", sep="\t", index=False)
    Ag = nes_agreement(Cn); Ag.to_csv(out / "pathway_level_agreement.tsv", sep="\t", index=False)
    K = G[G.term.isin(P["key_terms"])]
    missing = sorted(set(P["key_terms"]) - set(K.term))
    Kw = K.pivot_table(index="term", columns=["gene_filter", "ranking"], values=["NES", "fdr"])
    Kw.columns = [f"{v}|{f}|{r}" for v, f, r in Kw.columns]
    Kw.to_csv(out / "key_pathways_by_cohort.tsv", sep="\t")
    Kc = Cn[Cn.term.isin(P["key_terms"])].sort_values(["term", "gene_filter", "discovery", "kpmp"])
    Kc.to_csv(out / "key_pathways_consistency.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30, "display.max_colwidth", 50)
    print("key terms not testable (size filter / absent):", missing)
    print(K.pivot_table(index="term", columns=["gene_filter", "ranking"], values="NES").round(2).to_string())
    print(K.pivot_table(index="term", columns=["gene_filter", "ranking"], values="fdr").round(3).to_string())
    print(Kc[["term", "gene_filter", "discovery", "kpmp", "NES_disc", "fdr_disc", "NES_kpmp", "fdr_kpmp", "consistent"]].round(3).to_string())
    print(Ag.round(3).to_string())
    write_provenance(STAGE, [ROOT / p for p, *_ in SOURCES.values()] + [REF / f for f in m09.GMTS.values()],
                     sorted(out.glob("*.tsv")), seed,
                     extra={"key_terms_not_testable": missing, "n_technical_genes_excluded_list": len(ex)})


if __name__ == "__main__":
    main()
