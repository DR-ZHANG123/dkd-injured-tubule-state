"""12a_technical_robustness: is the DKD-vs-HKD iPT shift (and its ribosome/translation GSEA
signal) a cytoplasmic / ambient-RNA artefact of snRNA?

Per sample (discovery SN iPT pseudobulk; KPMP snRNA iPT pseudobulk): ribosomal-protein,
mitochondrial and cytoplasmic-marker read fractions (cytoplasmic markers = the genes most
enriched in scRNA over snRNA on the PT technical axis of 03b). DKD vs HKD tests on these; the
correlation of signature effects with the technical axis; concordance + preranked GSEA re-run
(i) excluding ribosomal / mitochondrial / translation genes and (ii) with DESeq2 adjusted for the
sample ribosomal fraction."""
from __future__ import annotations
import importlib, sys, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from stats import compare

STAGE = "12a_technical_robustness"
m07c = importlib.import_module("07c_kpmp_concordance")
m09 = importlib.import_module("09_signature_biology")
REF = ROOT / "data/raw/reference"


def read_gmt(path: Path) -> dict[str, list[str]]:
    out = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        out[f[0]] = f[2:]
    return out


def exclusion_genes(cfg) -> set:
    gm = {**read_gmt(REF / "c5.go.bp.v2024.1.Hs.symbols.gmt"), **read_gmt(REF / "c2.cp.reactome.v2024.1.Hs.symbols.gmt")}
    ex = set()
    for t in cfg["robustness"]["exclude_terms"]:
        ex |= set(gm[t])
    return ex


def is_technical(genes: pd.Index, ex: set) -> np.ndarray:
    genes = pd.Index(genes)
    return np.asarray(genes.str.match(r"^(RP[LS]\d|MRP[LS]\d|MT-)")) | np.asarray(genes.isin(ex))


def tech_fractions(pb: pd.DataFrame, cyto: list[str]) -> pd.DataFrame:
    g = pb.columns; tot = pb.sum(1)
    return pd.DataFrame({"frac_ribo": pb.loc[:, g.str.match(r"^RP[LS]\d")].sum(1) / tot,
                         "frac_mt": pb.loc[:, g.str.startswith("MT-")].sum(1) / tot,
                         "frac_cyto_markers": pb[[c for c in cyto if c in g]].sum(1) / tot,
                         "total_counts": tot})


def deseq_disc(pb: pd.DataFrame, md: pd.DataFrame, cfg, covars: list[str]) -> pd.DataFrame:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    p = cfg["pseudobulk"]
    md = md.copy()
    for c in covars:
        md[c] = (md[c] - md[c].mean()) / md[c].std(ddof=1)
    counts = pb.loc[:, (pb >= p["min_counts_gene"]).sum(0) >= p["min_samples_gene"]].round().astype(int)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(counts=counts, metadata=md[["group", *covars]], design_factors=[*covars, "group"],
                           continuous_factors=covars or None, ref_level=["group", "Control"], n_cpus=16, quiet=True)
        dds.deseq2()
        st = DeseqStats(dds, contrast=["group", "DKD", "HKD"], quiet=True); st.summary()
    return st.results_df


def gsea_block(z: pd.Series, name: str, cfg, seed) -> pd.DataFrame:
    rows = []
    for lib, f in m09.GMTS.items():
        r = m09.prerank(z, REF / f, seed, cfg); r.insert(0, "library", lib); r.insert(0, "ranking", name)
        rows.append(r)
    return pd.concat(rows)


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); R = cfg["robustness"]
    out = stage_dir(STAGE)
    ex = exclusion_genes(cfg)
    axis = pd.read_csv(ROOT / "results/03b_technical_axis/technical_axis_by_lineage.tsv", sep="\t", index_col=0)["PT"]
    cyto = axis.sort_values().head(R["n_cyto_markers"]).index.tolist()     # most SC_A-enriched = cytoplasmic
    # discovery SN iPT pseudobulk
    q = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    pbf = ad.read_h5ad(ROOT / cfg["paths"]["processed"] / "pb_fine.h5ad")
    o = pbf.obs
    sel = ((o.celltype == "iPT") & (o.library.map(q.batch) == "SN") & o.library.map(q.qc_pass).astype(bool)
           & (o.n_cells >= cfg["pseudobulk"]["min_cells_per_sample"])).values
    X = pbf.X[sel]; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    dpb = pd.DataFrame(X, index=o.index[sel], columns=pbf.var_names)
    dmd = o.loc[sel, ["donor", "group"]].copy()
    # KPMP sn iPT pseudobulk
    kpb = pd.read_csv(ROOT / "results/07c_kpmp_concordance/iPT_pseudobulk_sn_all_calls.tsv.gz", sep="\t", index_col=0)
    groups = importlib.import_module("00b_download_kpmp").participant_groups().set_index("Participant ID")
    kmd = pd.DataFrame({"group": groups.group.reindex(kpb.index).values,
                        "sex": groups["Sex"].reindex(kpb.index).astype(str).values}, index=kpb.index)
    # 1. technical fractions and group tests
    tf_rows, tests = [], []
    for name, pb, md in [("discovery_SN", dpb, dmd), ("KPMP_sn", kpb, kmd)]:
        tf = tech_fractions(pb, cyto).join(md[["group"]]); tf["cohort"] = name
        tf_rows.append(tf)
        for c in ["frac_ribo", "frac_mt", "frac_cyto_markers"]:
            for a, b in [("DKD", "HKD"), ("DKD", "HKD_withDM"), ("HKD_withDM", "HKD")]:
                r = compare(tf[c], tf.group, a, b)
                if r.get("n_a", 0) >= 2 and r.get("n_b", 0) >= 2:
                    tests.append({"cohort": name, "metric": c, "contrast": f"{a}-{b}", **r})
    TF = pd.concat(tf_rows); TF.to_csv(out / "technical_fractions_by_sample.tsv", sep="\t")
    T = pd.DataFrame(tests); T.to_csv(out / "technical_fraction_tests.tsv", sep="\t", index=False)
    # 2. signature effects vs technical axis
    comb = pd.read_csv(ROOT / "results/09_signature_biology/combined_ranking.tsv", sep="\t", index_col=0)
    comb = comb.join(axis.rename("tech_axis"), how="inner")
    comb["technical_gene"] = is_technical(comb.index, ex)
    ax_rows = []
    for subset, d in [("all", comb), ("replicated", comb[comb.replicated]),
                      ("replicated_non_technical", comb[comb.replicated & ~comb.technical_gene])]:
        for col in ["disc_lfc", "kpmp_lfc", "z_stouffer"]:
            rr = stats.spearmanr(d[col], d.tech_axis, nan_policy="omit")
            ax_rows.append({"subset": subset, "effect": col, "n": len(d), "spearman_vs_tech_axis": rr.correlation, "p": rr.pvalue})
    A = pd.DataFrame(ax_rows); A.to_csv(out / "signature_vs_technical_axis.tsv", sep="\t", index=False)
    # 3. concordance: (i) exclusion, (ii) ribo-adjusted DESeq2 in both cohorts
    disc = pd.read_csv(ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    kp = pd.read_csv(ROOT / "results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv", sep="\t", index_col=0)
    tech_d = pd.Series(is_technical(disc.index, ex), index=disc.index)
    masks = {"disc_dE_fdr10": disc.dE_padj < 0.10,
             "disc_dE_fdr10_non_technical": (disc.dE_padj < 0.10) & ~tech_d,
             "all_tested_non_technical": disc.dE_padj.notna() & ~tech_d}
    conc = [dict(r, analysis="original") for r in m07c.concordance(disc, kp, "dE_lfc", masks)]
    dmd2 = dmd.join(TF[TF.cohort == "discovery_SN"][["frac_ribo"]])
    disc_adj = deseq_disc(dpb, dmd2, cfg, ["frac_ribo"])
    disc_adj = disc_adj.rename(columns={"log2FoldChange": "dE_lfc", "lfcSE": "dE_se", "padj": "dE_padj", "pvalue": "dE_p"})
    disc_adj.to_csv(out / "discovery_iPT_DKD_vs_HKD_adj_frac_ribo.tsv", sep="\t")
    kmd2 = kmd.join(TF[TF.cohort == "KPMP_sn"][["frac_ribo"]])
    kp_adj = m07c.deseq_contrast(kpb, kmd2, cfg, "DKD", "HKD", ["frac_ribo"])
    kp_adj.to_csv(out / "kpmp_sn_iPT_DKD_vs_HKD_adj_frac_ribo.tsv", sep="\t")
    tech_a = pd.Series(is_technical(disc_adj.index, ex), index=disc_adj.index)
    masks_a = {"disc_dE_fdr10": disc_adj.dE_padj < 0.10,
               "disc_dE_fdr10_non_technical": (disc_adj.dE_padj < 0.10) & ~tech_a,
               "orig_disc_dE_fdr10": disc.dE_padj.reindex(disc_adj.index) < 0.10}
    conc += [dict(r, analysis="ribo_adjusted_both") for r in m07c.concordance(disc_adj, kp_adj, "dE_lfc", masks_a)]
    C = pd.DataFrame(conc); C.to_csv(out / "concordance_robustness.tsv", sep="\t", index=False)
    # 4. GSEA re-runs
    z0 = comb.z_stouffer
    z_ex = z0[~comb.technical_gene]
    za = m09.wald(disc_adj, "dE_lfc", "dE_se"); zk = m09.wald(kp_adj, "log2FoldChange", "lfcSE")
    g = za.index.intersection(zk.index); z_adj = (za[g] + zk[g]) / np.sqrt(2)
    G = pd.concat([gsea_block(z_ex, "stouffer_excluding_technical_genes", cfg, seed),
                   gsea_block(z_adj, "stouffer_ribo_adjusted", cfg, seed),
                   gsea_block(z_adj[~pd.Series(is_technical(g, ex), index=g)], "stouffer_ribo_adjusted_excluding_technical", cfg, seed)])
    G.to_csv(out / "gsea_robustness.tsv", sep="\t", index=False)
    orig = pd.read_csv(ROOT / "results/09_signature_biology/gsea_prerank.tsv", sep="\t")
    orig = orig[orig.ranking == "DKD_vs_HKD_stouffer"].assign(ranking="original")
    key = pd.concat([orig, G])
    key = key[key.term.isin(R["key_terms"])].pivot_table(index="term", columns="ranking", values="NES")
    fdr = pd.concat([orig, G]); fdr = fdr[fdr.term.isin(R["key_terms"])].pivot_table(index="term", columns="ranking", values="fdr")
    key.join(fdr, rsuffix="_fdr").to_csv(out / "key_pathways_robustness.tsv", sep="\t")
    pd.set_option("display.width", 250)
    print(T.round(4).to_string()); print(A.round(3).to_string()); print(C.round(4).to_string())
    print(key.round(2).to_string()); print(fdr.round(3).to_string())
    write_provenance(STAGE, [ROOT / "data/processed/pb_fine.h5ad",
                             ROOT / "results/07c_kpmp_concordance/iPT_pseudobulk_sn_all_calls.tsv.gz",
                             ROOT / "results/09_signature_biology/combined_ranking.tsv",
                             ROOT / "results/03b_technical_axis/technical_axis_by_lineage.tsv"],
                     sorted(out.glob("*.tsv")), seed, extra={"n_excluded_terms_genes": len(ex)})


if __name__ == "__main__":
    main()
