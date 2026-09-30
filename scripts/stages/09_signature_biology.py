"""09_signature_biology: biological annotation of the DKD-vs-HKD iPT shift.

Replicated signature = genes with discovery SN iPT dE FDR < 0.10 whose KPMP snRNA DKD-vs-HKD
log2FC has the same sign. Ranking = Stouffer combination of the discovery and KPMP snRNA Wald
statistics (both DKD vs HKD, iPT pseudobulk). Preranked GSEA (gseapy) on MSigDB Hallmark,
Reactome and GO BP; the same ranking for the HKD_withDM-vs-HKD contrast shows which pathways
follow diabetic exposure."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "09_signature_biology"
GMTS = {"hallmark": "h.all.v2024.1.Hs.symbols.gmt", "reactome": "c2.cp.reactome.v2024.1.Hs.symbols.gmt",
        "gobp": "c5.go.bp.v2024.1.Hs.symbols.gmt"}


def wald(df: pd.DataFrame, lfc: str, se: str) -> pd.Series:
    return (df[lfc] / df[se]).replace([np.inf, -np.inf], np.nan)


def prerank(rnk: pd.Series, gmt: Path, seed: int, cfg) -> pd.DataFrame:
    import gseapy as gp
    g = cfg["gsea"]
    r = gp.prerank(rnk=rnk.dropna().sort_values(ascending=False), gene_sets=str(gmt), min_size=g["min_size"],
                   max_size=g["max_size"], permutation_num=g["permutations"], seed=seed, threads=16,
                   outdir=None, verbose=False)
    res = r.res2d.rename(columns={"Term": "term", "NES": "NES", "FDR q-val": "fdr", "NOM p-val": "p", "Lead_genes": "lead_genes"})
    return res[["term", "ES", "NES", "p", "fdr", "lead_genes"]].astype({"NES": float, "fdr": float, "p": float})


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed)
    out = stage_dir(STAGE)
    disc = pd.read_csv(ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    kp = pd.read_csv(ROOT / "results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv", sep="\t", index_col=0)
    kdm = pd.read_csv(ROOT / "results/07c_kpmp_concordance/kpmp_iPT_HKD_withDM_vs_HKD_sn.tsv", sep="\t", index_col=0)
    zd, zk = wald(disc, "dE_lfc", "dE_se"), wald(kp, "log2FoldChange", "lfcSE")
    g = zd.index.intersection(zk.index)
    comb = pd.DataFrame({"z_discovery": zd[g], "z_kpmp_sn": zk[g]})
    comb["z_stouffer"] = (comb.z_discovery + comb.z_kpmp_sn) / np.sqrt(2)
    comb["z_kpmp_HKDwithDM_vs_HKD"] = wald(kdm, "log2FoldChange", "lfcSE").reindex(g)
    comb["disc_fdr"] = disc.dE_padj.reindex(g); comb["disc_lfc"] = disc.dE_lfc.reindex(g)
    comb["kpmp_lfc"] = kp.log2FoldChange.reindex(g)
    comb["replicated"] = (comb.disc_fdr < 0.10) & (np.sign(comb.disc_lfc) == np.sign(comb.kpmp_lfc))
    comb.sort_values("z_stouffer", ascending=False).to_csv(out / "combined_ranking.tsv", sep="\t")
    rows = []
    for rank_name, col in [("DKD_vs_HKD_stouffer", "z_stouffer"), ("HKDwithDM_vs_HKD_kpmp", "z_kpmp_HKDwithDM_vs_HKD")]:
        for lib, f in GMTS.items():
            res = prerank(comb[col], ROOT / "data/raw/reference" / f, seed, cfg)
            res.insert(0, "library", lib); res.insert(0, "ranking", rank_name)
            rows.append(res)
    R = pd.concat(rows); R.to_csv(out / "gsea_prerank.tsv", sep="\t", index=False)
    rep = comb[comb.replicated].sort_values("z_stouffer", ascending=False)
    rep.to_csv(out / "replicated_signature.tsv", sep="\t")
    pd.set_option("display.width", 250, "display.max_colwidth", 60)
    print("replicated genes:", len(rep), "up", int((rep.disc_lfc > 0).sum()), "down", int((rep.disc_lfc < 0).sum()))
    print(rep.head(30)[["disc_lfc", "kpmp_lfc", "z_stouffer"]].round(2))
    print(rep.tail(20)[["disc_lfc", "kpmp_lfc", "z_stouffer"]].round(2))
    for (rk, lib), d in R.groupby(["ranking", "library"]):
        d = d.sort_values("NES")
        print(f"== {rk} {lib}"); print(pd.concat([d.head(8), d.tail(8)])[["term", "NES", "fdr"]].round(3).to_string())
    write_provenance(STAGE, [ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv",
                             ROOT / "results/07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv",
                             ROOT / "results/07c_kpmp_concordance/kpmp_iPT_HKD_withDM_vs_HKD_sn.tsv",
                             *[ROOT / "data/raw/reference" / f for f in GMTS.values()]],
                     sorted(out.glob("*.tsv")), seed)


if __name__ == "__main__":
    main()
