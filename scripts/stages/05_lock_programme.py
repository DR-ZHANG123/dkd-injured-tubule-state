"""05_lock_programme: apply the marker-gene selection rules and fix the
gene sets used by every validation stage. Nothing downstream may alter these lists.

Gene sets written to results/05_lock_programme/gene_sets.tsv (set, gene, direction, rank):
  DKDiPT_up / DKDiPT_down   DL-supported DKD-biased iPT programme
  DEonly_up / DEonly_down   baseline: same size, ranked by pseudobulk dE p-value only
  SharedInjury_up / _down   dD and dN significant and concordant, dE not significant (rule 10)
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "05_lock_programme"


def main():
    cfg = load_config(); set_global_seed(cfg["seed"]); L = cfg["lock"]
    f_ckg = ROOT / "results/04_scdisinfact/primary/ckg_scores.tsv"
    f_de = ROOT / "results/03_pseudobulk_de/fine_SN_only/iPT.tsv"
    f_lodo = ROOT / "results/03c_lodo_sn/fine_iPT.tsv"
    f_lin = ROOT / "results/03_pseudobulk_de/lineage_SN_only/PT.tsv"
    ckg = pd.read_csv(f_ckg, sep="\t").set_index("gene")
    de = pd.read_csv(f_de, sep="\t").set_index("gene")
    lodo = pd.read_csv(f_lodo, sep="\t").set_index("gene")
    t = de.join(lodo[["dE_sign_stable", "dD_sign_stable", "dE_frac_fdr10"]]).join(ckg[["pct_rank", "mean_rank"]])
    t["ckg_top"] = t.pct_rank <= L["ckg_top_frac"]
    t["de_pass"] = (t.dE_padj < L["dE_fdr"]) & (np.sign(t.dE_lfc) == np.sign(t.dD_lfc))
    t["lodo_pass"] = t.dE_sign_stable >= L["lodo_sign_stable"]
    t["direction"] = np.where(t.dE_lfc > 0, "up", "down")
    hg = pd.read_csv(ROOT / cfg["paths"]["raw"] / "reference/hgnc_complete_set.txt", sep="\t",
                     usecols=["symbol", "location", "status"], low_memory=False)
    hg = hg[hg.status == "Approved"].set_index("symbol")
    sexg = set(hg.index[hg.location.astype(str).str.match(r"^[XY](p|q|$| )")])
    t["excluded"] = (t.index.isin(sexg) | t.index.str.match(L["exclude_regex"])
                     | (L["require_hgnc_approved"] & ~t.index.isin(hg.index)))
    t["programme"] = t.ckg_top & t.de_pass & t.lodo_pass & ~t.excluded
    rows, report = [], {}
    for d in ["up", "down"]:
        c = t[t.programme & (t.direction == d)].copy()
        c["abs_lfc"] = c.dE_lfc.abs()
        c = c.sort_values("abs_lfc", ascending=False).head(L["max_genes"])
        report[f"DKDiPT_{d}_n_eligible"] = int((t.programme & (t.direction == d)).sum())
        report[f"DKDiPT_{d}_n"] = len(c)
        report[f"DKDiPT_{d}_reported"] = len(c) >= L["min_genes"]
        for i, g in enumerate(c.index):
            rows.append(("DKDiPT", d, g, i + 1))
        # baseline: same size, pseudobulk evidence only (no model, no LODO), among the CKG input genes
        b = t[t.pct_rank.notna() & t.de_pass & ~t.excluded & (t.direction == d)].sort_values("dE_p").head(len(c))
        for i, g in enumerate(b.index):
            rows.append(("DEonly", d, g, i + 1))
        report[f"DEonly_{d}_overlap_with_DKDiPT"] = len(set(b.index) & set(c.index))
    # shared injury programme at PT-lineage level (iPT expansion is shared, 03d)
    lin = pd.read_csv(f_lin, sep="\t").set_index("gene")
    lin = lin[~(lin.index.isin(sexg) | lin.index.str.match(L["exclude_regex"])
                | (L["require_hgnc_approved"] & ~lin.index.isin(hg.index)))]
    si = lin[(lin.dD_padj < L["dE_fdr"]) & (lin.dN_padj < L["dE_fdr"]) &
             (np.sign(lin.dD_lfc) == np.sign(lin.dN_lfc)) & (lin.dE_padj > L["shared_dE_min_padj"])].copy()
    si["mean_lfc"] = (si.dD_lfc + si.dN_lfc) / 2
    for d, asc in [("up", False), ("down", True)]:
        s = si[(si.mean_lfc > 0) == (d == "up")].sort_values("mean_lfc", ascending=asc).head(L["max_genes"])
        report[f"SharedInjury_{d}_n"] = len(s)
        for i, g in enumerate(s.index):
            rows.append(("SharedInjury", d, g, i + 1))
    gs = pd.DataFrame(rows, columns=["set", "direction", "gene", "rank"])
    gs["set_id"] = gs.set + "_" + gs.direction
    out = stage_dir(STAGE)
    gs.to_csv(out / "gene_sets.tsv", sep="\t", index=False)
    t.to_csv(out / "gene_evidence_iPT.tsv", sep="\t")
    report["n_ckg_input_genes"] = int(ckg.shape[0])
    report["n_de_pass"] = int(t.de_pass.sum()); report["n_ckg_and_de"] = int((t.ckg_top & t.de_pass).sum())
    (out / "lock_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=1))
    for sid, g in gs.groupby("set_id"):
        print(sid, len(g), g.gene.tolist()[:50])
    write_provenance(STAGE, [f_ckg, f_de, f_lodo, f_lin, ROOT / cfg["paths"]["raw"] / "reference/hgnc_complete_set.txt"], [out / "gene_sets.tsv", out / "gene_evidence_iPT.tsv",
                     out / "lock_report.json"], cfg["seed"], extra={"rules": L})


if __name__ == "__main__":
    main()
