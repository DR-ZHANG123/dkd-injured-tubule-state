"""08_validate_ercb_arrays: locked programmes in microdissected tubulointerstitium (GSE104954).

Platforms are analysed separately (GPL24120: DKD vs hypertensive nephropathy, the direct
cross-disease contrast; GPL22945: DKD vs living donor). Scores = mean per-gene z over samples
of the platform. Bulk tissue mixes compartments, so the DKD-biased score is additionally
modelled conditional on the shared-injury score."""
from __future__ import annotations
import gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from stats import compare, matched_random_null

STAGE = "08_validate_ercb_arrays"
DX = {"Diabetic nephropathy": "DKD", "Hypertensive nephropathy": "HTN", "Living donor": "LD"}


def read_matrix(path: Path) -> pd.DataFrame:
    lines = gzip.open(path, "rt").read().split("!series_matrix_table_begin\n")[1].split("!series_matrix_table_end")[0]
    from io import StringIO
    return pd.read_csv(StringIO(lines), sep="\t", index_col=0)


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed)
    raw = ROOT / cfg["paths"]["raw"]
    hg = pd.read_csv(raw / "reference/hgnc_complete_set.txt", sep="\t", usecols=["symbol", "entrez_id", "status"], low_memory=False)
    hg = hg[(hg.status == "Approved") & hg.entrez_id.notna()]
    e2s = dict(zip(hg.entrez_id.astype(int).astype(str) + "_at", hg.symbol))
    samples = pd.read_csv(ROOT / "results/01_cohort_audit/sample_table_GSE104954.tsv", sep="\t").set_index("sample")
    gs = pd.read_csv(ROOT / "results/05_lock_programme/gene_sets.tsv", sep="\t")
    sets = {k: v.gene.tolist() for k, v in gs.groupby("set_id")}
    out = stage_dir(STAGE)
    rows, adj = [], []
    inputs = []
    for plat, contrasts in [("GPL24120", [("DKD", "HTN"), ("DKD", "LD"), ("HTN", "LD")]), ("GPL22945", [("DKD", "LD")])]:
        f = raw / f"GSE104954/GSE104954-{plat}_series_matrix.txt.gz"; inputs.append(f)
        m = read_matrix(f)
        m.index = m.index.map(lambda x: e2s.get(x)); m = m[m.index.notna()]
        m = m.groupby(level=0).mean()
        grp = samples.diagnosis.reindex(m.columns).map(DX)
        keep = grp.notna(); m, grp = m.loc[:, keep], grp[keep]
        expr = m.T
        expr = expr.loc[:, expr.mean() > np.quantile(expr.mean(), 0.25)]   # drop the lowest-expressed quartile
        z = (expr - expr.mean()) / expr.std(ddof=1)
        sc_ = pd.DataFrame({k: z[[g for g in v if g in z.columns]].mean(1) for k, v in sets.items()})
        sc_.join(grp.rename("group")).to_csv(out / f"scores_{plat}.tsv", sep="\t")
        for a, b in contrasts:
            for k in sets:
                r = compare(sc_[k], grp, a, b)
                r.update(matched_random_null(z, sets[k], grp, a, b, expr.mean(), seed=seed))
                rows.append({"platform": plat, "set": k, "contrast": f"{a}-{b}",
                             "n_measured": len([g for g in sets[k] if g in z.columns]), **r})
            d = sc_.join(grp.rename("group"))
            d = d[d.group.isin([a, b])].copy(); d["case"] = (d.group == a).astype(int)
            for k in ["DKDiPT_up", "DEonly_up"]:
                mdl = smf.ols(f"{k} ~ case + SharedInjury_up", d).fit(cov_type="HC3")
                adj.append({"platform": plat, "set": k, "contrast": f"{a}-{b}", "n": len(d),
                            "beta_case_adj_injury": mdl.params["case"], "se": mdl.bse["case"], "p": mdl.pvalues["case"],
                            "beta_injury": mdl.params["SharedInjury_up"]})
        print(plat, grp.value_counts().to_dict())
    T = pd.DataFrame(rows); T.to_csv(out / "programme_tests.tsv", sep="\t", index=False)
    A = pd.DataFrame(adj); A.to_csv(out / "injury_adjusted.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(T[["platform", "set", "contrast", "n_measured", "n_a", "n_b", "hedges_g", "welch_p", "p_emp_greater"]].round(4).to_string())
    print(A.round(4).to_string())
    write_provenance(STAGE, inputs + [ROOT / "results/05_lock_programme/gene_sets.tsv"],
                     [out / "programme_tests.tsv", out / "injury_adjusted.tsv", *sorted(out.glob("scores_*.tsv"))], seed)


if __name__ == "__main__":
    main()
