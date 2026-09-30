"""03d_composition: within-lineage state fractions per donor (plan rule 3, composition arm).

For each state in config state_fractions the fraction state/lineage is computed per library and
tested with a binomial GLM (quasi-likelihood scale, cluster-robust SE by donor):
logit(frac) ~ batch + group. Reported per stratum: all QC libraries and SN only."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "03d_composition"


def main():
    cfg = load_config(); set_global_seed(cfg["seed"])
    comp = pd.read_csv(ROOT / "results/02_discovery_pseudobulk/composition_fine.tsv", sep="\t", index_col=0)
    q = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    q = q[q.qc_pass]
    comp = comp.loc[comp.index.intersection(q.index)]
    rows, fr = [], []
    for lin, state in cfg["state_fractions"].items():
        cts = [c for c in cfg["lineage_map"][lin] if c in comp.columns]
        d = q[["donor", "group", "batch"]].copy()
        d["k"] = comp.loc[d.index, state]; d["n"] = comp.loc[d.index, cts].sum(1)
        d = d[d.n >= cfg["pseudobulk"]["min_cells_per_sample"]]
        d["frac"] = d.k / d.n; d["lineage"] = lin; d["state"] = state
        fr.append(d.reset_index())
        for stratum, dd in [("all_libraries", d), ("SN_only", d[d.batch == "SN"])]:
            dd = dd.copy(); dd["group"] = pd.Categorical(dd.group, ["Control", "HKD", "DKD"])
            form = "frac ~ C(batch) + group" if dd.batch.nunique() > 1 else "frac ~ group"
            m = smf.glm(form, dd, family=sm.families.Binomial(), var_weights=dd.n).fit(
                scale="X2", cov_type="cluster", cov_kwds={"groups": pd.factorize(dd.donor)[0]})
            ref = dd.group.cat.categories
            for a, b in [("DKD", "Control"), ("HKD", "Control"), ("DKD", "HKD")]:
                L = np.zeros(len(m.params)); names = list(m.params.index)
                if a != "Control": L[names.index(f"group[T.{a}]")] += 1
                if b != "Control": L[names.index(f"group[T.{b}]")] -= 1
                t = m.t_test(L)
                rows.append({"lineage": lin, "state": state, "stratum": stratum, "contrast": f"{a}-{b}",
                             "log_odds": float(t.effect), "se": float(t.sd), "p": float(t.pvalue),
                             "n_libraries": len(dd), "n_donors": dd.donor.nunique()})
    out = stage_dir(STAGE)
    res = pd.DataFrame(rows); res.to_csv(out / "state_fraction_tests.tsv", sep="\t", index=False)
    pd.concat(fr).to_csv(out / "state_fractions_by_library.tsv", sep="\t", index=False)
    print(res.round(4).to_string())
    print(pd.concat(fr).groupby(["state", "group", "batch"]).frac.median().round(3))
    write_provenance(STAGE, [ROOT / "results/02_discovery_pseudobulk/composition_fine.tsv",
                             ROOT / "results/02c_library_qc/library_qc.tsv"], sorted(out.glob("*.tsv")), cfg["seed"])


if __name__ == "__main__":
    main()
