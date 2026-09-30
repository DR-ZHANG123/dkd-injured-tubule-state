"""14g_exposure_trend: ordered trend of the marker-gene score across HKD < HKD with diabetes < DKD.

Participant-level scores from 07b; Jonckheere–Terpstra-type test implemented as Kendall tau
between score and the ordinal group code, with a label-permutation P (seeded)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "14g_exposure_trend"
ORDER = {"HKD": 0, "HKD_withDM": 1, "DKD": 2}


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); n_perm = cfg["exposure_trend"]["n_perm"]
    rows, ins = [], []
    for mod in ["sn", "sc"]:
        f = ROOT / f"results/07b_validate_kpmp/scores_{mod}.tsv"; ins.append(f)
        s = pd.read_csv(f, sep="\t", index_col=0)
        s = s[s.group.isin(ORDER)].dropna(subset=["DKDiPT_up"])
        x = s.group.map(ORDER).values
        for col in ["DKDiPT_up", "SharedInjury_up"]:
            y = s[col].values
            tau = stats.kendalltau(x, y).statistic
            rng = np.random.default_rng(seed)
            null = np.array([stats.kendalltau(rng.permutation(x), y).statistic for _ in range(n_perm)])
            rows.append({"modality": mod, "score": col, "n": len(s),
                         **{f"n_{k}": int((s.group == k).sum()) for k in ORDER},
                         **{f"mean_{k}": float(s.loc[s.group == k, col].mean()) for k in ORDER},
                         "kendall_tau": tau, "perm_p_one_sided": float((np.sum(null >= tau) + 1) / (n_perm + 1))})
    out = stage_dir(STAGE); T = pd.DataFrame(rows); T.to_csv(out / "exposure_trend.tsv", sep="\t", index=False)
    print(T.round(4).to_string())
    write_provenance(STAGE, ins, [out / "exposure_trend.tsv"], seed)


if __name__ == "__main__":
    main()
