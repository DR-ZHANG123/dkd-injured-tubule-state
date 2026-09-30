"""03c_lodo_sn: leave-one-donor-out stability of the SN-stratum contrasts (plan rule 5).

For each lineage / fine type the SN-only DESeq2 fit is repeated leaving out each donor in turn
(designs that drop a group below 2 donors are skipped). Per gene: fraction of replicates whose
dE and dD log2FC keep the sign of the full fit, and fraction with dE FDR < 0.10."""
from __future__ import annotations
import importlib, sys
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
run_deseq = importlib.import_module("03_pseudobulk_de").run_deseq

STAGE = "03c_lodo_sn"


def fit(pb, sub, cfg):
    idx = [pb.obs_names.get_loc(i) for i in sub.index]
    X = pb.X[idx]; X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    return run_deseq(pd.DataFrame(X, index=sub.index, columns=pb.var_names), sub, cfg, n_cpus=2)


def main():
    cfg = load_config(); set_global_seed(cfg["seed"])
    proc = ROOT / cfg["paths"]["processed"]
    q = pd.read_csv(ROOT / "results/02c_library_qc/library_qc.tsv", sep="\t").set_index("library")
    out = stage_dir(STAGE)
    targets = cfg["lodo"]["targets"]
    jobs, subs, pbs = [], {}, {}
    for level, cts in targets.items():
        pb = ad.read_h5ad(proc / f"pb_{level}.h5ad")
        pb.obs["batch"] = pb.obs.library.map(q.batch).astype(str)
        pbs[level] = pb
        for ct in cts:
            o = pb.obs
            sub = o[(o.celltype == ct) & (o.batch == "SN") & o.library.map(q.qc_pass).astype(bool)
                    & (o.n_cells >= cfg["pseudobulk"]["min_cells_per_sample"])]
            subs[(level, ct)] = sub
            for d in sorted(sub.donor.unique()):
                s2 = sub[sub.donor != d]
                if s2.groupby("group").donor.nunique().reindex(cfg["cohort"]["groups"]).fillna(0).min() >= 2:
                    jobs.append((level, ct, d))
    res = Parallel(n_jobs=32)(delayed(fit)(pbs[l], subs[(l, c)][subs[(l, c)].donor != d], cfg) for l, c, d in jobs)
    summ = []
    for level, cts in targets.items():
        for ct in cts:
            full = pd.read_csv(ROOT / f"results/03_pseudobulk_de/{level}_SN_only/{ct.replace('/', '_')}.tsv",
                               sep="\t").set_index("gene")
            reps = [r for (l, c, d), r in zip(jobs, res) if l == level and c == ct]
            sE = pd.concat([np.sign(r.dE_lfc).reindex(full.index) == np.sign(full.dE_lfc) for r in reps], axis=1)
            sD = pd.concat([np.sign(r.dD_lfc).reindex(full.index) == np.sign(full.dD_lfc) for r in reps], axis=1)
            fq = pd.concat([(r.dE_padj < 0.10).reindex(full.index) for r in reps], axis=1)
            t = pd.DataFrame({"gene": full.index, "n_reps": len(reps), "dE_sign_stable": sE.mean(1).values,
                              "dD_sign_stable": sD.mean(1).values, "dE_frac_fdr10": fq.mean(1).values})
            t.to_csv(out / f"{level}_{ct.replace('/', '_')}.tsv", sep="\t", index=False)
            sig = full.dE_padj < 0.10
            summ.append({"level": level, "celltype": ct, "n_reps": len(reps), "n_dE_full": int(sig.sum()),
                         "n_dE_sign_stable80": int((t.set_index("gene").dE_sign_stable[sig] >= 0.8).sum()),
                         "median_frac_fdr10_among_sig": float(t.set_index("gene").dE_frac_fdr10[sig].median())
                         if sig.any() else np.nan})
    s = pd.DataFrame(summ); s.to_csv(out / "lodo_summary.tsv", sep="\t", index=False); print(s.to_string())
    write_provenance(STAGE, [proc / "pb_lineage.h5ad", proc / "pb_fine.h5ad"], sorted(out.glob("*.tsv")), cfg["seed"])


if __name__ == "__main__":
    main()
