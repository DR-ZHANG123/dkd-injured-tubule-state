"""Supplementary figures S1–S4 (plotting only).

S1 technical content of iPT pseudobulk and concordance after ribosomal adjustment (12a)
S2 pathway enrichment before / after technical adjustment (09, 12a)
S3 leave-one-donor-out stability in the SN stratum (03c)
S4 reference-based label-transfer accuracy (06, 07b)"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from _style import setup, panel_label, source, save, RES, GROUP_COL

GL = {"Control": "Ctrl", "Reference": "Ref", "HKD": "HKD", "HKD_withDM": "HKD+DM", "DKD": "DKD", "DM_noCKD": "DM"}


def jit(n, w=0.13, seed=0):
    return np.random.default_rng(int(seed * 10)).uniform(-w, w, n)


def figS1():
    t = pd.read_csv(RES / "12a_technical_robustness/technical_fractions_by_sample.tsv", sep="\t", index_col=0)
    tt = pd.read_csv(RES / "12a_technical_robustness/technical_fraction_tests.tsv", sep="\t")
    c = pd.read_csv(RES / "12a_technical_robustness/concordance_robustness.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 3.3))
    gs = GridSpec(1, 3, figure=fig, wspace=0.55, left=0.08, right=0.98, top=0.84, bottom=0.25)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for ax, metric, lab in [(axes[0], "frac_ribo", "Ribosomal read fraction"), (axes[1], "frac_cyto_markers", "Cytoplasmic-marker fraction")]:
        pos, xt, xl = 0, [], []
        for coh, groups in [("discovery_SN", ["Control", "HKD", "DKD"]), ("KPMP_sn", ["HKD", "HKD_withDM", "DKD"])]:
            for i, g in enumerate(groups):
                v = t[(t.cohort == coh) & (t.group == g)][metric].dropna().values * 100
                ax.scatter(pos + jit(len(v), seed=pos), v, s=6, lw=0, c=GROUP_COL[g], zorder=3)
                if len(v):
                    ax.hlines(np.median(v), pos - 0.3, pos + 0.3, color="k", lw=1)
                xt.append(pos); xl.append(GL[g]); pos += 1
            ax.text(pos - 2, 1.03, "Discovery SN" if coh.startswith("disc") else "KPMP snRNA",
                    transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.5)
            pos += 0.7
        ax.set_yscale("log")
        ax.set_xticks(xt); ax.set_xticklabels(xl, rotation=90, fontsize=6.5)
        ax.set_ylabel(f"{lab} (%)")
        r = tt[(tt.cohort == "KPMP_sn") & (tt.metric == metric) & (tt.contrast == "HKD_withDM-HKD")].iloc[0]
        ax.text(0.98, 0.03, f"KPMP HKD+DM vs HKD\ng = {r.hedges_g:.2f}, P = {r.welch_p:.3f}", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=6)
    ax = axes[2]
    rows = [("original", "disc_dE_fdr10", "Original"), ("original", "disc_dE_fdr10_non_technical", "Technical genes\nexcluded"),
            ("ribo_adjusted_both", "disc_dE_fdr10", "Ribosome-adjusted\n(both cohorts)"),
            ("ribo_adjusted_both", "orig_disc_dE_fdr10", "Original discovery vs\nadjusted KPMP")]
    vals = []
    for i, (an, gsn, lab) in enumerate(rows):
        r = c[(c.analysis == an) & (c.gene_set == gsn)].iloc[0]
        ax.barh(len(rows) - 1 - i, r.frac_agree * 100, color="#4C72B0", height=0.6, edgecolor="white", lw=0.3)
        ax.text(r.frac_agree * 100 + 1.5, len(rows) - 1 - i, f"{r.n_sign_agree}/{r.n_genes}", va="center", fontsize=6)
        vals.append({"analysis": an, "gene_set": gsn, "label": lab.replace("\n", " "), "n_genes": r.n_genes,
                     "n_sign_agree": r.n_sign_agree, "frac_agree": r.frac_agree, "binom_p": r.binom_p})
    ax.axvline(50, color="k", ls="--", lw=0.6)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[2] for r in rows][::-1], fontsize=6)
    ax.set_xlim(0, 100); ax.set_xlabel("Same-sign genes (%)")
    for a, l in zip(axes, "abc"):
        panel_label(a, l); a.set_title(l, fontweight="bold", fontsize=11, loc="left", pad=16)
    source(t.reset_index(), "FigS1_panel_ab"); source(pd.DataFrame(vals), "FigS1_panel_c")
    save(fig, "FigS1_technical_robustness", rows=[axes])
    plt.close(fig)


ACR = {"myc": "MYC", "il6": "IL6", "jak": "JAK", "stat3": "STAT3", "tnfa": "TNFα", "nfkb": "NF-κB", "kras": "KRAS",
       "bmp": "BMP", "dn": "down", "v1": "v1", "v2": "v2"}


def term_label(t: str) -> str:
    pre = {"HALLMARK": "Hallmark", "REACTOME": "Reactome", "GOBP": "GO BP"}[t.split("_")[0]]
    words = [ACR.get(w, w) for w in t.split("_", 1)[1].lower().split("_")]
    body = " ".join(words)
    return f"{pre}: {body[0].upper() + body[1:]}"


def figS2():
    k = pd.read_csv(RES / "12a_technical_robustness/key_pathways_robustness.tsv", sep="\t")
    cols = [("original", "Original"), ("stouffer_excluding_technical_genes", "Technical\ngenes\nexcluded"),
            ("stouffer_ribo_adjusted", "Ribosome-\nadjusted"), ("stouffer_ribo_adjusted_excluding_technical", "Adjusted\n+ excluded")]
    k = k.sort_values("original")
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    y = np.arange(len(k))
    for j, (cn, lab) in enumerate(cols):
        v = k[cn].values
        f = k[f"{cn}_fdr"].values if f"{cn}_fdr" in k else np.full(len(k), np.nan)
        for yi, vi, fi in zip(y, v, f):
            if np.isnan(vi):
                continue
            ax.scatter(j, yi, s=(abs(vi) * 22) ** 1.2, c="#C44E52" if vi > 0 else "#4C72B0",
                       alpha=1.0 if (not np.isnan(fi) and fi < 0.05) else 0.35, lw=0)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([c[1] for c in cols], fontsize=6.5)
    ax.set_yticks(y); ax.set_yticklabels([term_label(t) for t in k.term], fontsize=6)
    ax.set_xlim(-0.6, len(cols) - 0.4); ax.set_ylim(-0.8, len(k) - 0.2)
    hs = [plt.scatter([], [], s=(n * 22) ** 1.2, c="#8C8C8C", lw=0, label=f"|NES| = {n}") for n in (1, 2)]
    hs += [plt.scatter([], [], s=40, c="#C44E52", lw=0, label="up in DKD"), plt.scatter([], [], s=40, c="#4C72B0", lw=0, label="down in DKD")]
    ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=6, labelspacing=1.1,
              title="opaque: FDR < 0.05", title_fontsize=6)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    panel_label(ax, "")
    source(k, "FigS2")
    fig.subplots_adjust(left=0.36, right=0.78, top=0.95, bottom=0.1)
    save(fig, "FigS2_pathways_robustness")
    plt.close(fig)


def figS3():
    s = pd.read_csv(RES / "03c_lodo_sn/lodo_summary.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 2.9))
    gs = GridSpec(1, 2, figure=fig, wspace=0.45, left=0.1, right=0.98, top=0.86, bottom=0.26)
    a1, a2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    rows = []
    for i, r in enumerate(s.itertuples()):
        full = pd.read_csv(RES / f"03_pseudobulk_de/{r.level}_SN_only/{r.celltype.replace('/', '_')}.tsv", sep="\t").set_index("gene")
        lo = pd.read_csv(RES / f"03c_lodo_sn/{r.level}_{r.celltype.replace('/', '_')}.tsv", sep="\t").set_index("gene")
        v = lo.loc[full.index[full.dE_padj < 0.10], "dE_frac_fdr10"].dropna().values
        a2.boxplot(v, positions=[i], widths=0.6, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="white", lw=0.6), medianprops=dict(color="k", lw=0.9),
                   whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        rows.append({"level": r.level, "celltype": r.celltype, "n_dE_full": r.n_dE_full,
                     "n_sign_stable80": r.n_dE_sign_stable80, "median_frac_fdr10": np.median(v) if len(v) else np.nan})
    d = pd.DataFrame(rows)
    x = np.arange(len(d))
    a1.bar(x, d.n_dE_full, color="#4C72B0", width=0.6, edgecolor="white", lw=0.3, label="DKD vs HKD, FDR < 0.10")
    a1.scatter(x, d.n_sign_stable80, marker="_", s=120, c="k", lw=1.2, label="same sign in ≥ 80% of replicates", zorder=3)
    lab = [f"{c}" + ("*" if l == "lineage" else "") for l, c in zip(d.level, d.celltype)]
    for a in (a1, a2):
        a.set_xticks(x); a.set_xticklabels(lab, rotation=90, fontsize=6)
    a2.text(1.0, -0.42, "* lineage-level pseudobulk", transform=a2.transAxes, ha="right", va="top", fontsize=6)
    a1.set_ylabel("Genes (SN stratum)"); a1.set_ylim(0, d.n_dE_full.max() * 1.35)
    a1.legend(loc="upper right", fontsize=6)
    a2.set_ylabel("Fraction of leave-one-donor-out\nreplicates with FDR < 0.10")
    a2.set_ylim(0, 1.05)
    panel_label(a1, "a"); panel_label(a2, "b")
    source(d, "FigS3")
    save(fig, "FigS3_lodo_stability", rows=[[a1, a2]])
    plt.close(fig)


def figS4():
    g6 = pd.read_csv(RES / "06_validate_gse195460/reference_cv_accuracy.tsv", sep="\t", index_col=0).iloc[:, 0].rename("GSE195460 (SN ref.)")
    k = pd.read_csv(RES / "07b_validate_kpmp/reference_cv_accuracy.tsv", sep="\t", index_col=0)
    k.columns = ["KPMP snRNA (SN ref.)", "KPMP scRNA (SC ref.)"]
    m = pd.concat([g6, k], axis=1)
    m = m.loc[m.index.sort_values()]
    fig, ax = plt.subplots(figsize=(3.6, 6.2))
    im = ax.imshow(m.values, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_xticks(range(m.shape[1])); ax.set_xticklabels(m.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(m.shape[0])); ax.set_yticklabels(m.index, fontsize=5.5)
    for i, ct in enumerate(m.index):
        if ct in ("iPT", "PT_S1", "PT_S2", "PT_S3"):
            ax.get_yticklabels()[i].set_fontweight("bold")
        for j in range(m.shape[1]):
            v = m.iloc[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.5, color="white" if v < 0.6 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.03); cb.set_label("Donor-grouped CV accuracy", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    for s_ in ax.spines.values():
        s_.set_visible(False)
    source(m.reset_index().rename(columns={"index": "celltype"}), "FigS4")
    fig.subplots_adjust(left=0.3, right=0.92, top=0.98, bottom=0.2)
    save(fig, "FigS4_label_transfer_accuracy")
    plt.close(fig)


if __name__ == "__main__":
    setup()
    figS1(); figS2(); figS3(); figS4()
