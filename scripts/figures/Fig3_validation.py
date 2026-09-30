"""Fig3 validation (plotting only).

a  GSE195460: DKD-biased iPT score by donor (iPT and all-PT pseudobulk)
b  GSE104954: DKD-biased and shared-injury scores in tubulointerstitium, both platforms
c  KPMP: Hedges g (95% CI) of programme scores, per modality and contrast
d  KPMP: transcriptome-wide sign concordance with the discovery DKD-vs-HKD iPT effects
e  discovery vs KPMP snRNA log2FC for discovery FDR<0.10 genes
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from _style import RES, setup, panel_label, source, save, GROUP_COL, SET_COL, SET_LABEL, hedges_ci

rng = np.random.default_rng(0)


def strip(ax, groups, data, col, ylab):
    for i, g in enumerate(groups):
        v = data[data.group == g][col].dropna().values
        ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color="black", lw=0.8),
                   boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=9, color=GROUP_COL.get(g, "#888888"),
                   edgecolor="white", lw=0.3, zorder=3)
    ax.set_ylabel(ylab)


def main():
    setup()
    fig = plt.figure(figsize=(7.2, 8.1))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1.35, 1.15], hspace=0.8, wspace=1.1)

    # a
    ax = fig.add_subplot(gs[0, 0:2])
    t6 = pd.read_csv(RES / "06_validate_gse195460/programme_tests.tsv", sep="\t").set_index(["compartment", "set"])
    rows = []
    for k, comp in enumerate(["iPT", "PT_all"]):
        d = pd.read_csv(RES / f"06_validate_gse195460/scores_{comp}.tsv", sep="\t").rename(columns={"diagnosis": "group"})
        for i, g in enumerate(["Control", "DKD"]):
            v = d[d.group == g].DKDiPT_up.values
            x = k * 2.6 + i
            ax.boxplot(v, positions=[x], widths=0.55, showfliers=False, medianprops=dict(color="black", lw=0.8),
                       boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
            ax.scatter(x + rng.uniform(-0.12, 0.12, len(v)), v, s=10, color=GROUP_COL[g], edgecolor="white", lw=0.3, zorder=3)
        r = t6.loc[(comp, "DKDiPT_up")]
        ax.text(k * 2.6 + 0.5, 1.02, f"g = {r.hedges_g:.2f}\nP = {r.welch_p:.2g}", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=5.5)
        rows.append(d.assign(compartment=comp)[["donor", "group", "compartment", "DKDiPT_up"]])
    ax.set_xticks([0, 1, 2.6, 3.6]); ax.set_xticklabels(["Ctrl", "DKD", "Ctrl", "DKD"])
    ax.set_xlim(-0.6, 4.2)
    for k, lab in enumerate(["iPT", "All PT"]):
        ax.text(k * 2.6 + 0.5, -0.2, lab, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.5)
    ax.set_ylabel("DKD marker-gene score")
    panel_label(ax, "a")
    source(pd.concat(rows), "Fig3_panel_a")

    # b
    ax = fig.add_subplot(gs[0, 2:])
    t8 = pd.read_csv(RES / "08_validate_ercb_arrays/programme_tests.tsv", sep="\t")
    rows, x0, ticks, labs = [], 0, [], []
    for plat, grps in [("GPL24120", ["HTN", "DKD"]), ("GPL22945", ["LD", "DKD"])]:
        d = pd.read_csv(RES / f"08_validate_ercb_arrays/scores_{plat}.tsv", sep="\t")
        for i, g in enumerate(grps):
            v = d[d.group == g].DKDiPT_up.values; x = x0 + i
            ax.boxplot(v, positions=[x], widths=0.55, showfliers=False, medianprops=dict(color="black", lw=0.8),
                       boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
            ax.scatter(x + rng.uniform(-0.12, 0.12, len(v)), v, s=8, color=GROUP_COL[g], edgecolor="white", lw=0.3, zorder=3)
            ticks.append(x); labs.append(g)
        r = t8[(t8.platform == plat) & (t8.set == "DKDiPT_up") & (t8.contrast == f"DKD-{grps[0]}")].iloc[0]
        ax.text(x0 + 0.5, 1.02, f"g = {r.hedges_g:.2f}\nP = {r.welch_p:.2g}", transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=5.5)
        ax.text(x0 + 0.5, -0.2, plat, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.5)
        rows.append(d.assign(platform=plat)[["group", "platform", "DKDiPT_up", "SharedInjury_up"]]); x0 += 2.6
    ax.set_xticks(ticks); ax.set_xticklabels(labs); ax.set_xlim(-0.6, 4.2)
    ax.set_ylabel("DKD marker-gene score\n(5 of 12 genes on array)")
    panel_label(ax, "b")
    source(pd.concat(rows), "Fig3_panel_b")

    # c: forest
    ax = fig.add_subplot(gs[1, 0:2])
    T = pd.read_csv(RES / "07b_validate_kpmp/programme_tests.tsv", sep="\t")
    contrasts = [("DKD-HKD", "DKD vs HKD"), ("DKD-HKD_withDM", "DKD vs HKD+DM"), ("HKD_withDM-HKD", "HKD+DM vs HKD"),
                 ("DKD-DM_noCKD", "DKD vs DM, no CKD")]
    sets = ["DKDiPT_up", "DEonly_up", "SharedInjury_up"]
    y, yt, yl, rows = 0, [], [], []
    for c, cl in contrasts:
        for mod, mk in [("sn", "o"), ("sc", "s")]:
            for j, sname in enumerate(sets):
                r = T[(T.modality == mod) & (T.contrast == c) & (T.set == sname)]
                if r.empty or pd.isna(r.hedges_g.iloc[0]):
                    continue
                r = r.iloc[0]; lo, hi = hedges_ci(r.hedges_g, r.n_a, r.n_b)
                yy = y + (j - 1) * 0.22
                ax.plot([lo, hi], [yy, yy], color=SET_COL[sname], lw=0.8)
                ax.plot(r.hedges_g, yy, marker=mk, color=SET_COL[sname], ms=3.2, mec="white", mew=0.3)
                rows.append({"modality": mod, "contrast": c, "set": sname, "n_a": r.n_a, "n_b": r.n_b,
                             "hedges_g": r.hedges_g, "ci_lo": lo, "ci_hi": hi, "welch_p": r.welch_p})
            yt.append(y); yl.append(f"{cl} ({mod}RNA)"); y += 1
        y += 0.4
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=5.8); ax.invert_yaxis()
    ax.set_xlabel("Hedges g (95% CI)")
    hs = [plt.Line2D([], [], color=SET_COL[s], marker="o", ms=3, lw=0.8, label=SET_LABEL[s]) for s in sets]
    ax.legend(handles=hs, fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, handlelength=1.2,
              columnspacing=0.8)
    panel_label(ax, "c")
    source(pd.DataFrame(rows), "Fig3_panel_c")

    # d: concordance
    ax = fig.add_subplot(gs[1, 2:])
    C = pd.read_csv(RES / "07c_kpmp_concordance/transcriptome_concordance.tsv", sep="\t")
    items = [("DKD_vs_HKD", "disc_dE_fdr10", "DKD vs HKD"), ("DKD_vs_HKD_adj_egfr_lower", "disc_dE_fdr10", "+ eGFR category"),
             ("DKD_vs_HKD_adj_if_pct", "disc_dE_fdr10", "+ fibrosis %"), ("DKD_vs_HKD_withDM", "disc_dE_fdr10", "DKD vs HKD+DM"),
             ("HKD_withDM_vs_HKD", "disc_dE_fdr10", "HKD+DM vs HKD"), ("DKD_vs_HKD", "shared_injury_genes", "Common injury")]
    P = pd.read_csv(RES / "14a_concordance_permutation/permutation_concordance.tsv", sep="\t")
    P = P[P.primary.astype(str) == "True"]
    perm_key = {("DKD_vs_HKD", "disc_dE_fdr10"): ("DKD_vs_HKD", "disc_dE_fdr10"),
                ("DKD_vs_HKD_adj_egfr_lower", "disc_dE_fdr10"): ("DKD_vs_HKD_adj_egfr_lower", "disc_dE_fdr10"),
                ("DKD_vs_HKD", "shared_injury_genes"): ("DKD_vs_HKD", "shared_injury_genes")}
    rows = []
    for i, (c, gset, lab) in enumerate(items):
        for j, (mod, col) in enumerate([("sn", "#4C72B0"), ("sc", "#DD8452")]):
            r = C[(C.modality == mod) & (C.contrast == c) & (C.gene_set == gset)].iloc[0]
            y = i + (j - 0.5) * 0.36
            ax.barh(y, r.frac_agree * 100, height=0.34, color=col, edgecolor="white", lw=0.3,
                    label=f"KPMP {mod}RNA" if i == 0 else None)
            row = {"modality": mod, "contrast": c, "gene_set": gset,
                   **r[["n_genes", "n_sign_agree", "frac_agree", "spearman_lfc"]].to_dict()}
            k = perm_key.get((c, gset))
            if k is not None:
                pr = P[(P.modality == mod) & (P.contrast == k[0]) & (P.gene_set == k[1])].iloc[0]
                ax.plot([pr.perm_null_q95 * 100] * 2, [y - 0.17, y + 0.17], color="black", lw=0.9)
                star = "***" if pr.perm_p < 1e-3 else "**" if pr.perm_p < 1e-2 else "*" if pr.perm_p < 0.05 else "n.s."
                ax.text(max(r.frac_agree, pr.perm_null_q95) * 100 + 1.5, y, star, va="center", fontsize=5.5)
                row.update({"perm_null_q95": pr.perm_null_q95, "perm_p": pr.perm_p})
            rows.append(row)
    ax.set_yticks(range(len(items))); ax.set_yticklabels([l for *_, l in items], fontsize=6); ax.invert_yaxis()
    ax.axhline(len(items) - 1.5, color="#BBBBBB", lw=0.5)
    ax.set_xlim(0, 100); ax.set_xlabel("Same-sign genes (%)")
    ax.legend(fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2)
    panel_label(ax, "d")
    source(pd.DataFrame(rows), "Fig3_panel_d")

    # e: lfc scatter
    ax = fig.add_subplot(gs[2, 0:2])
    disc = pd.read_csv(RES / "03_pseudobulk_de/fine_SN_only/iPT.tsv", sep="\t").set_index("gene")
    kp = pd.read_csv(RES / "07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv", sep="\t", index_col=0)
    j = disc[disc.dE_padj < 0.10][["dE_lfc"]].join(kp[["log2FoldChange"]], how="inner").dropna()
    agree = np.sign(j.dE_lfc) == np.sign(j.log2FoldChange)
    ax.scatter(j.dE_lfc[agree], j.log2FoldChange[agree], s=6, color="#4C72B0", lw=0, alpha=0.8, label="Same sign")
    ax.scatter(j.dE_lfc[~agree], j.log2FoldChange[~agree], s=6, color="#BBBBBB", lw=0, label="Opposite sign")
    ax.axhline(0, color="black", lw=0.4); ax.axvline(0, color="black", lw=0.4)
    rho = stats.spearmanr(j.dE_lfc, j.log2FoldChange).correlation
    ax.text(0.03, 0.97, f"n = {len(j)}; {agree.sum()} same sign\nSpearman ρ = {rho:.2f}", transform=ax.transAxes,
            ha="left", va="top", fontsize=6)
    ax.set_xlabel("Discovery log2FC (DKD vs HKD)"); ax.set_ylabel("KPMP snRNA log2FC")
    ax.legend(fontsize=5.5, loc="lower right")
    panel_label(ax, "e")
    source(j.reset_index().rename(columns={"index": "gene", "dE_lfc": "discovery_log2FC", "log2FoldChange": "kpmp_sn_log2FC"}), "Fig3_panel_e")

    # f: per-gene replication of the locked programme (snRNA)
    ax = fig.add_subplot(gs[2, 2:])
    G = pd.read_csv(RES / "07c_kpmp_concordance/per_gene_replication.tsv", sep="\t")
    G = G[(G.modality == "sn") & G.kpmp_lfc.notna()].copy()
    G = G.sort_values(["set", "disc_dE_lfc"], ascending=[False, False])
    yy = np.arange(len(G))
    cols = np.where(G.set == "DKDiPT_up", SET_COL["DKDiPT_up"], SET_COL["DKDiPT_down"])
    ax.barh(yy, G.kpmp_lfc, color=cols, height=0.65, edgecolor="white", lw=0.3)
    ax.set_yticks(yy); ax.set_yticklabels(G.gene, fontsize=5.3, fontstyle="italic"); ax.invert_yaxis()
    ax.axvline(0, color="black", lw=0.5)
    ax.set_xlabel("KPMP snRNA log2FC\n(DKD vs HKD)")
    hs = [plt.Rectangle((0, 0), 1, 1, color=SET_COL["DKDiPT_up"], label="Marker genes, up"), plt.Rectangle((0, 0), 1, 1, color=SET_COL["DKDiPT_down"], label="Marker genes, down")]
    ax.legend(handles=hs, fontsize=5.5, loc="upper left")
    panel_label(ax, "f")
    source(G, "Fig3_panel_f")

    A = fig.axes
    save(fig, "Fig3_validation", rows=[[A[0], A[1]], [A[2], A[3]], [A[4], A[5]]])


if __name__ == "__main__":
    main()
