"""Fig2 discovery (plotting only).

a  library QC: mitochondrial vs ribosomal read fraction, coloured by protocol batch
b  DKD-HKD effect vs nuclear/cytoplasmic technical axis (r over FDR<0.10 genes), SC-A vs SN strata
c  DKD-HKD genes (FDR<0.10) per lineage, SC-A vs SN strata
d  leave-one-donor-out stability in the SN stratum
e  iPT fraction of proximal tubule by group (SN stratum)
f  scDisInFact seed-to-seed rank agreement of condition-associated gene scores
g  locked-programme gene evidence: CKG rank vs pseudobulk dE significance in iPT
"""
from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from _style import RES, setup, panel_label, source, save, GROUP_COL, BATCH_COL, SET_COL

LINEAGES = ["PT", "TAL", "DCT_CNT", "Stroma", "Endo", "Myeloid", "Lymphoid"]


def main():
    setup()
    fig = plt.figure(figsize=(7.2, 8.4))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 1.25], hspace=0.75, wspace=0.62)

    # a
    q = pd.read_csv(RES / "02c_library_qc/library_qc.tsv", sep="\t")
    ax = fig.add_subplot(gs[0, 0])
    mk = {"Control": "o", "HKD": "s", "DKD": "^"}
    for (b, g), d in q[q.qc_pass].groupby(["batch", "group"]):
        ax.scatter(d.frac_mt * 100, d.frac_ribo * 100, s=14, marker=mk[g], color=BATCH_COL[b], edgecolor="white", lw=0.3)
    f = q[~q.qc_pass]
    ax.scatter(f.frac_mt * 100, f.frac_ribo * 100, s=16, marker="x", color="black", lw=0.8)
    ax.axvline(15, color="#999999", lw=0.6, ls="--")
    ax.set_xlabel("Mitochondrial reads (%)"); ax.set_ylabel("Ribosomal reads (%)")
    h = [plt.Line2D([], [], ls="", marker="o", color=BATCH_COL[b], ms=3.5, label=b.replace("_", "-")) for b in BATCH_COL]
    h += [plt.Line2D([], [], ls="", marker=mk[g], color="#555555", ms=3.5, label=g) for g in mk]
    h += [plt.Line2D([], [], ls="", marker="x", color="black", ms=3.5, label="QC fail")]
    ax.legend(handles=h, fontsize=5.5, loc="upper right", ncol=2, handletextpad=0.2, columnspacing=0.6,
              borderpad=0.3, labelspacing=0.2)
    ax.set_ylim(-1, 42)
    panel_label(ax, "a")
    source(q[["library", "donor", "group", "tech", "batch", "qc_pass", "frac_mt", "frac_ribo"]], "Fig2_panel_a")

    # b
    t = pd.read_csv(RES / "03b_technical_axis/dE_vs_technical_axis.tsv", sep="\t")
    t = t[t.lineage.isin(LINEAGES) & t["mode"].isin(["SC_A_only", "SN_only"])]
    ax = fig.add_subplot(gs[0, 1:])
    x = np.arange(len(LINEAGES))
    for j, (mode, lab) in enumerate([("SC_A_only", "SC-A stratum"), ("SN_only", "SN stratum")]):
        d = t[t["mode"] == mode].set_index("lineage").reindex(LINEAGES)
        ax.bar(x + (j - 0.5) * 0.36, d.r_sig.fillna(0), width=0.34, color=BATCH_COL[mode.replace("_only", "")],
               edgecolor="white", lw=0.3, label=lab)
        for xi, (v, n) in enumerate(zip(d.r_sig, d.n_sig)):
            if np.isnan(v):
                ax.text(xi + (j - 0.5) * 0.36, 0.02, "n.a.", ha="center", va="bottom", fontsize=5, rotation=90)
    ax.axhline(0, color="black", lw=0.5)
    ax.set_xticks(x); ax.set_xticklabels([l.replace("_", "/") for l in LINEAGES], fontsize=6.5)
    ax.set_ylabel("r (DKD-HKD vs\ntechnical axis)"); ax.set_ylim(-0.45, 1.05)
    ax.legend(fontsize=6, loc="upper right", ncol=2)
    panel_label(ax, "b")
    source(t, "Fig2_panel_b")

    # c
    s = pd.read_csv(RES / "03_pseudobulk_de/screen_summary.tsv", sep="\t")
    s = s[(s.level == "lineage") & s.celltype.isin(LINEAGES) & s["mode"].isin(["SC_A_only", "SN_only"])]
    ax = fig.add_subplot(gs[1, :2])
    for j, mode in enumerate(["SC_A_only", "SN_only"]):
        d = s[s["mode"] == mode].set_index("celltype").reindex(LINEAGES)
        ax.bar(x + (j - 0.5) * 0.36, d.dE_nDE_fdr10.fillna(0), width=0.34, color=BATCH_COL[mode.replace("_only", "")],
               edgecolor="white", lw=0.3, label=mode.replace("_only", "").replace("_", "-") + " stratum")
    ax.set_xticks(x); ax.set_xticklabels([l.replace("_", "/") for l in LINEAGES], fontsize=6.5)
    ax.set_ylabel("DKD vs HKD genes\n(FDR < 0.10)"); ax.set_ylim(0, s.dE_nDE_fdr10.max() * 1.2)
    ax.legend(fontsize=6, loc="upper right", ncol=2)
    panel_label(ax, "c")
    source(s[["celltype", "mode", "n_donors_Control", "n_donors_HKD", "n_donors_DKD", "dE_nDE_fdr10", "dD_nDE_fdr10", "dN_nDE_fdr10"]], "Fig2_panel_c")

    # d
    lo = pd.read_csv(RES / "03c_lodo_sn/lodo_summary.tsv", sep="\t")
    names = {"DCT_CNT": "DCT/CNT", "ThinLimb": "Thin limb", "C_TAL": "Cortical TAL", "Fibroblast_1": "Fibroblast 1"}
    lo["name"] = [("  " if l == "fine" else "") + names.get(c, c) for l, c in zip(lo.level, lo.celltype)]
    ax = fig.add_subplot(gs[1, 2])
    y = np.arange(len(lo))
    ax.barh(y, lo.median_frac_fdr10_among_sig, color="#8172B3", height=0.6, edgecolor="white", lw=0.3)
    ax.set_yticks(y); ax.set_yticklabels(lo.name, fontsize=6); ax.invert_yaxis()
    ax.set_xlim(0, 1.0); ax.set_xlabel("Median fraction of\nLODO fits with FDR < 0.10")
    panel_label(ax, "d")
    source(lo, "Fig2_panel_d")

    # e
    fr = pd.read_csv(RES / "03d_composition/state_fractions_by_library.tsv", sep="\t")
    fr = fr[(fr.state == "iPT") & (fr.batch == "SN")]
    tests = pd.read_csv(RES / "03d_composition/state_fraction_tests.tsv", sep="\t")
    tests = tests[(tests.state == "iPT") & (tests.stratum == "SN_only")]
    ax = fig.add_subplot(gs[2, 0])
    grps = ["Control", "HKD", "DKD"]
    rng = np.random.default_rng(0)
    for i, g in enumerate(grps):
        v = fr[fr.group == g].frac.values
        ax.boxplot(v, positions=[i], widths=0.5, showfliers=False, medianprops=dict(color="black", lw=0.8),
                   boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=12, color=GROUP_COL[g], edgecolor="white", lw=0.3, zorder=3)
    ax.set_xticks(range(3)); ax.set_xticklabels(grps); ax.set_ylabel("iPT fraction of PT")
    ax.set_ylim(0, 1.25)
    pv = tests.set_index("contrast").p
    for (a, b, yy) in [(0, 1, 1.02), (0, 2, 1.13), (1, 2, 1.02)]:
        c = f"{grps[b]}-{grps[a]}" if f"{grps[b]}-{grps[a]}" in pv else f"{grps[a]}-{grps[b]}"
        off = 0.04 if (a, b) == (1, 2) else 0
        ax.plot([a + 0.05 + off, b - 0.05], [yy, yy], color="black", lw=0.5)
        ax.text((a + b) / 2, yy + 0.01, f"P = {pv[c]:.2g}", ha="center", va="bottom", fontsize=5.5)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    panel_label(ax, "e")
    source(fr[["library", "donor", "group", "k", "n", "frac"]], "Fig2_panel_e")

    # f
    sr = pd.read_csv(RES / "04_scdisinfact/primary/seed_rank_spearman.tsv", sep="\t", index_col=0)
    ax = fig.add_subplot(gs[2, 1])
    im = ax.imshow(sr.values, vmin=0.8, vmax=1.0, cmap="Blues", aspect="auto")
    ax.set_xticks(range(5)); ax.set_yticks(range(5))
    ax.set_xticklabels([c.replace("seed", "") for c in sr.columns]); ax.set_yticklabels([c.replace("seed", "") for c in sr.index])
    ax.set_xlabel("Seed"); ax.set_ylabel("Seed")
    for i in range(5):
        for j in range(5):
            if i != j:
                ax.text(j, i, f"{sr.values[i, j]:.2f}", ha="center", va="center", fontsize=5,
                        color="white" if sr.values[i, j] > 0.9 else "black")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cb.ax.tick_params(labelsize=6)
    cb.set_label("Spearman ρ", fontsize=6.5)
    ax.spines[["left", "bottom"]].set_visible(False); ax.tick_params(length=0)
    panel_label(ax, "f")
    source(sr.reset_index().rename(columns={"index": "seed"}), "Fig2_panel_f")

    # g
    ev = pd.read_csv(RES / "05_lock_programme/gene_evidence_iPT.tsv", sep="\t")
    ev = ev[ev.pct_rank.notna() & ev.dE_p.notna()].copy()
    ev["neglog10p"] = -np.log10(ev.dE_p)
    ax = fig.add_subplot(gs[2, 2])
    bg = ev[~ev.programme]
    ax.scatter(bg.pct_rank * 100, bg.neglog10p, s=2, color="#CCCCCC", lw=0, rasterized=True)
    for d, col in [("up", SET_COL["DKDiPT_up"]), ("down", SET_COL["DKDiPT_down"])]:
        pg = ev[ev.programme & (ev.direction == d)]
        ax.scatter(pg.pct_rank * 100, pg.neglog10p, s=10, color=col, edgecolor="black", lw=0.3, zorder=3,
                   label=f"Marker genes {d} ({len(pg)})")
    ax.axvline(10, color="#999999", lw=0.6, ls="--")
    ax.set_xlabel("CKG rank percentile"); ax.set_ylabel("−log10 P (DKD vs HKD)")
    ax.set_xlim(0, 100)
    ax.legend(fontsize=5.5, loc="upper right", handletextpad=0.2, borderpad=0.3)
    panel_label(ax, "g")
    source(ev[["gene", "pct_rank", "dE_lfc", "dE_p", "dE_padj", "dE_sign_stable", "excluded", "programme", "direction"]], "Fig2_panel_g")

    axes = fig.axes
    save(fig, "Fig2_discovery", rows=[[axes[0], axes[1]], [axes[2], axes[3]], [axes[4], axes[5], axes[7]]])


if __name__ == "__main__":
    main()
