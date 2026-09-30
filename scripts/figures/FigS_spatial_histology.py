"""Supplementary figures S5–S7 (plotting only).

S5 KPMP iPT fraction and programme-score distributions by group and modality (07b)
S6 Visium PT-rich spot selection: positive controls and programme gene coverage (10a, 10b)
S7 H&E–Visium registration QC and tissue filtering (11a)"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from _style import setup, panel_label, source, save, RES, ROOT, GROUP_COL, SET_LABEL

ORDER = ["Reference", "DM_noCKD", "HKD", "HKD_withDM", "DKD"]
GL = {"Reference": "Ref", "DM_noCKD": "DM", "HKD": "HKD", "HKD_withDM": "HKD+DM", "DKD": "DKD"}


def jit(n, w=0.14, seed=0):
    return np.random.default_rng(int(seed * 10)).uniform(-w, w, n)


def strip(ax, df, col, x0, seed=0):
    xs = []
    for i, g in enumerate(ORDER):
        v = df.loc[df.group == g, col].dropna().values
        x = x0 + i
        ax.boxplot(v, positions=[x], widths=0.6, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="white", lw=0.6), medianprops=dict(color="k", lw=0.9),
                   whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        ax.scatter(x + jit(len(v), seed=seed + i), v, s=5, lw=0, c=GROUP_COL[g], zorder=3)
        xs.append(x)
    return xs


def figS5():
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = GridSpec(2, 3, figure=fig, hspace=0.75, wspace=0.45, left=0.09, right=0.98, top=0.93, bottom=0.12)
    axes, out = [], []
    for r, tag in enumerate(["sn", "sc"]):
        p = pd.read_csv(RES / f"07b_validate_kpmp/participants_{tag}.tsv", sep="\t", index_col=0)
        s = pd.read_csv(RES / f"07b_validate_kpmp/scores_{tag}.tsv", sep="\t", index_col=0)
        for c, (df, col, lab) in enumerate([(p, "iPT_frac", "iPT fraction of PT"), (s, "DKDiPT_up", SET_LABEL["DKDiPT_up"]),
                                            (s, "SharedInjury_up", SET_LABEL["SharedInjury_up"])]):
            ax = fig.add_subplot(gs[r, c]); axes.append(ax)
            xs = strip(ax, df, col, 0, seed=r * 10 + c)
            ax.set_xticks(xs); ax.set_xticklabels([GL[g] for g in ORDER], rotation=90, fontsize=6.5)
            ax.set_ylabel(f"{lab}\n(KPMP {'snRNA' if tag == 'sn' else 'scRNA'})" if c == 0 else f"{lab} score", fontsize=7)
            out.append(df[["group", col]].assign(modality=tag, variable=col).rename(columns={col: "value"}).reset_index())
    for a, l in zip(axes, "abcdef"):
        panel_label(a, l)
    source(pd.concat(out), "FigS5")
    save(fig, "FigS5_kpmp_composition_scores", rows=[axes[:3], axes[3:]])
    plt.close(fig)


def figS6():
    pc = pd.read_csv(RES / "10a_visium_gse211785/deconvolution_positive_control.tsv", sep="\t")
    c1 = pd.read_csv(RES / "10a_visium_gse211785/set_coverage.tsv", sep="\t")
    c2 = pd.read_csv(RES / "10b_visium_kpmp/set_coverage.tsv", sep="\t")
    fig = plt.figure(figsize=(7.2, 2.9))
    gs = GridSpec(1, 2, figure=fig, wspace=0.35, left=0.33, right=0.98, top=0.86, bottom=0.18, width_ratios=[1.1, 1])
    a1, a2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    y = np.arange(len(pc))[::-1]
    a1.barh(y, pc.auc, color="#4C72B0", height=0.6, edgecolor="white", lw=0.3)
    for yi, v in zip(y, pc.auc):
        a1.text(v + 0.01, yi, f"{v:.2f}", va="center", fontsize=6)
    a1.axvline(0.5, color="k", ls="--", lw=0.6)
    import textwrap
    a1.set_yticks(y); a1.set_yticklabels([textwrap.fill(t.replace("_", " "), 40) for t in pc.check], fontsize=5.5)
    a1.set_xlim(0, 1.1); a1.set_xlabel("AUC (or fraction of spots)\nvs author spot labels, GSE211785")
    sets = ["DKDiPT_up", "DKDiPT_down", "DEonly_up", "SharedInjury_up", "SharedInjury_down"]
    m = pd.DataFrame({"GSE211785 (FFPE)": c1.set_index("set").n_measured.reindex(sets),
                      "KPMP (fresh frozen)": c2.set_index("set").n_measured.reindex(sets)})
    nset = c1.set_index("set").n_set.reindex(sets)
    x = np.arange(len(sets))
    for k, (col, colr) in enumerate(zip(m.columns, ["#8C8C8C", "#4C72B0"])):
        a2.bar(x + (k - 0.5) * 0.36, m[col] / nset * 100, width=0.34, color=colr, edgecolor="white", lw=0.3, label=col)
        for xi, v, n in zip(x, m[col], nset):
            a2.text(xi + (k - 0.5) * 0.36, v / n * 100 + 2, f"{int(v)}", ha="center", fontsize=5.5)
    a2.set_xticks(x); a2.set_xticklabels([SET_LABEL[s_].replace(" (DL + DE)", "") for s_ in sets], rotation=90, fontsize=6)
    a2.set_ylabel("Set genes measured (%)"); a2.set_ylim(0, 125)
    a2.legend(loc="upper right", fontsize=6)
    panel_label(a1, "a"); panel_label(a2, "b")
    source(pc, "FigS6_panel_a"); source(m.assign(n_set=nset).reset_index().rename(columns={"index": "set"}), "FigS6_panel_b")
    save(fig, "FigS6_visium_spot_selection", rows=[[a1, a2]])
    plt.close(fig)


def figS7():
    idx = pd.read_csv(RES / "11a_histology_embed/section_index.tsv", sep="\t")
    auc = idx.orient_tissue.fillna("{}").map(lambda s: sorted(json.loads(s).values()))
    idx["best_auc"] = auc.map(lambda v: v[-1] if v else np.nan)
    idx["second_auc"] = auc.map(lambda v: v[-2] if len(v) > 1 else np.nan)
    idx["tissue_pass_frac"] = idx.n_tissue_pass / idx.n_patches
    fig = plt.figure(figsize=(7.2, 3.0))
    gs = GridSpec(1, 3, figure=fig, wspace=0.5, left=0.08, right=0.98, top=0.86, bottom=0.2, width_ratios=[1, 1, 1.1])
    a1, a2, a3 = (fig.add_subplot(gs[0, i]) for i in range(3))
    for k, (coh, col) in enumerate([("gse211785", "#8C8C8C"), ("kpmp", "#4C72B0")]):
        d = idx[idx.cohort == coh]
        if d.best_auc.notna().any():
            a1.scatter(d.second_auc, d.best_auc, s=8, lw=0, c=col, label=f"GSE211785 (n = {d.best_auc.notna().sum()})")
        a2.scatter(k + jit(len(d), seed=k), d.tissue_pass_frac * 100, s=8, lw=0, c=col)
    a1.plot([0.4, 1], [0.4, 1], color="#999999", lw=0.6, ls="--")
    a1.set_xlim(0.4, 1.02); a1.set_ylim(0.4, 1.02)
    a1.set_xlabel("Runner-up orientation AUC"); a1.set_ylabel("Chosen orientation AUC\n(in- vs off-tissue spots)")
    a1.legend(loc="lower right", fontsize=6)
    a1.text(0.03, 0.45, "KPMP: Space Ranger\ncoordinates used directly", transform=a1.transAxes, ha="left", va="top", fontsize=5.5)
    a2.set_xticks([0, 1]); a2.set_xticklabels(["GSE211785", "KPMP"], fontsize=6.5)
    a2.set_xlim(-0.6, 1.6); a2.set_ylabel("Spots passing tissue QC (%)")
    # example: KPMP low-res image with tissue-pass flags
    f = pd.read_csv(RES / "11a_histology_embed/features_kpmp.tsv.gz", sep="\t")
    sec = f.groupby("section").tissue_pass.mean().sort_values().index[len(f.section.unique()) // 2]
    pid, sample = sec.split("__")
    d = ROOT / "data/processed/kpmp_visium" / pid / sample / "outs/spatial"
    sf = json.loads((d / "scalefactors_json.json").read_text())
    img = plt.imread(d / "tissue_lowres_image.png")
    ff = f[f.section == sec]
    s = sf["tissue_lowres_scalef"]
    a3.imshow(img)
    a3.set_anchor("N")   # top-align the image so its panel letter shares the row baseline
    for ok, colr, lab in [(True, "#55A868", "tissue QC pass"), (False, "#C44E52", "excluded")]:
        g = ff[ff.tissue_pass == ok]
        a3.scatter(g.px_col * s, g.px_row * s, s=0.8, lw=0, c=colr, label=lab)
    a3.set_xticks([]); a3.set_yticks([])
    for sp_ in a3.spines.values():
        sp_.set_visible(False)
    a3.legend(loc="upper left", bbox_to_anchor=(0.0, -0.02), fontsize=6, markerscale=5, ncol=2, frameon=False)
    a3.text(1.0, 1.01, f"KPMP {pid}", transform=a3.transAxes, ha="right", va="bottom", fontsize=6)
    for a, l in zip((a1, a2, a3), "abc"):
        panel_label(a, l)
    source(idx[["cohort", "section", "participant", "um_per_px", "n_spots", "n_tissue_pass", "orientation", "best_auc",
                "second_auc", "tissue_pass_frac"]], "FigS7_panel_ab")
    source(ff[["spot_id", "px_row", "px_col", "tissue_pass"]].assign(section=sec), "FigS7_panel_c")
    save(fig, "FigS7_histology_registration_qc", rows=[[a1, a2, a3]])
    plt.close(fig)


if __name__ == "__main__":
    setup()
    figS5(); figS6(); figS7()
