"""Fig4 — tissue, histology and clinical correspondence (plotting only).

Reads results/12b (Zenodo gene level), 10b (KPMP Visium), 11c (spot-score reliability),
11b (H&E leave-one-participant-out prediction), 12c (KPMP clinical) and the KPMP Visium
hires H&E image + spot table for the example crop. Rerunning regenerates every panel."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from _style import (setup, panel_label, source, save, RES, ROOT, GROUP_COL, SET_COL, SET_LABEL)

ZEN_COL = {"Control": "#8C8C8C", "DM_noDKD": "#CCB974", "other_CKD": "#937860", "DKD": "#C44E52"}
GROUP_ORDER = ["Reference", "DM_noCKD", "HKD", "HKD_withDM", "DKD"]
GROUP_LABEL = {"Reference": "Ref", "DM_noCKD": "DM", "HKD": "HKD", "HKD_withDM": "HKD+DM", "DKD": "DKD",
               "Control": "Ctrl", "DM_noDKD": "DM", "other_CKD": "CKD"}


def jitter(n, w=0.12, seed=0):
    return np.random.default_rng(int(seed * 10)).uniform(-w, w, n)


def panel_a(ax):
    d = pd.read_csv(RES / "12b_zenodo_gene_level/donor_celltype_expression.tsv", sep="\t")
    loc = pd.read_csv(RES / "12b_zenodo_gene_level/iPT_vs_PT_localisation.tsv", sep="\t")
    rows = []
    for (plat, gene), _ in loc.groupby(["platform", "gene"], sort=False):
        col = f"lcpm_{gene}"
        s = d[d.tech == plat].pivot_table(index=["donor", "group"], columns="celltype", values=col).dropna()
        s = s.reset_index()
        for r in s.itertuples():
            rows.append({"platform": plat, "gene": gene, "donor": r.donor, "group": r.group,
                         "iPT_minus_PT": r.iPT - r.PT})
    t = pd.DataFrame(rows)
    order = [(p, g) for p, g in loc[["platform", "gene"]].itertuples(index=False)]
    for i, (p, g) in enumerate(order):
        s = t[(t.platform == p) & (t.gene == g)]
        ax.scatter(i + jitter(len(s), seed=i), s.iPT_minus_PT, s=7, lw=0,
                   c=[ZEN_COL.get(x, "#8C8C8C") for x in s.group])
        ax.hlines(s.iPT_minus_PT.median(), i - 0.25, i + 0.25, color="k", lw=1)
        fr = loc[(loc.platform == p) & (loc.gene == g)].frac_donors_iPT_higher.iloc[0]
        ax.text(i, 1.02, f"{fr * 100:.0f}%", transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6)
    ax.text(-0.6, 1.02, "donors iPT > PT:", transform=ax.get_xaxis_transform(), ha="right", va="bottom", fontsize=6)
    ax.axhline(0, color="#999999", lw=0.6, ls="--")
    hs = [plt.Line2D([], [], ls="", marker="o", ms=3, mfc=ZEN_COL[g], mec="none", label=GROUP_LABEL[g])
          for g in ["Control", "DM_noDKD", "other_CKD", "DKD"]]
    ax.legend(handles=hs, loc="upper right", fontsize=6, ncol=2, handletextpad=0.1, columnspacing=0.6, borderpad=0.3)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_xticks(range(len(order)), labels=[f"{g}\n{p}" for p, g in order], fontsize=6.5)
    for lab in ax.get_xticklabels():
        lab.set_fontstyle("italic")
    ax.set_ylabel("iPT − PT log CPM\n(per donor)")
    source(t, "Fig4_panel_a")


def box_strip(ax, df, col, groups, x0, colour_by_group=True, width=0.6):
    for i, g in enumerate(groups):
        v = df.loc[df.group == g, col].dropna().values
        x = x0 + i
        ax.boxplot(v, positions=[x], widths=width, showfliers=False, patch_artist=True,
                   boxprops=dict(facecolor="white", lw=0.6), medianprops=dict(color="k", lw=0.9),
                   whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
        ax.scatter(x + jitter(len(v), 0.15, i), v, s=6, lw=0, c=GROUP_COL[g], zorder=3)


def panel_b(ax):
    s = pd.read_csv(RES / "10b_visium_kpmp/participant_scores.tsv", sep="\t")
    tst = pd.read_csv(RES / "10b_visium_kpmp/participant_tests.tsv", sep="\t")
    groups = [g for g in GROUP_ORDER if g in set(s.group)]
    xt, xl = [], []
    for k, sname in enumerate(["DKDiPT_up", "SharedInjury_up"]):
        x0 = k * (len(groups) + 1)
        box_strip(ax, s, sname, groups, x0)
        r = tst[(tst.metric == sname) & (tst.contrast == "DKD-HKD")].iloc[0]
        ax.text(x0 + (len(groups) - 1) / 2, 1.02, f"{SET_LABEL[sname]}\nDKD vs HKD g = {r.hedges_g:.2f}, P = {r.welch_p:.2f}",
                transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6)
        xt += [x0 + i for i in range(len(groups))]; xl += [GROUP_LABEL[g] for g in groups]
    ax.set_xticks(xt); ax.set_xticklabels(xl, rotation=90, fontsize=6)
    ax.set_ylabel("Visium score, PT-rich spots\n(participant mean)")
    source(s[["participant", "group", "n_PT_rich_spots", "DKDiPT_up", "SharedInjury_up"]], "Fig4_panel_b")


def panel_c(ax):
    r = pd.read_csv(RES / "11c_histology_reliability/spot_score_reliability.tsv", sep="\t")
    r = r[r.set.isin(["DKDiPT_up", "SharedInjury_up"])]
    pos, lab = 0, []
    for c, cname in [("gse211785", "GSE211785"), ("kpmp", "KPMP")]:
        for sname in ["DKDiPT_up", "SharedInjury_up"]:
            v = r[(r.cohort == c) & (r.set == sname)].reliability_sb.dropna().values
            ax.boxplot(v, positions=[pos], widths=0.55, showfliers=False, patch_artist=True,
                       boxprops=dict(facecolor="white", lw=0.6), medianprops=dict(color="k", lw=0.9),
                       whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
            ax.scatter(pos + jitter(len(v), 0.14, pos), v, s=5, lw=0, c=SET_COL[sname], zorder=3)
            lab.append(f"{'DKD state' if sname == 'DKDiPT_up' else 'Common inj.'}\n{cname}")
            pos += 1
        pos += 0.5
    ax.set_xticks([0, 1, 2.5, 3.5]); ax.set_xticklabels(lab, fontsize=6.5)
    ax.axhline(0, color="#999999", lw=0.6, ls="--")
    ax.set_ylabel("Split-half reliability\n(Spearman–Brown, per section)")
    source(r, "Fig4_panel_c")


def panel_d(ax):
    rows = []
    for c, cname in [("gse211785", "GSE211785"), ("kpmp", "KPMP")]:
        p = pd.read_csv(RES / f"11b_histology_association/prediction_lopo_{c}.tsv", sep="\t")
        p["cohort"] = cname; rows.append(p)
    p = pd.concat(rows)
    cats = [("SharedInjury_up", "encoder"), ("SharedInjury_up", "morphology"),
            ("DKDiPT_up", "encoder"), ("DKDiPT_up", "morphology")]
    pos, xt, xl = 0, [], []
    for cname in ["GSE211785", "KPMP"]:
        for tgt, feat in cats:
            v = p[(p.cohort == cname) & (p.target == tgt) & (p.features == feat)].pearson_r.dropna().values
            col = SET_COL[tgt]
            ax.scatter(pos + jitter(len(v), 0.14, pos), v, s=5, lw=0, c=col if feat == "encoder" else "#BBBBBB", zorder=3)
            ax.hlines(np.mean(v), pos - 0.3, pos + 0.3, color="k", lw=1)
            xt.append(pos); xl.append("Enc." if feat == "encoder" else "Morph.")
            pos += 1
        ax.text(pos - 2.5, 1.02, cname, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=7)
        pos += 0.6
    ax.axhline(0, color="#999999", lw=0.6, ls="--")
    ax.set_xticks(xt); ax.set_xticklabels(xl, fontsize=6, rotation=60, ha="right", rotation_mode="anchor")
    hs = [plt.Line2D([], [], ls="", marker="o", ms=3, mfc=SET_COL[t], mec="none", label=f"{SET_LABEL[t]} target")
          for t in ["SharedInjury_up", "DKDiPT_up"]]
    hs.append(plt.Line2D([], [], ls="", marker="o", ms=3, mfc="#BBBBBB", mec="none", label="morphology-only features"))
    ax.legend(handles=hs, loc="upper right", fontsize=5.5, handletextpad=0.1, borderpad=0.3)
    ax.set_ylim(-0.4, 1.05)
    ax.set_ylabel("Held-out Pearson r\n(one point per participant fold)")
    source(p, "Fig4_panel_d")


def panel_e(ax):
    c = pd.read_csv(RES / "12c_kpmp_clinical/clinical_associations.tsv", sep="\t")
    c = c[(c.modality == "sn") & (c.analysis == "group_adjusted_rank_ols")]
    var = [("albuminuria", "Albuminuria band"), ("proteinuria", "Proteinuria band"), ("egfr_band", "eGFR band"),
           ("if_pct", "Interstitial fibrosis %"), ("ta_pct", "Tubular atrophy %"),
           ("tubular_injury_pct", "Tubular injury %"), ("interstitial_wbc_pct", "Interstitial WBC %"),
           ("arteriolar_hyalinosis_grade", "Arteriolar hyalinosis")]
    y = np.arange(len(var))[::-1]
    out = []
    for k, (sname, dy) in enumerate([("DKDiPT_up", 0.17), ("SharedInjury_up", -0.17)]):
        for yi, (v, lab) in zip(y, var):
            r = c[(c.score == sname) & (c.variable == v)]
            if r.empty:
                continue
            r = r.iloc[0]
            filled = r.fdr_bh < 0.05
            ax.scatter(r.estimate, yi + dy, s=18, marker="o", c=SET_COL[sname] if filled else "white",
                       edgecolors=SET_COL[sname], lw=0.8, zorder=3, label=SET_LABEL[sname] if yi == y[0] else None)
            out.append({"score": sname, "variable": v, "n": r.n, "estimate": r.estimate, "p": r.p, "fdr_bh": r.fdr_bh})
    ax.axvline(0, color="#999999", lw=0.6, ls="--")
    ax.set_yticks(y); ax.set_yticklabels([l for _, l in var], fontsize=6.5)
    ax.set_xlabel("Rank-regression coefficient\n(adjusted for diagnosis)")
    ax.legend(loc="upper left", fontsize=6, handletextpad=0.3, borderpad=0.4, title="filled: FDR < 0.05",
              title_fontsize=6, ncol=2, columnspacing=0.8)
    ax.set_ylim(-0.6, len(var) + 1.2)
    ax.set_xlim(-1.0, 1.25)
    source(pd.DataFrame(out), "Fig4_panel_e")


def pick_example():
    sp = pd.read_csv(RES / "10b_visium_kpmp/spot_scores.tsv.gz", sep="\t")
    cand = sp[sp.group == "DKD"].groupby(["participant", "sample"]).PT_rich.sum().sort_values(ascending=False)
    for (pid, sample), _ in cand.items():
        d = ROOT / "data/processed/kpmp_visium" / pid / sample / "outs/spatial"
        if (d / "tissue_hires_image.png").exists():
            return sp[(sp.participant == pid) & (sp["sample"] == sample)], d, pid
    raise FileNotFoundError("no example section")


def panel_f(ax, cax):
    sp, d, pid = pick_example()
    sf = json.loads((d / "scalefactors_json.json").read_text())
    img = plt.imread(d / "tissue_hires_image.png")
    s = sf["tissue_hires_scalef"]
    x, y = sp.pxl_col_fullres.values * s, sp.pxl_row_fullres.values * s
    half = 260
    grid = [(gx, gy) for gx in np.linspace(x.min(), x.max(), 25) for gy in np.linspace(y.min(), y.max(), 25)]
    cx, cy = max(grid, key=lambda c: np.sum((np.abs(x - c[0]) < half * 0.8) & (np.abs(y - c[1]) < half * 0.8)))
    x0, x1 = int(max(0, cx - half)), int(min(img.shape[1], cx + half))
    y0, y1 = int(max(0, cy - half)), int(min(img.shape[0], cy + half))
    ax.imshow(img[y0:y1, x0:x1], interpolation="nearest")
    m = (x >= x0) & (x < x1) & (y >= y0) & (y < y1)
    z = sp.SharedInjury_up.values
    z = (z - np.nanmean(z)) / np.nanstd(z)
    dia = sf["spot_diameter_fullres"] * s
    sc = ax.scatter(x[m] - x0, y[m] - y0, c=np.clip(z[m], -2.5, 2.5), cmap="magma", s=(dia * 0.38) ** 2 / 4,
                    lw=0, alpha=0.55, vmin=-2.5, vmax=2.5)
    # 200 µm scale bar: spot diameter corresponds to 55 µm
    px_per_um = dia / 55.0
    L = 200 * px_per_um
    ax.plot([12, 12 + L], [(y1 - y0) - 14] * 2, color="k", lw=2)
    ax.text(12 + L / 2, (y1 - y0) - 22, "200 µm", ha="center", va="bottom", fontsize=6)
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    cb = plt.colorbar(sc, cax=cax)
    cb.set_label("Common injury score (z)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    ax.text(0.0, -0.04, f"KPMP DKD participant {pid}", transform=ax.transAxes, ha="left", va="top", fontsize=6)
    source(sp.loc[m, ["barcode", "pxl_row_fullres", "pxl_col_fullres", "SharedInjury_up", "DKDiPT_up", "PT_rich"]]
           .assign(participant=pid), "Fig4_panel_f")


def main():
    setup()
    fig = plt.figure(figsize=(7.2, 8.4))
    gs = GridSpec(3, 2, figure=fig, hspace=0.75, wspace=0.42, left=0.11, right=0.97, top=0.95, bottom=0.07)
    A = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]),
         fig.add_subplot(gs[1, 1]), fig.add_subplot(gs[2, 0])]
    sub = gs[2, 1].subgridspec(1, 2, width_ratios=[1, 0.05], wspace=0.05)
    axf, cax = fig.add_subplot(sub[0, 0]), fig.add_subplot(sub[0, 1])
    panel_a(A[0]); panel_b(A[1]); panel_c(A[2]); panel_d(A[3]); panel_e(A[4]); panel_f(axf, cax)
    for ax, l in zip(A + [axf], "abcdef"):
        panel_label(ax, l)
    for ax in A:
        ax.title.set_position((0, 1.0))
    # extra headroom for the annotation lines above a, b, d
    for ax in [A[0], A[1], A[3]]:
        ax.set_title(ax._left_title.get_text(), fontweight="bold", fontsize=11, loc="left", pad=22)
    for ax in [A[2], A[4], axf]:
        ax.set_title(ax._left_title.get_text(), fontweight="bold", fontsize=11, loc="left", pad=22)
    save(fig, "Fig4_tissue_histology_clinical", rows=[[A[0], A[1]], [A[2], A[3]], [A[4], axf]])


if __name__ == "__main__":
    main()
