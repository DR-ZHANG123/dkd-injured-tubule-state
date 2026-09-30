"""Fig1 study design and cohorts (plotting only).

a  study framework (vector icons + editable labels; drawn by Fig1a_framework.draw_framework)
b  discovery donor x modality matrix (GSE211785 + Zenodo overlap)
c  units per cohort and diagnostic group
d  discovery library flow: QC -> protocol batch -> SN stratum by group
"""
from __future__ import annotations
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from _style import RES, setup, panel_label, source, save, GROUP_COL
from Fig1a_framework import draw_framework, check_layout

UNIFIED = ["Control/reference", "HKD", "HKD with DM", "DM without CKD", "DKD", "Other CKD"]
UCOL = {"Control/reference": "#8C8C8C", "HKD": "#4C72B0", "HKD with DM": "#8172B3",
        "DM without CKD": "#CCB974", "DKD": "#C44E52", "Other CKD": "#64B5CD"}


def cohort_counts() -> pd.DataFrame:
    a = json.loads((RES / "01_cohort_audit/audit_summary.json").read_text())
    rows = []
    add = lambda cohort, unit, grp, n: rows.append({"cohort": cohort, "unit": unit, "group": grp, "n": int(n)})
    m = {"Control": "Control/reference", "HKD": "HKD", "DKD": "DKD"}
    for k, v in a["GSE211785_rna_donors_by_dx"].items():
        add("GSE211785 sc/snRNA", "donors", m[k], v)
    vm = {"Control": "Control/reference", "HKD": "HKD", "DKD": "DKD", "CKD": "Other CKD"}
    for k, v in a["GSE211785_visium_sections_by_dx"].items():
        add("GSE211785 Visium", "sections", vm[k], v)
    for k, v in a["GSE195460_rna_libraries_by_dx"].items():
        add("GSE195460 snRNA", "donors", m[k], v)
    em = {"Diabetic nephropathy": "DKD", "Hypertensive nephropathy": "HKD", "Living donor": "Control/reference"}
    for k, v in a["GSE104954_by_platform_dx"].items():
        plat, dx = k.split("|")
        if dx in em:
            add(f"GSE104954 {plat}", "biopsies", em[dx], v)
    km = {"DKD": "DKD", "HKD": "HKD", "HKD_withDM": "HKD with DM", "DM_noCKD": "DM without CKD", "Reference": "Control/reference"}
    for tag, lab in [("sn", "KPMP snRNA"), ("sc", "KPMP scRNA")]:
        p = pd.read_csv(RES / f"07b_validate_kpmp/participants_{tag}.tsv", sep="\t")
        for k, v in p.group.value_counts().items():
            add(lab, "participants", km[k], v)
    zm = {"DKD": "DKD", "Control": "Control/reference", "DM": "DM without CKD", "DM/HTN": "DM without CKD"}
    for k, v in a["zenodo_donors_by_dx"].items():
        add("Zenodo CosMx/Xenium", "donors", zm.get(k, "Other CKD"), v)
    return pd.DataFrame(rows).groupby(["cohort", "unit", "group"], as_index=False).n.sum()


def main():
    setup()
    fig = plt.figure(figsize=(7.2, 8.0))
    gs = fig.add_gridspec(3, 2, height_ratios=[2.75, 2.6, 1.9], width_ratios=[1.1, 1], hspace=0.55, wspace=0.75)
    # a: study framework
    axa = fig.add_subplot(gs[0, :])
    pos = axa.get_position()                   # full figure width, keep the framework aspect (no distortion)
    wf = 0.97; hf = wf * fig.get_figwidth() * 31 / 100 / fig.get_figheight()
    axa.set_position([0.015, pos.y1 - hf, wf, hf])
    fw = draw_framework(axa, scale=wf * fig.get_figwidth() / 7.2)
    panel_label(axa, "a")

    # b: donor x modality matrix
    mat = pd.read_csv(RES / "01_cohort_audit/donor_modality_matrix.tsv", sep="\t")
    mat["dx"] = mat.dx_GSE211785_sc.fillna(mat.dx_GSE211785_visium).fillna(mat.dx_zenodo)
    mat = mat[mat.dx.isin(["Control", "HKD", "DKD", "CKD"])].copy()
    order = {"Control": 0, "HKD": 1, "DKD": 2, "CKD": 3}
    mat["o"] = mat.dx.map(order)
    mat = mat.sort_values(["o", "donor"])
    mods = [("scRNA", "scRNA"), ("snRNA", "snRNA"), ("visium", "Visium"), ("visium_HE_image", "H&E image"),
            ("cosmx", "CosMx (Zenodo)"), ("xenium", "Xenium (Zenodo)")]
    M = np.array([[int(mat[c].iloc[i] > 0) for i in range(len(mat))] for c, _ in mods])
    axb = fig.add_subplot(gs[1, 0])
    code = M * (mat.o.values[None, :] + 1)
    cmap = ListedColormap(["#F2F2F2", GROUP_COL["Control"], GROUP_COL["HKD"], GROUP_COL["DKD"], "#64B5CD"])
    axb.imshow(code, aspect="auto", cmap=cmap, vmin=0, vmax=4, interpolation="nearest")
    axb.set_yticks(range(len(mods))); axb.set_yticklabels([l for _, l in mods], fontsize=6.5)
    axb.set_xticks([]); axb.set_xlabel(f"GSE211785 and Zenodo donors (n = {len(mat)})")
    bounds = np.cumsum(mat.groupby("o").size().values)[:-1]
    for bnd in bounds:
        axb.axvline(bnd - 0.5, color="white", lw=1.2)
    starts = np.r_[0, bounds]; ends = np.r_[bounds, len(mat)]
    names = {0: "Control", 1: "HKD", 2: "DKD", 3: "Other CKD"}
    for (o, g_), s, e in zip(mat.groupby("o"), starts, ends):
        axb.text((s + e - 1) / 2, -0.75, names[o], ha="center", va="bottom", fontsize=6.5)
    ov = mat[(mat.scRNA + mat.snRNA + mat.visium > 0) & (mat.cosmx + mat.xenium > 0)]
    for i, d in enumerate(mat.donor):
        if d in set(ov.donor):
            axb.plot(i, len(mods) - 0.35, marker="^", ms=3, color="black", clip_on=False)
    axb.set_ylim(len(mods) - 0.5 + 0.35, -0.5)
    axb.spines[["left", "bottom"]].set_visible(False)
    axb.tick_params(length=0)
    panel_label(axb, "b")
    src_b = mat[["donor", "dx"] + [c for c, _ in mods]].assign(shared_with_zenodo=mat.donor.isin(ov.donor))
    source(src_b, "Fig1_panel_b")

    # c: cohort sizes
    cc = cohort_counts()
    source(cc, "Fig1_panel_c")
    cohorts = ["GSE211785 sc/snRNA", "GSE211785 Visium", "GSE195460 snRNA", "GSE104954 GPL24120",
               "GSE104954 GPL22945", "KPMP snRNA", "KPMP scRNA", "Zenodo CosMx/Xenium"]
    axc = fig.add_subplot(gs[1, 1])
    left = np.zeros(len(cohorts))
    piv = cc.pivot_table(index="cohort", columns="group", values="n", aggfunc="sum").reindex(cohorts).fillna(0)
    for g_ in UNIFIED:
        v = piv.get(g_, pd.Series(0, index=cohorts)).values
        axc.barh(range(len(cohorts)), v, left=left, color=UCOL[g_], edgecolor="white", lw=0.3, label=g_, height=0.7)
        left += v
    units = cc.drop_duplicates("cohort").set_index("cohort").unit.reindex(cohorts)
    for i, (t, u) in enumerate(zip(left, units)):
        axc.text(t + 3, i, f"{int(t)} {u}", va="center", fontsize=6)
    axc.set_yticks(range(len(cohorts))); axc.set_yticklabels(cohorts, fontsize=6.5); axc.invert_yaxis()
    axc.set_xlim(0, left.max() * 1.55); axc.set_xlabel("Units analysed")
    axc.legend(fontsize=6, loc="upper center", bbox_to_anchor=(0.35, -0.2), ncol=3, handlelength=1,
               borderpad=0.3, labelspacing=0.25, columnspacing=0.8)
    panel_label(axc, "c")

    # d: discovery library flow
    q = pd.read_csv(RES / "02c_library_qc/library_qc.tsv", sep="\t")
    steps = [("All RNA libraries", q), ("QC pass", q[q.qc_pass])]
    for b in ["SC_A", "SC_B", "SN"]:
        steps.append((f"Batch {b.replace('_', '-')}", q[q.qc_pass & (q.batch == b)]))
    axd = fig.add_subplot(gs[2, :])
    rows = []
    for i, (lab, d) in enumerate(steps):
        left = 0
        for g_ in ["Control", "HKD", "DKD"]:
            n = int((d.group == g_).sum())
            axd.barh(i, n, left=left, color=GROUP_COL[g_], edgecolor="white", lw=0.3, height=0.65,
                     label=g_ if i == 0 else None)
            if n >= 3:
                axd.text(left + n / 2, i, str(n), ha="center", va="center", fontsize=6, color="white")
            left += n; rows.append({"step": lab, "group": g_, "n_libraries": n})
        axd.text(left + 0.6, i, f"{len(d)} libraries / {d.donor.nunique()} donors", va="center", fontsize=6)
    axd.set_yticks(range(len(steps))); axd.set_yticklabels([s for s, _ in steps], fontsize=6.5); axd.invert_yaxis()
    axd.set_xlim(0, len(q) * 1.45); axd.set_xlabel("Discovery RNA libraries")
    axd.axhspan(len(steps) - 1.5, len(steps) - 0.5, color="#4C72B0", alpha=0.07, lw=0)
    axd.legend(fontsize=6, loc="center right", ncol=1, handlelength=1, borderpad=0.3)
    panel_label(axd, "d")
    source(pd.DataFrame(rows), "Fig1_panel_d")
    # panel letter a shares the left edge of letters b and d (left column)
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    xb = fig.transFigure.inverted().transform((axb._left_title.get_window_extent(r).x0, 0))[0]
    pa = axa.get_position()
    axa._left_title.set_x((xb - pa.x0) / pa.width)
    fw_issues = check_layout(fig, axa, *fw)
    print("Fig1a layout issues:", len(fw_issues)); [print("  ", i) for i in fw_issues]
    save(fig, "Fig1_cohorts", rows=[[axb, axc]])


if __name__ == "__main__":
    main()
