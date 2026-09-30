"""Shared figure style, save and layout checks (plotting only; never recomputes analyses).

Figure style: sans-serif 8 pt (Arial; Liberation Sans is the
metric-compatible substitute on this host), top/right spines off, bold lowercase panel letters,
PDF with editable text (fonttype 42) + 600-dpi PNG."""
from __future__ import annotations
from pathlib import Path
import itertools
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"
FIG = ROOT / "figures"
SRC = FIG / "source_data"
FIG.mkdir(exist_ok=True); SRC.mkdir(exist_ok=True)

PALETTE = ['#4C72B0', '#DD8452', '#55A868', '#C44E52', '#8172B3', '#937860', '#DA8BC3', '#8C8C8C', '#CCB974', '#64B5CD']
GROUP_COL = {"Control": "#8C8C8C", "Reference": "#8C8C8C", "HKD": "#4C72B0", "HTN": "#4C72B0",
             "HKD_withDM": "#8172B3", "DM_noCKD": "#CCB974", "DKD": "#C44E52", "LD": "#8C8C8C"}
BATCH_COL = {"SC_A": "#DD8452", "SC_B": "#55A868", "SN": "#4C72B0"}
SET_COL = {"DKDiPT_up": "#C44E52", "DEonly_up": "#DD8452", "SharedInjury_up": "#55A868",
           "DKDiPT_down": "#937860", "SharedInjury_down": "#64B5CD"}
SET_LABEL = {"DKDiPT_up": "DKD-associated state", "DEonly_up": "DE-only gene set",
             "SharedInjury_up": "Common injury", "DKDiPT_down": "DKD-associated state, down",
             "SharedInjury_down": "Common injury, down"}


def setup():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7, "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.bbox": "tight", "savefig.pad_inches": 0.05, "legend.frameon": True,
        "legend.fancybox": False, "legend.edgecolor": "#CCCCCC", "axes.labelweight": "bold"})


def panel_label(ax, letter: str):
    ax.set_title(letter, fontweight="bold", fontsize=11, loc="left", pad=8)


def source(df: pd.DataFrame, name: str):
    df.to_csv(SRC / f"{name}.tsv", sep="\t", index=False)


def _texts(fig):
    r = fig.canvas.get_renderer()
    out = []
    for ax in fig.axes:
        items = [ax.title, ax._left_title, ax._right_title, ax.xaxis.label, ax.yaxis.label, *ax.texts]
        items += [t for t in ax.get_xticklabels() + ax.get_yticklabels() if t.get_visible()]
        leg = ax.get_legend()
        for t in items:
            if t.get_visible() and t.get_text().strip():
                out.append((ax, t, t.get_window_extent(r)))
        if leg is not None:   # the whole legend frame counts as one text block
            out.append((ax, leg.get_texts()[0] if leg.get_texts() else ax.title, leg.get_window_extent(r)))
    out += [(None, t, t.get_window_extent(r)) for t in fig.texts if t.get_text().strip()]
    return out


def check_overlaps(fig, tol: float = 0.5) -> list[str]:
    """Pairwise text-bbox overlaps and text crossing into another axes' plotting area."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    T = _texts(fig)
    issues = []
    for (a1, t1, b1), (a2, t2, b2) in itertools.combinations(T, 2):
        ov_w = min(b1.x1, b2.x1) - max(b1.x0, b2.x0); ov_h = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
        if ov_w > tol and ov_h > tol:
            issues.append(f"text overlap: '{t1.get_text()[:25]}' x '{t2.get_text()[:25]}'")
    axes_boxes = [(ax, ax.get_window_extent(r)) for ax in fig.axes if ax.get_visible() and ax.axison]
    for ax_t, t, b in T:
        for ax, ab in axes_boxes:
            if ax is ax_t:
                continue
            ov_w = min(b.x1, ab.x1) - max(b.x0, ab.x0); ov_h = min(b.y1, ab.y1) - max(b.y0, ab.y0)
            if ov_w > tol and ov_h > tol:
                issues.append(f"text '{t.get_text()[:25]}' enters another axes")
    return issues


def check_panel_alignment(fig, rows: list[list], tol_px: float = 2.0) -> list[str]:
    """Panel letters of axes placed in the same row must share a baseline (y)."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    issues = []
    for row in rows:
        ys = [ax._left_title.get_window_extent(r).y0 for ax in row]
        if max(ys) - min(ys) > tol_px:
            issues.append(f"panel letters misaligned in row {[ax._left_title.get_text() for ax in row]}: {ys}")
    return issues


def save(fig, name: str, rows: list[list] | None = None):
    issues = check_overlaps(fig) + (check_panel_alignment(fig, rows) if rows else [])
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=600)
    print(name, "layout issues:", len(issues))
    for i in issues:
        print("  ", i)
    return issues


def hedges_ci(g, na, nb):
    import numpy as np
    se = np.sqrt((na + nb) / (na * nb) + g ** 2 / (2 * (na + nb)))
    return g - 1.96 * se, g + 1.96 * se
