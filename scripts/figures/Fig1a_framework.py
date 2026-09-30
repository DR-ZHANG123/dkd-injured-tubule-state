"""Fig1a study framework (vector): traced icons (Fig1a_vectorize.py) + editable text labels.

draw_framework(ax) draws into any axes (used by Fig1_cohorts.py for panel a); running this file
saves the stand-alone panel figures/Fig1a_framework.{svg,pdf,png}. Numbers in the labels come
from the manuscript / results tables (see figures/source_data/Fig1a_labels.tsv)."""
from __future__ import annotations
import re
from functools import lru_cache
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, PathPatch, FancyArrowPatch
from matplotlib.transforms import Affine2D
from svgpath2mpl import parse_path
from _style import FIG, setup, source

ART = FIG / "Fig1a_art"
W, H = 100.0, 31.0                      # panel coordinate frame (equal aspect)
STAGE_COL = {"disc": "#EEF2F8", "lock": "#F4F0F8", "val": "#F3F6F1"}
EDGE = "#9AA3AE"
RED, BLUE, GREEN, GREY = "#C44E52", "#4C72B0", "#3E8E57", "#6E6E6E"
PATH_RE = re.compile(r'<path d="([^"]+)" fill="([^"]+)"(?: fill-opacity="([^"]+)")? transform="translate\(([-\d.]+),([-\d.]+)\)"')


@lru_cache(maxsize=None)
def load_icon(name: str):
    s = (ART / f"icon_{name}.svg").read_text()
    w = float(re.search(r'width="([\d.]+)"', s).group(1)); h = float(re.search(r'height="([\d.]+)"', s).group(1))
    paths = []
    for d, fill, op, tx, ty in PATH_RE.findall(s):
        p = parse_path(d).transformed(Affine2D().translate(float(tx), float(ty)))
        paths.append((p, fill, float(op) if op else 1.0))
    return paths, w, h


def draw_icon(ax, name: str, cx: float, cy: float, height: float):
    """Place icon centred at (cx, cy) with the given height in panel units; returns its extent."""
    paths, w, h = load_icon(name)
    s = height / h
    tr = Affine2D().scale(s, -s).translate(cx - w * s / 2, cy + height / 2)
    for p, fill, op in paths:
        ax.add_patch(PathPatch(p.transformed(tr), facecolor=fill, edgecolor="none", alpha=op, zorder=3))
    return (cx - w * s / 2, cx + w * s / 2, cy - height / 2, cy + height / 2)


def box(ax, x0, y0, x1, y1, fc, ec=EDGE, lw=0.6, r=1.2, z=1):
    ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0, boxstyle=f"round,pad=0,rounding_size={r}",
                                fc=fc, ec=ec, lw=lw, zorder=z))
    return (x0, x1, y0, y1)


def arrow(ax, x0, x1, y):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>,head_length=3,head_width=2",
                                 color="#5B6470", lw=1.0, zorder=4, shrinkA=0, shrinkB=0))


LABELS = [
    # key, box, text, fontsize, weight, colour, y (centre)
    ("disc_title", "disc", "Discovery atlas (GSE211785)", 7, "bold", "#2B2B2B", 28.3),
    ("disc_1", "disc", "16 control · 13 HKD · 7 DKD donors", 6, "normal", "#2B2B2B", 9.6),
    ("disc_2", "disc", "Technically robust single-nucleus contrasts", 6, "normal", "#2B2B2B", 6.4),
    ("disc_3", "disc", "iPT expansion common to DKD and HKD", 6, "normal", "#2B2B2B", 3.2),
    ("lock_title", "lock", "DKD-associated state", 7, "bold", "#2B2B2B", 28.3),
    ("lock_1", "lock", "scDisInFact + donor pseudobulk", 6, "normal", "#2B2B2B", 9.6),
    ("lock_2", "lock", "DKD vs HKD within iPT", 6, "normal", "#2B2B2B", 6.4),
    ("lock_3", "lock", "12 marker genes", 6, "bold", RED, 3.2),
    ("val_title", "val", "Independent patients and tissue", 7, "bold", "#2B2B2B", 28.3),
    ("val_rep_t", "rep", "Present in KPMP", 6.5, "bold", GREEN, 10.2),
    ("val_rep_1", "rep", "KPMP snRNA (40 DKD, 16 HKD)", 5.6, "normal", "#2B2B2B", 7.6),
    ("val_rep_2", "rep", "79% of genes same direction", 5.6, "normal", "#2B2B2B", 5.3),
    ("val_rep_3", "rep", "Accompanies albuminuria", 5.6, "normal", "#2B2B2B", 3.0),
    ("val_nd_t", "nd", "Confined to injured cells", 6.5, "bold", GREY, 10.2),
    ("val_nd_1", "nd", "Tissue and histology", 5.6, "normal", "#2B2B2B", 7.6),
    ("val_nd_2", "nd", "reflect common injury", 5.6, "normal", "#2B2B2B", 5.3),
    ("val_nd_3", "nd", "State resolved in iPT cells", 5.6, "normal", "#2B2B2B", 3.0),
]


def draw_framework(ax, scale: float = 1.0):
    """scale = panel width / 7.2 in; font sizes follow so labels keep their fit."""
    ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    boxes = {"disc": box(ax, 0.4, 0.8, 34.0, 30.4, STAGE_COL["disc"]),
             "lock": box(ax, 37.6, 0.8, 58.4, 30.4, STAGE_COL["lock"]),
             "val": box(ax, 62.0, 0.8, 99.6, 30.4, STAGE_COL["val"])}
    boxes["rep"] = box(ax, 63.0, 1.6, 80.6, 12.0, "white", ec=GREEN, lw=0.6, r=0.8, z=2)
    boxes["nd"] = box(ax, 81.4, 1.6, 98.8, 12.0, "white", ec="#AAAAAA", lw=0.6, r=0.8, z=2)
    arrow(ax, 34.4, 37.2, 18.5); arrow(ax, 58.8, 61.6, 18.5)
    icons = {"kidney": draw_icon(ax, "kidney", 6.2, 18.5, 11.0), "tubule": draw_icon(ax, "tubule", 16.8, 18.5, 9.5),
             "nuclei": draw_icon(ax, "nuclei", 27.6, 18.5, 8.0), "network": draw_icon(ax, "network", 45.4, 18.8, 8.6),
             "padlock": draw_icon(ax, "padlock", 54.4, 18.8, 7.0), "slide": draw_icon(ax, "slide", 69.4, 18.5, 7.0),
             "he": draw_icon(ax, "he", 81.2, 18.5, 8.6), "clinical": draw_icon(ax, "clinical", 92.4, 18.5, 8.6)}
    texts = {}
    for key, b, txt, fs, fw, col, y in LABELS:
        x0, x1, _, _ = boxes[b]
        texts[key] = (ax.text((x0 + x1) / 2, y, txt, ha="center", va="center", fontsize=fs * scale, fontweight=fw,
                              color=col, zorder=5), b)
    return boxes, icons, texts


def check_layout(fig, ax, boxes, icons, texts, tol=0.5) -> list[str]:
    """Each label inside its box, horizontally centred, not overlapping other labels or icons."""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); T = ax.transData
    disp = lambda x0, x1, y0, y1: (*T.transform((x0, y0)), *T.transform((x1, y1)))
    issues = []
    bb = {k: t.get_window_extent(r) for k, (t, _) in texts.items()}
    for k, (t, b) in texts.items():
        X0, Y0, X1, Y1 = disp(*boxes[b]); e = bb[k]
        if e.x0 < X0 + tol or e.x1 > X1 - tol or e.y0 < Y0 or e.y1 > Y1:
            issues.append(f"label '{k}' exits box '{b}'")
        if abs((e.x0 + e.x1) / 2 - (X0 + X1) / 2) > 1.0:
            issues.append(f"label '{k}' not centred")
        for name, ext in icons.items():
            I0, J0, I1, J1 = disp(*ext)
            if min(e.x1, I1) - max(e.x0, I0) > tol and min(e.y1, J1) - max(e.y0, J0) > tol:
                issues.append(f"label '{k}' overlaps icon '{name}'")
    keys = list(bb)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, c = bb[keys[i]], bb[keys[j]]
            if min(a.x1, c.x1) - max(a.x0, c.x0) > tol and min(a.y1, c.y1) - max(a.y0, c.y0) > tol:
                issues.append(f"labels overlap: {keys[i]} x {keys[j]}")
    for name, ext in icons.items():                           # icons stay inside their stage box
        stage = "disc" if ext[1] < 34.5 else "lock" if ext[1] < 58.5 else "val"
        x0, x1, y0, y1 = boxes[stage]
        if ext[0] < x0 or ext[1] > x1 or ext[2] < y0 or ext[3] > y1:
            issues.append(f"icon '{name}' exits stage box")
    return issues


def main():
    setup()
    fig, ax = plt.subplots(figsize=(7.2, 7.2 * H / W))
    fig.subplots_adjust(0, 0, 1, 1)
    boxes, icons, texts = draw_framework(ax)
    issues = check_layout(fig, ax, boxes, icons, texts)
    for ext in ["svg", "pdf"]:
        fig.savefig(FIG / f"Fig1a_framework.{ext}")
    fig.savefig(FIG / "Fig1a_framework.png", dpi=600)
    source(pd.DataFrame([{"key": k, "box": b, "text": t, "fontsize": fs} for k, b, t, fs, *_ in LABELS]), "Fig1a_labels")
    print("Fig1a layout issues:", len(issues)); [print("  ", i) for i in issues]


if __name__ == "__main__":
    main()
