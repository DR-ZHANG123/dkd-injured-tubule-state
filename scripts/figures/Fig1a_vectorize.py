"""Fig1a: crop the text-free icon sheet into single icons and trace each to colour SVG (vtracer).

Input  figures/Fig1a_raw/icons_v2.png (4 x 2 icon grid; see docs/fig1a_prompts.md)
Output figures/Fig1a_art/icon_<name>.svg, figures/Fig1a_art.svg (whole traced sheet),
       figures/source_data/Fig1a_vectorize_params.tsv"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage
import vtracer

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "figures/Fig1a_raw/icons_v2.png"
OUT = ROOT / "figures/Fig1a_art"
NAMES = [["kidney", "tubule", "nuclei", "network"], ["padlock", "slide", "he", "clinical"]]
PARAMS = dict(colormode="color", hierarchical="stacked", mode="spline", filter_speckle=6, color_precision=6,
              layer_difference=16, corner_threshold=60, length_threshold=4.0, max_iterations=10,
              splice_threshold=45, path_precision=3)
UPSCALE = 2


def crop_cells(img: np.ndarray):
    h, w = img.shape[:2]
    for r in range(2):
        for c in range(4):
            cell = img[r * h // 2:(r + 1) * h // 2, c * w // 4:(c + 1) * w // 4]
            ink = (cell[..., :3].min(2) < 235)
            # keep the ink cluster(s) of this cell's icon; fragments of neighbouring icons that
            # touch the cell's left/right edge are dropped
            lab, n = ndimage.label(ndimage.binary_dilation(ink, iterations=12))
            keep = np.zeros_like(ink)
            for k in range(1, n + 1):
                comp = lab == k
                cols = np.where(comp.any(0))[0]
                touches = cols.min() == 0 or cols.max() == comp.shape[1] - 1
                if not touches or comp.sum() > 0.3 * (lab > 0).sum():
                    keep |= comp
            ink &= keep
            ys, xs = np.where(ink)
            pad = 6
            y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad, cell.shape[0])
            x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad, cell.shape[1])
            crop = cell[y0:y1, x0:x1].copy()
            crop[~ink[y0:y1, x0:x1] & (crop.min(2) >= 200)] = 255   # blank stray light pixels outside the icon
            yield NAMES[r][c], crop


def transparent_background(rgb: np.ndarray, thr: int = 238) -> np.ndarray:
    """RGBA with alpha 0 on near-white regions connected to the crop border (interior whites kept)."""
    white = rgb.min(2) >= thr
    lab, _ = ndimage.label(white)
    border = set(np.unique(np.r_[lab[0], lab[-1], lab[:, 0], lab[:, -1]])) - {0}
    bg = np.isin(lab, list(border))
    return np.dstack([rgb, np.where(bg, 0, 255).astype(np.uint8)])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    img = np.asarray(Image.open(SRC).convert("RGB"))
    rows = []
    for name, crop in crop_cells(img):
        im = Image.fromarray(transparent_background(crop), "RGBA")
        im = im.resize((im.width * UPSCALE, im.height * UPSCALE), Image.LANCZOS)
        p = OUT / f"icon_{name}.png"; im.save(p)
        vtracer.convert_image_to_svg_py(str(p), str(OUT / f"icon_{name}.svg"), **PARAMS)
        p.unlink()
        rows.append({"icon": name, "crop_w": crop.shape[1], "crop_h": crop.shape[0],
                     "svg_bytes": (OUT / f"icon_{name}.svg").stat().st_size})
    sheet = ROOT / "figures/Fig1a_art/_sheet.png"
    Image.fromarray(transparent_background(img), "RGBA").save(sheet)
    vtracer.convert_image_to_svg_py(str(sheet), str(ROOT / "figures/Fig1a_art.svg"), **PARAMS)
    sheet.unlink()
    t = pd.DataFrame(rows)
    for k, v in {**PARAMS, "upscale": UPSCALE, "source": SRC.name}.items():
        t[k] = v
    t.to_csv(ROOT / "figures/source_data/Fig1a_vectorize_params.tsv", sep="\t", index=False)
    print(t[["icon", "crop_w", "crop_h", "svg_bytes"]])


if __name__ == "__main__":
    main()
