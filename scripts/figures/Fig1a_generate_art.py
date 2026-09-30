"""Fig1a art generation via an OpenAI-compatible image endpoint (text-free illustrations).

Endpoint and key are read from env IMAGE_API_BASE / IMAGE_API_KEY (never stored). Model / prompts are logged to
docs/fig1a_prompts.md. Usage: IMAGE_API_BASE=... IMAGE_API_KEY=... python Fig1a_generate_art.py <name> <prompt_key> [model]"""
from __future__ import annotations
import base64, datetime, json, os, sys, urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = os.environ["IMAGE_API_BASE"]   # OpenAI-compatible endpoint, supplied at run time
STYLE = ("Clean BioRender-style scientific illustration, flat vector style, white background, thin dark-grey outlines, "
         "soft medical palette: muted red #C44E52, muted blue #4C72B0, light grey #BBBBBB, pale teal and pale pink accents. "
         "ABSOLUTELY NO TEXT, NO LETTERS, NO NUMBERS, NO LABELS, NO WATERMARK anywhere in the image.")
PROMPTS = {
    "workflow": STYLE + " A left-to-right three-stage research workflow on one wide canvas: stage 1 a human kidney with a "
              "magnified nephron and proximal tubule epithelial cells and single-nucleus droplets; stage 2 a small neural network "
              "with two separated latent spaces (a padlock symbol meaning frozen); stage 3 validation icons: a spatial "
              "transcriptomics slide grid, a haematoxylin-eosin tissue patch, a simple clinical line chart, a microarray chip. "
              "Stages separated by large empty white gaps with simple arrows; generous empty space above each stage.",
    "icons": STYLE + " An icon sheet: eight separate small icons arranged in a regular 4 columns x 2 rows grid, each icon centred "
              "in its own cell with lots of white space between cells, no frames: (1) human kidney, (2) proximal tubule segment "
              "with cuboidal epithelial cells, (3) cluster of single nuclei droplets, (4) neural network with two separated "
              "latent clouds, (5) padlock, (6) spatial transcriptomics slide with a spot grid, (7) pink-purple haematoxylin-eosin "
              "tissue section patch with tubules, (8) clinical chart with a rising line and a urine test tube.",
}


def main():
    name, key = sys.argv[1], sys.argv[2]
    model = sys.argv[3] if len(sys.argv) > 3 else "gpt-image-2"
    body = {"model": model, "prompt": PROMPTS[key], "size": "1536x1024", "n": 1}
    req = urllib.request.Request(f"{BASE}/images/generations", data=json.dumps(body).encode(),
                                 headers={"Authorization": f"Bearer {os.environ['IMAGE_API_KEY']}",
                                          "Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=600) as r:
        d = json.load(r)
    item = d["data"][0]
    img = base64.b64decode(item["b64_json"]) if item.get("b64_json") else urllib.request.urlopen(item["url"]).read()
    out = ROOT / "figures/Fig1a_raw" / f"{name}.png"
    out.write_bytes(img)
    log = ROOT / "docs/fig1a_prompts.md"
    if not log.exists():
        log.write_text("# Fig1a image-generation log\n\nEndpoint: OpenAI-compatible images/generations. Key supplied via environment, not stored.\n\n")
    with open(log, "a") as fh:
        fh.write(f"## {name}\n- date: {datetime.date.today()}\n- model: {model}\n- size: 1536x1024\n- prompt ({key}): {PROMPTS[key]}\n- output: figures/Fig1a_raw/{name}.png\n\n")
    print(out, len(img))


if __name__ == "__main__":
    main()
