"""10a_visium_gse211785: tissue localisation of the locked programmes on GSE211785 Visium (FFPE)
sections with same-section H&E.

Several Visium donors also contributed discovery sc/snRNA (see 01_cohort_audit), so this stage is
reported as tissue localisation, not as independent replication. Spot pixel coordinates come from
the Loupe alignment json of each section (oligo row/col/imageX/imageY) joined to spot barcodes via
the Visium v1 barcode -> array map; a section-json assignment is accepted only if the json's
in-tissue spot count equals the section's spot count in the author metadata.

Outputs results/10a_visium_gse211785/: spot_scores.tsv.gz, set_coverage.tsv, ipt_markers.tsv,
section_table.tsv, celltype_scores.tsv, section_summary.tsv, spatial_association.tsv,
diagnosis_contrast.tsv
"""
from __future__ import annotations
import gzip, importlib, json, re, sys, tarfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance
from stats import compare
vc = importlib.import_module("10_visium_common")

STAGE = "10a_visium_gse211785"
RAW = ROOT / "data/raw/GSE211785"


def _key(name: str) -> str:
    n = name.replace("_", "-").replace("V11Y24-79-", "V11Y24-079-")
    return re.search(r"(V\d+[A-Z]\d+-\d+-[A-D]1)", n).group(1)


def loupe_jsons() -> dict[str, dict]:
    """All alignment jsons (series level + RAW.tar) keyed by slide-area, with the H&E tif of the
    same slide-area and the donor id carried by the tif / json file name."""
    out = {}
    for p in RAW.glob("GSE211785_V*.json.gz"):
        out[_key(p.name)] = {"json": json.load(gzip.open(p)), "src": p.name}
    tifs = {_key(p.name): p.name for p in RAW.glob("*.tif.gz")}
    with tarfile.open(RAW / "GSE211785_RAW.tar") as t:
        for m in t.getmembers():
            if m.name.endswith(".json.gz"):
                out[_key(m.name)] = {"json": json.loads(gzip.decompress(t.extractfile(m).read())), "src": m.name}
            elif m.name.endswith(".tif.gz"):
                tifs[_key(m.name)] = m.name
    for k, v in out.items():
        v["tif"] = tifs.get(k)
        ids = re.findall(r"HK\d+", (v["tif"] or "") + v["src"])
        v["file_donor"] = ids[0] if ids else None
        v["n_tissue"] = sum(1 for o in v["json"]["oligo"] if o.get("tissue"))
    return out


def main():
    cfg = load_config(); seed = cfg["seed"]; set_global_seed(seed); V = cfg["visium"]
    out = stage_dir(STAGE)
    proc = ROOT / cfg["paths"]["processed"]
    M = scipy.io.mmread(proc / "GSE211785_ST_counts.mtx").tocsr().T.tocsr()      # spots x genes
    genes = pd.Index(open(proc / "GSE211785_ST_genes.txt").read().split())
    bcs = pd.Index(open(proc / "GSE211785_ST_barcodes.txt").read().split())
    meta = pd.read_csv(RAW / "GSE211785_ST_metadata.txt.gz", sep="\t", index_col=0).reindex(bcs)
    meta["donor"] = meta["orig.ident"].str.extract(r"(HK\d+)")[0]
    meta["raw_barcode"] = [b.split("_")[0] for b in bcs]
    sec = pd.read_csv(ROOT / "results/01_cohort_audit/sample_table_GSE211785_visium.tsv", sep="\t").set_index("donor")
    meta["diagnosis"] = meta.donor.map(sec.diagnosis_geo)

    # --- pixel coordinates ---
    amap = vc.barcode_array_map()
    meta = meta.join(amap, on="raw_barcode")
    js = loupe_jsons()
    n_spots = meta.groupby("donor").size()
    assign, srows = {}, []
    for k, v in js.items():
        # file-named donor: accept if spot counts agree within 1% (author QC removed a few spots);
        # otherwise (one json without image) require an exact spot-count match
        if v["file_donor"]:
            d = v["file_donor"]
            ok = [d] if d in n_spots.index and abs(n_spots[d] - v["n_tissue"]) <= 0.01 * n_spots[d] else []
            how = "file name + spot count (within 1%)"
        else:
            ok = [d for d in n_spots.index if n_spots[d] == v["n_tissue"] and d not in
                  {x["file_donor"] for x in js.values()}]
            how = "exact spot count"
        srows.append({"slide_area": k, "json": v["src"], "tif": v["tif"], "n_tissue_json": v["n_tissue"],
                      "n_spots_metadata": n_spots.get(ok[0]) if len(ok) == 1 else None,
                      "donor": ok[0] if len(ok) == 1 else None, "assignment": how if len(ok) == 1 else "unassigned"})
        if len(ok) == 1:
            assign[ok[0]] = v
    st = pd.DataFrame(srows)
    meta["pxl_row_fullres"] = np.nan; meta["pxl_col_fullres"] = np.nan; meta["in_tissue_json"] = np.nan
    for d, v in assign.items():
        ol = pd.DataFrame(v["json"]["oligo"]).set_index(["row", "col"])
        idx = meta.index[meta.donor == d]
        key = list(zip(meta.loc[idx, "array_row"], meta.loc[idx, "array_col"]))
        sub = ol.reindex(key)
        meta.loc[idx, "pxl_row_fullres"] = sub.imageY.values
        meta.loc[idx, "pxl_col_fullres"] = sub.imageX.values
        meta.loc[idx, "in_tissue_json"] = sub["tissue"].fillna(False).astype(float).values
    st["frac_spots_in_tissue_json"] = st.donor.map(meta.groupby("donor").in_tissue_json.mean())
    st = st.merge(sec[["diagnosis_geo", "he_image"]].rename(columns={"diagnosis_geo": "diagnosis"}),
                  left_on="donor", right_index=True, how="left")
    st.to_csv(out / "section_table.tsv", sep="\t", index=False)

    # --- scores ---
    mk = vc.ipt_markers(cfg); mk.to_csv(out / "ipt_markers.tsv", sep="\t", index=False)
    sets = {k: v for k, v in vc.locked_sets().items() if k in vc.SETS_USED}
    sets["iPT_marker"] = mk[mk.set == "iPT_marker"].gene.tolist()
    sets["healthyPT_marker"] = mk[mk.set == "healthyPT_marker"].gene.tolist()
    Z, _ = vc.lognorm_z(M, genes); Z.index = bcs
    P, _ = vc.reference_signatures(cfg, V["deconv_markers"], V["deconv_min_cp10k"])
    P = P.loc[P.index.intersection(genes)]
    P, dmk = vc.reference_signatures_from(P, V["deconv_markers"], V["deconv_min_cp10k"])
    dmk.to_csv(out / "deconvolution_markers.tsv", sep="\t", index=False)
    W, _ = vc.nnls_deconvolve(M, genes, P, dmk); W.index = bcs; W.columns = [f"prop_{c}" for c in W.columns]
    S, cov = vc.set_scores(Z, sets, V["min_genes"])
    cov.to_csv(out / "set_coverage.tsv", sep="\t", index=False)
    evaluable = cov.set[cov.evaluable].tolist()
    spot = meta[["orig.ident", "donor", "diagnosis", "celltype", "raw_barcode", "pxl_row_fullres",
                 "pxl_col_fullres", "in_tissue_json"]].join(S).join(W)
    ptm = vc.pt_lineage_markers(P, V["deconv_markers"], V["deconv_min_cp10k"])
    pd.Series(ptm, name="gene").to_csv(out / "pt_lineage_markers.tsv", sep="\t", index=False)
    spot["PT_lineage_score"] = Z[[g for g in ptm if g in Z.columns]].mean(1)
    spot["PT_share"] = spot.prop_PT_healthy + spot.prop_iPT
    q = 1 - V["pt_rich_top_frac"]
    spot["PT_rich"] = spot.groupby("orig.ident").PT_lineage_score.transform(lambda x: x >= x.quantile(q))
    spot["injury_index"] = spot.iPT_marker - spot.healthyPT_marker
    spot = spot.rename(columns={"orig.ident": "section", "raw_barcode": "barcode_raw", "diagnosis": "group"})
    spot.index.name = "barcode"

    # --- localisation by author spot label (section means, then across-section mean) ---
    cs = spot.groupby(["section", "celltype"])[evaluable].mean()
    ct = cs.groupby("celltype").agg(["mean", "count"])
    ct.columns = [f"{a}_{b}" for a, b in ct.columns]
    ct.to_csv(out / "celltype_scores.tsv", sep="\t")
    # --- spatial association within each section ---
    ar = []
    for s, d0 in spot.groupby("section"):
        d = d0[d0.PT_rich]
        if len(d) < V["min_spots_section"]:
            continue
        for k in [x for x in evaluable if x.startswith(("DKDiPT", "DEonly"))]:
            ar.append({"section": s, "group": d.group.iloc[0], "set": k, "n_PT_rich_spots": len(d),
                       "r_with_injury_index": stats.spearmanr(d[k], d.injury_index).correlation,
                       "partial_r_injury_given_SharedInjury_up": vc.partial_r(d[k], d.injury_index, d.SharedInjury_up),
                       "r_with_SharedInjury_up": stats.spearmanr(d[k], d.SharedInjury_up).correlation,
                       "r_with_author_iPT_label": stats.spearmanr(d[k], d.celltype.eq("iPT")).correlation})
    A = pd.DataFrame(ar); A.to_csv(out / "spatial_association.tsv", sep="\t", index=False)
    # --- section-level summaries (unit = donor section) ---
    from sklearn.metrics import roc_auc_score
    ptl = spot.celltype.isin(["iPT", "PT_S1", "PT_S2", "PT_S3"])
    pc = [{"check": "PT_lineage_score: PT/iPT-labelled vs other spots", "auc": roc_auc_score(ptl, spot.PT_lineage_score)},
          {"check": "NNLS PT_share: PT/iPT-labelled vs other spots", "auc": roc_auc_score(ptl, spot.PT_share)},
          {"check": "prop_iPT: iPT vs PT_S1-3 spots", "auc": roc_auc_score(spot.celltype[ptl] == "iPT", spot.prop_iPT[ptl])},
          {"check": "injury_index: iPT vs PT_S1-3 spots", "auc": roc_auc_score(spot.celltype[ptl] == "iPT", spot.injury_index[ptl])},
          {"check": "SharedInjury_up: iPT vs PT_S1-3 spots", "auc": roc_auc_score(spot.celltype[ptl] == "iPT", spot.SharedInjury_up[ptl])},
          {"check": "fraction of PT_rich spots with PT/iPT author label", "auc": float(ptl[spot.PT_rich].mean())},
          {"check": "fraction of other spots with PT/iPT author label", "auc": float(ptl[~spot.PT_rich].mean())}]
    pd.DataFrame(pc).to_csv(out / "deconvolution_positive_control.tsv", sep="\t", index=False)
    spot.to_csv(out / "spot_scores.tsv.gz", sep="\t")
    spot["iPT_rich"] = spot.PT_rich            # section summaries use PT-rich spots
    rows = []
    for s, d in spot.groupby("section"):
        r = {"section": s, "donor": d.donor.iloc[0], "group": d.group.iloc[0], "n_spots": len(d),
             "frac_author_iPT": (d.celltype == "iPT").mean()}
        for k in evaluable:
            r[f"{k}_iPTrich"] = d.loc[d.iPT_rich, k].mean(); r[f"{k}_all"] = d[k].mean()
            r[f"{k}_authoriPT"] = d.loc[d.celltype == "iPT", k].mean()
        rows.append(r)
    SS = pd.DataFrame(rows); SS.to_csv(out / "section_summary.tsv", sep="\t", index=False)
    SS["grp2"] = np.where(SS.group == "DKD", "DKD", np.where(SS.group.isin(["HKD", "CKD"]), "nonDM_CKD", SS.group))
    dc = []
    for k in [c for c in SS.columns if c.endswith(("_iPTrich", "_authoriPT"))]:
        for a, b in [("DKD", "nonDM_CKD"), ("DKD", "Control"), ("nonDM_CKD", "Control")]:
            dc.append({"metric": k, "contrast": f"{a}-{b}", **compare(SS[k], SS.grp2, a, b)})
    A["grp2"] = np.where(A.group == "DKD", "DKD", np.where(A.group.isin(["HKD", "CKD"]), "nonDM_CKD", A.group))
    for k in ["r_with_injury_index", "partial_r_injury_given_SharedInjury_up"]:
        for a, b in [("DKD", "nonDM_CKD"), ("DKD", "Control")]:
            dc.append({"metric": f"DKDiPT_up_{k}_PTrich", "contrast": f"{a}-{b}", **compare(A[k], A.grp2, a, b)})
    DC = pd.DataFrame(dc); DC.to_csv(out / "diagnosis_contrast.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(st.to_string()); print(cov[["set", "n_set", "n_measured", "evaluable"]].to_string())
    print(ct.round(3).to_string())
    print(A.groupby("set")[["r_with_injury_index", "partial_r_injury_given_SharedInjury_up", "r_with_SharedInjury_up", "r_with_author_iPT_label"]].median().round(3))
    print(pd.DataFrame(pc).round(3).to_string())
    print(DC.round(3).to_string())
    write_provenance(STAGE, [proc / "GSE211785_ST_counts.mtx", RAW / "GSE211785_ST_metadata.txt.gz", RAW / "GSE211785_RAW.tar",
                             ROOT / "results/05_lock_programme/gene_sets.tsv", proc / "pb_fine.h5ad"],
                     sorted(out.glob("*.tsv")) + [out / "spot_scores.tsv.gz"], seed, extra={"visium": V})


if __name__ == "__main__":
    main()
