"""01_cohort_audit: donor x diagnosis x modality pairing table for every data source.

Outputs (results/01_cohort_audit/):
  donor_modality_matrix.tsv   one row per donor, one column per modality
  sample_table_<source>.tsv   one row per library / section / array
  donor_overlap.tsv           donors present in more than one source
  audit_summary.json          counts per source x diagnosis, flags
"""
from __future__ import annotations
import gzip, json, re, sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
from repro import ROOT, load_config, set_global_seed, stage_dir, write_provenance

STAGE = "01_cohort_audit"


def series_matrix(path: Path) -> pd.DataFrame:
    """Parse !Sample_ rows of a GEO series matrix; repeated keys get _1, _2 ..."""
    d: dict[str, list[str]] = {}
    for line in gzip.open(path, "rt"):
        if not line.startswith("!Sample_"):
            continue
        row = line.rstrip("\n").split("\t")
        key, k = row[0], 0
        while key in d:
            k += 1
            key = f"{row[0]}_{k}"
        d[key] = [v.strip('"') for v in row[1:]]
    return pd.DataFrame(d)


def characteristics(df: pd.DataFrame) -> pd.DataFrame:
    """Turn 'key: value' characteristic columns into named columns."""
    out = {}
    for c in [c for c in df.columns if c.startswith("!Sample_characteristics")]:
        for i, v in df[c].items():
            if ": " in v:
                k, val = v.split(": ", 1)
                out.setdefault(k.strip().lower(), {})[i] = val.strip()
    return pd.DataFrame(out, index=df.index)


def donor_id(x: str) -> str:
    m = re.match(r"(HK\d+)", str(x))
    return m.group(1) if m else str(x)


def audit_gse211785(raw: Path):
    m = pd.read_csv(raw / "GSE211785/GSE211785_scRNA-seq_snRNA-seq_snATAC-seq_metadata.txt.gz",
                    index_col=0, low_memory=False)
    lib = (m.groupby("orig_ident")
             .agg(tech=("tech", "first"), sex=("sex", "first"), age=("age", "first"),
                  diagnosis=("group", "first"), n_cells=("tech", "size"))
             .reset_index().rename(columns={"orig_ident": "library"}))
    lib["donor"] = lib.library.map(donor_id)
    lib["source"] = "GSE211785"
    # ST sections: diagnosis from series matrix; donor id from image file names / ST metadata
    sm = series_matrix(raw / "GSE211785/GSE211785_series_matrix.txt.gz")
    ch = characteristics(sm)
    sm = pd.concat([sm, ch], axis=1)
    st = sm[sm["!Sample_title"].str.endswith("_ST")].copy()
    sup = [c for c in st.columns if "supplementary_file" in c]
    st["donor"] = st[sup].apply(lambda r: next((donor_id(re.search(r"HK\d+", v).group(0))
                                                for v in r if re.search(r"HK\d+", v)), None), axis=1)
    stmeta = pd.read_csv(raw / "GSE211785/GSE211785_ST_metadata.txt.gz", sep="\t", low_memory=False)
    sections = (stmeta.groupby("orig.ident").agg(status=("Status", "first"), n_spots=("Status", "size"))
                .reset_index())
    sections["donor"] = sections["orig.ident"].map(donor_id)
    # H&E images shipped at series level or in RAW.tar (names listed in filelist.txt)
    names = [p.name for p in (raw / "GSE211785").glob("*.tif.gz")]
    fl = raw / "GSE211785/filelist.txt"
    if fl.exists():
        names += [l.split("\t")[1] for l in fl.read_text().splitlines()[1:] if l.split("\t")[1].endswith(".tif.gz")]
    img_donors = {donor_id(re.search(r"HK\d+", n).group(0)) for n in names if re.search(r"HK\d+", n)}
    sections["he_image"] = sections.donor.isin(img_donors)
    # diagnosis for ST sections: GEO rows whose file names carry a donor id are exact; the rest
    # are matched on (age, sex, diagnosis) against donor tables of the sc/sn metadata and of the
    # Zenodo atlas, then by elimination on the Control/Disease status of the ST metadata.
    dx = dict(zip(st.donor.dropna(), st.loc[st.donor.notna(), "disease"]))
    z = audit_zenodo(raw).drop_duplicates("donor")
    ref = pd.concat([lib.drop_duplicates("donor")[["donor", "age", "sex", "diagnosis"]],
                     z.rename(columns={"Age": "age", "Sex": "sex"})[["donor", "age", "sex", "diagnosis"]]])
    ref = ref[ref.donor.isin(sections.donor)].drop_duplicates("donor")
    pending = st[st.donor.isna()]
    for i, r in pending.iterrows():
        hit = ref[(ref.age.astype(float) == float(r["age"])) & (ref.sex.str[0] == r["gender"][0])
                  & (ref.diagnosis == r["disease"]) & ~ref.donor.isin(dx)]
        if len(hit) == 1:
            dx[hit.donor.iloc[0]] = r["disease"]; st.loc[i, "donor"] = hit.donor.iloc[0]
    status = dict(zip(sections.donor, sections.status))
    left = st[st.donor.isna()]
    for lab in ["Control", "Disease"]:
        rows = left[(left.disease == "Control") == (lab == "Control")]
        free = [d for d in sections.donor if d not in dx and status[d] == lab]
        if len(rows) >= 1 and len(set(rows.disease)) == 1 and len(free) == len(rows):
            for d in free:
                dx[d] = rows.disease.iloc[0]
    sections["diagnosis_geo"] = sections.donor.map(dx)
    sections["source"] = "GSE211785_Visium"
    return lib, sections, st


def audit_zenodo(raw: Path) -> pd.DataFrame:
    z = pd.read_excel(raw / "zenodo_dkd_spatial/Diagnosis.xlsx")
    z.columns = [c.strip() for c in z.columns]
    z["platform"] = z["Sample ID"].str.extract(r"_(CosMx|Xenium)$")[0]
    z["donor"] = z["Sample ID"].str.replace(r"(_2)?_(CosMx|Xenium)$", "", regex=True).map(donor_id)
    z["source"] = "Zenodo19868428"
    return z.rename(columns={"Condition": "diagnosis"})


def audit_gse195460(raw: Path) -> pd.DataFrame:
    sm = series_matrix(raw / "GSE195460/GSE195460_series_matrix.txt.gz")
    ch = characteristics(sm)
    df = pd.concat([sm[["!Sample_title", "!Sample_geo_accession"]], ch], axis=1)
    df["library"] = df["!Sample_title"].str.replace(r"_sn(RNA|ATAC)seq$", "", regex=True)
    df["assay"] = df["!Sample_title"].str.extract(r"_sn(RNA|ATAC)seq$")[0]
    rna_files = sorted((raw / "GSE195460").glob("*_filtered_feature_bc_matrix.h5"))
    rna = pd.DataFrame({"library": [re.search(r"_(Control\d|DN\d)_", p.name).group(1) for p in rna_files],
                        "rna_file": [p.name for p in rna_files]})
    clin = df.drop_duplicates("library").drop(columns=["assay"])
    out = rna.merge(clin, on="library", how="left")
    out["diagnosis"] = out.library.str.replace(r"\d+$", "", regex=True).map({"Control": "Control", "DN": "DKD"})
    out["source"] = "GSE195460"
    return out


def audit_gse104954(raw: Path) -> pd.DataFrame:
    rows = []
    for f in sorted((raw / "GSE104954").glob("*series_matrix.txt.gz")):
        sm = series_matrix(f)
        ch = characteristics(sm)
        d = pd.DataFrame({"sample": sm["!Sample_geo_accession"], "title": sm["!Sample_title"],
                          "platform": re.search(r"GPL\d+", f.name).group(0),
                          "diagnosis": ch.get("diagnosis", pd.Series(index=sm.index, dtype=str))})
        # living donors carry no diagnosis field; their titles contain "-LD"
        d.loc[d.title.str.contains("-LD") & d.diagnosis.isna(), "diagnosis"] = "Living donor"
        rows.append(d)
    out = pd.concat(rows, ignore_index=True)
    out["source"] = "GSE104954"
    return out


def main():
    cfg = load_config()
    set_global_seed(cfg["seed"])
    raw = ROOT / cfg["paths"]["raw"]
    out = stage_dir(STAGE)

    lib, sections, _ = audit_gse211785(raw)
    zen = audit_zenodo(raw)
    val = audit_gse195460(raw)
    ercb = audit_gse104954(raw)

    lib.to_csv(out / "sample_table_GSE211785_sc_sn.tsv", sep="\t", index=False)
    sections.to_csv(out / "sample_table_GSE211785_visium.tsv", sep="\t", index=False)
    zen.to_csv(out / "sample_table_zenodo_spatial.tsv", sep="\t", index=False)
    val.to_csv(out / "sample_table_GSE195460.tsv", sep="\t", index=False)
    ercb.to_csv(out / "sample_table_GSE104954.tsv", sep="\t", index=False)

    # donor x modality matrix (GSE211785 donors + Zenodo donors that carry HK ids)
    rna = lib[lib.tech.isin(["SC_RNA", "SN_RNA"])]
    donors = sorted(set(lib.donor) | set(sections.donor) | set(zen.donor[zen.donor.str.startswith("HK")]))
    rows = []
    for d in donors:
        l = lib[lib.donor == d]
        s = sections[sections.donor == d]
        z = zen[zen.donor == d]
        dx_sc = l.diagnosis.iloc[0] if len(l) else None
        dx_st = s.diagnosis_geo.iloc[0] if len(s) else None
        dx_z = z.diagnosis.iloc[0] if len(z) else None
        labels = {x for x in [dx_sc, dx_st, dx_z] if isinstance(x, str)}
        rows.append({
            "donor": d, "dx_GSE211785_sc": dx_sc, "dx_GSE211785_visium": dx_st, "dx_zenodo": dx_z,
            "diagnosis_concordant": len(labels) <= 1,
            "age": l.age.iloc[0] if len(l) else (z.Age.iloc[0] if len(z) else None),
            "sex": l.sex.iloc[0] if len(l) else (z.Sex.iloc[0] if len(z) else None),
            "scRNA": int((rna[rna.donor == d].tech == "SC_RNA").any()),
            "snRNA": int((rna[rna.donor == d].tech == "SN_RNA").any()),
            "n_rna_cells": int(rna[rna.donor == d].n_cells.sum()),
            "snATAC": int((l.tech == "SN_ATAC").any()),
            "visium": int(len(s) > 0), "visium_HE_image": int(s.he_image.any()) if len(s) else 0,
            "cosmx": int((z.platform == "CosMx").any()), "xenium": int((z.platform == "Xenium").any()),
            "gfr_band_zenodo": z.GFR.iloc[0] if len(z) else None,
        })
    mat = pd.DataFrame(rows)
    mat["n_sources"] = ((mat.scRNA | mat.snRNA | mat.snATAC).astype(int) + mat.visium.clip(0, 1)
                        + (mat.cosmx | mat.xenium).astype(int))
    mat.to_csv(out / "donor_modality_matrix.tsv", sep="\t", index=False)
    overlap = mat[((mat.scRNA | mat.snRNA | mat.visium) > 0) & ((mat.cosmx | mat.xenium) > 0)]
    overlap.to_csv(out / "donor_overlap.tsv", sep="\t", index=False)

    rna_don = rna.drop_duplicates("donor")
    summary = {
        "GSE211785_rna_donors_by_dx": rna_don.diagnosis.value_counts().to_dict(),
        "GSE211785_rna_libraries_by_dx_tech": rna.groupby(["diagnosis", "tech"]).size()
                                                  .to_dict(),
        "GSE211785_visium_sections_by_dx": sections.diagnosis_geo.fillna("unresolved").value_counts().to_dict(),
        "GSE211785_visium_with_HE": int(sections.he_image.sum()),
        "zenodo_samples_by_dx_platform": zen.groupby(["diagnosis", "platform"]).size().to_dict(),
        "zenodo_donors_by_dx": zen.drop_duplicates("donor").diagnosis.value_counts().to_dict(),
        "GSE195460_rna_libraries_by_dx": val.diagnosis.value_counts().to_dict(),
        "GSE104954_by_platform_dx": ercb.groupby(["platform", "diagnosis"]).size().to_dict(),
        "donors_shared_discovery_and_zenodo": overlap.donor.tolist(),
        "diagnosis_discordant_donors": mat.loc[~mat.diagnosis_concordant, "donor"].tolist(),
    }
    fix = lambda d: {("|".join(map(str, k)) if isinstance(k, tuple) else str(k)): v for k, v in d.items()}
    summary = {k: (fix(v) if isinstance(v, dict) else v) for k, v in summary.items()}
    (out / "audit_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    inputs = [raw / "GSE211785/GSE211785_scRNA-seq_snRNA-seq_snATAC-seq_metadata.txt.gz",
              raw / "GSE211785/GSE211785_series_matrix.txt.gz", raw / "GSE211785/GSE211785_ST_metadata.txt.gz",
              raw / "zenodo_dkd_spatial/Diagnosis.xlsx", raw / "GSE195460/GSE195460_series_matrix.txt.gz",
              *sorted((raw / "GSE104954").glob("*series_matrix.txt.gz"))]
    write_provenance(STAGE, inputs, sorted(out.glob("*.tsv")) + [out / "audit_summary.json"], cfg["seed"])
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
