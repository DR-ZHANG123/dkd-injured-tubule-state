"""Build a participant-level open-data manifest from the harvested KPMP repository index,
the OpenAccessClinicalData.csv and the TIV descriptor-score workbook."""
import pandas as pd, sys, os
D = sys.argv[1] if len(sys.argv) > 1 else "."
m = pd.read_csv(os.path.join(D, "kpmp_repository_all_files.tsv"), sep="\t", low_memory=False)
clin = pd.read_csv(os.path.join(D, "OpenAccessClinicalData_20260803.csv"))
desc = pd.read_excel(os.path.join(D, "KPMP_TIV_Descriptor_Scores_02.19.2026.xlsx"),
                     sheet_name="Data Table", header=None).iloc[2:, :3]
desc.columns = ["pid", "release", "consensus"]
desc_ids = set(desc.pid.astype(str).str.strip())

single = m[~m.redcap_id.fillna("").str.contains(r"\|") & m.redcap_id.notna()].copy()
op = single[single.access == "open"]
def count(mask, name):
    return op[mask].groupby("redcap_id").size().rename(name)
es, wf = op.experimental_strategy.fillna(""), op.workflow_type.fillna("")
lm = es.str.contains("Light Microscopic")
cols = [count(lm, "wsi_total")]
for st, nm in [("H&E stain", "wsi_HE"), ("Frozen H&E stain", "wsi_frozenHE"), ("PAS stain", "wsi_PAS"),
               ("TRI stain", "wsi_TRI"), ("SIL stain", "wsi_SIL"), ("TOL stain", "wsi_TOL"),
               ("Other stain", "wsi_other")]:
    cols.append(count(lm & (wf == st), nm))
cols.append(op[lm].groupby("redcap_id").file_size.sum().div(1e9).round(2).rename("wsi_total_GB"))
cols.append(count(es.str.contains("Segmentation Masks"), "seg_masks_files"))
for s, nm in [("Single-nucleus RNA-Seq", "snRNA_expr_matrix"), ("Single-cell RNA-Seq", "scRNA_expr_matrix"),
              ("10x Multiome", "multiome_expr_matrix")]:
    cols.append(count((es == s) & wf.str.contains("Expression Matrix"), nm))
cols.append(count((es == "Spatial Transcriptomics") & (wf == "Expression Matrix"), "visium_expr_matrix"))
cols.append(count((es == "Spatial Transcriptomics") & (wf == "Spatial Imaging Data"), "visium_spatial_images"))
for s, nm in [("CODEX", "codex_files"), ("3D Tissue Imaging and Cytometry", "3D_tissue_cytometry_files"),
              ("Imaging Mass Cytometry", "IMC_files"), ("Regional Transcriptomics", "regional_transcriptomics_files"),
              ("Spatial Metabolomics", "spatial_metabolomics_files"), ("Spatial Lipidomics", "spatial_lipidomics_files")]:
    cols.append(count(es == s, nm))
ctrl = single[single.access == "controlled"].groupby("redcap_id").size().rename("controlled_files")
tab = pd.concat(cols + [ctrl], axis=1).fillna(0)
c = clin.rename(columns={"Participant ID": "redcap_id"}).set_index("redcap_id")
keep = ["Enrollment Category", "Primary Adjudicated Category", "Diabetes History", "Diabetes Duration (Years)",
        "Hypertension History", "Baseline eGFR (ml/min/1.73m2) (Binned)", "Proteinuria (mg) (Binned)",
        "Albuminuria (mg) (Binned)", "A1c (%) (Binned)", "KDIGO Stage", "Age (Years) (Binned)", "Sex", "Race",
        "On RAAS Blockade", "Sample Type", "Tissue Source", "Protocol"]
out = c[keep].join(tab, how="outer")
out.insert(0, "in_clinical_csv", out.index.isin(c.index))
out.insert(3, "tiv_descriptor_scores", out.index.isin(desc_ids))
out = out.fillna({k: 0 for k in tab.columns})
out = out.astype({k: int for k in tab.columns if k != "wsi_total_GB"})
out.index.name = "participant_id"
out.to_csv(os.path.join(D, "kpmp_open_file_manifest.tsv"), sep="\t")
print(out.shape)
