# 14d_contamination_controls

Script: `scripts/stages/14d_contamination_controls.py`, helpers in `scripts/lib/contamination.py`. Config: `contamination_controls` in `config/run.yaml`. Seed 20260929.

**Question.** The DKD-vs-HKD iPT shift includes CD44 and the other `DKDiPT_up` genes. Could it be leukocyte or ambient RNA carried into iPT nuclei, rather than a tubular cell state?

## Data and units

The unit is always the sample: the participant for KPMP, the SN library for discovery. Cells are never units.

- **KPMP sn / sc.** Input is the QC counts (`data/processed/KPMP_{sn,sc}_qc_counts.h5ad`) with the 07b labels. Raw counts are summed per participant × compartment. Compartments are iPT; healthy PT (PT_S1-3); immune (Myeloid + Lymphoid labels of `lineage_map`); stroma (Stroma lineage).
- **Discovery.** SN libraries that pass QC (`02c_library_qc`). The compartments are built from the stage-02 fine pseudobulk `pb_fine.h5ad`, which holds the raw counts of `GSE211785_rna_counts.h5ad` summed per library × cell type.
- **Sample filter.** An iPT sample is kept when it has ≥ 20 iPT cells. This leaves 3 DKD / 7 HKD discovery libraries, 40 / 16 KPMP sn participants and 68 / 20 KPMP sc participants.
- **Scores.** Programme scores are the mean z (the same function as 07b/07c). They reproduce the 07b primary result: KPMP sn `DKDiPT_up` g = 0.646, `DKDiPT_down` g = −0.660.

## Contamination indices, all on the iPT pseudobulk

| index | definition |
|---|---|
| `frac_leukocyte_genes` | Fraction of iPT counts from PTPRC, CD74, HLA-DRA, LYZ, CD3E, CD2, MS4A1, CD52, CORO1A, LCP1, TYROBP, FCER1G, C1QA, NKG7. |
| `frac_leukocyte_strict` | Same, but without MHC-II and C1QA/FCER1G. CD74 and HLA-DRA are inducible in injured tubular cells. |
| `frac_stroma_strict` | DCN, LUM, COL1A2, COL3A1, PDGFRB, C7. |
| `frac_ig_genes` | IGKC, IGHG1, IGHA1, IGLC2, JCHAIN. Kept out of the leukocyte index. All zero in KPMP sc, so not estimable there. |
| `frac_nonpt_markers` | Fraction from genes with pooled immune+stroma CPM ≥ 4 × pooled PT CPM and ≥ 10 CPM. There are 945 such genes (KPMP sn), 1026 (sc) and 1492 (discovery). Programme genes are excluded from this list. |
| `c_immune`, `c_stroma` | Ratio estimate of the contaminating share: iPT CPM of the strict lineage genes divided by their CPM in that compartment. This is an upper bound, because it assumes iPT never expresses those genes. |
| `r_ambient_all`, `r_ambient_markers` | Pearson r between the iPT log2CPM and the same sample's own immune+stroma log2CPM, across all expressed genes or across non-PT markers. Needs ≥ 20 non-PT cells. |

## Results

**(1) Is DKD iPT more contaminated than HKD iPT?** (`index_DKD_vs_HKD.tsv`)
- **KPMP sn.** DKD iPT carries slightly more leukocyte-gene signal: strict fraction 4.8e-5 vs 3.0e-5, g = 0.46, Welch p = 0.040, MWU p = 0.20. The canonical-list fraction gives g = 0.45, p = 0.046. The estimated immune share is 6.0% vs 3.6%, an upper bound. The stromal, Ig, non-PT-marker and profile-correlation indices do not differ (|g| ≤ 0.20, p ≥ 0.42). The absolute leukocyte fraction is tiny: < 0.02% of iPT counts.
- **KPMP sc.** No excess in DKD. The leukocyte indices point slightly the other way (g −0.24 to −0.32). `r_ambient_markers` is lower in DKD (g = −0.47, p = 0.048).
- **Discovery SN.** The indices are numerically higher in DKD (g 0.2–0.97), but none is significant with 3 vs 7 libraries (p ≥ 0.33).

**(2) Do the programme scores survive adjustment?** (`score_models_adjusted.tsv`, OLS HC3, score ~ DKD + standardised index)
- **KPMP sn `DKDiPT_up`.**
  - Unadjusted: β = 0.298 (p = 0.061).
  - Adjusted for the strict leukocyte fraction: β = 0.259 (p = 0.12), about 13% attenuation.
  - Adjusted for stromal, Ig, non-PT-marker or ambient-correlation indices: β = 0.282–0.295 (p = 0.061–0.071).
  - Adjusted for all four indices jointly: β = 0.259 (p = 0.13).
  - Residualising each gene on the leukocyte fraction first (the conservative option) gives g = 0.545 (p = 0.12) against the original 0.646.
- **KPMP sn `DKDiPT_down`.**
  - Unadjusted: β = −0.365 (p = 0.036).
  - Leukocyte-adjusted: β = −0.325 (p = 0.067).
  - Stroma- or ambient-adjusted: β = −0.360 to −0.383 (p = 0.030–0.050).
  - Joint: β = −0.350 (p = 0.060).
- **Discovery SN.** The DKD coefficients are essentially unchanged by any single index: `DKDiPT_up` β = 1.58 → 1.55–1.70, and `DKDiPT_down` β = −1.72 → −1.67 to −1.92. They stay at p < 0.06 with n = 10. The joint model for `DKDiPT_up` gives p = 0.23 because of over-parameterisation with n = 10.
- **KPMP sc.** No DKD effect before or after adjustment.
- **Reference sets, KPMP sn.** `DEonly_up` keeps β = 0.34 (p = 0.038) after leukocyte adjustment. `SharedInjury_up` goes from β = 0.154 to 0.007.

**(3) Where are the marker genes expressed?** (`marker_specificity.tsv`: participant-mean CPM in DKD+HKD, count share by compartment; `expected_contamination_programme_genes.tsv`)
- **CD44 is mostly a leukocyte/stromal transcript in KPMP sn.**
  - Mean CPM: iPT 18.7, healthy PT 12.3, immune 146.7, stroma 149.1. That is 7.8× and 8.0× the iPT level.
  - 53% of CD44 counts from DKD+HKD participants come from immune nuclei; 2.2% come from iPT.
  - The pattern is the same in discovery SN (immune/iPT 12.4×, stroma/iPT 10.4×) and in KPMP sc (13.7×).
  - Contamination could supply a median 35% (DKD) and 28% (HKD) of the iPT CD44 signal in KPMP sn (upper bound). The DKD−HKD difference in iPT CD44 is small (+4.0 CPM), and up to 70% of it could be explained by the higher immune contamination in DKD.
  - As a single gene, CD44 does not replicate in KPMP sn (β = 0.09 z, p = 0.73), and it turns negative after leukocyte adjustment (−0.08).
  - In discovery, iPT CD44 is 7.5× higher in DKD (51.8 vs 6.9 CPM). Contamination could explain at most 23% of that difference. The per-gene effect keeps 87% of its size after adjustment (β 1.45 → 1.25), but loses significance at n = 10 (p 0.024 → 0.28).
- **Other genes enriched in immune or stromal cells relative to iPT (KPMP sn):**
  - SYTL3: immune 2.1×.
  - RIMBP2: immune 3.1×.
  - DLC1: stroma 7.0×.
  - NDRG1: about 1.3× in both.
  - EPB41L4B, CASZ1, FAT4, PDZRN3: 1.5–2.8×.

  None of these replicates as a single gene in KPMP sn, except FAT4 (down, p = 0.02, which survives every adjustment).
- **The KPMP sn single-gene replications are iPT-enriched and not contamination.**
  - RIMS2: β = 0.94, p = 0.013; iPT CPM 141 vs immune 33 and stroma 28.
  - SYT14: β = 0.68, p = 0.037; ratio ≤ 0.30.
  - SULT1C4: β = 0.51, p = 0.12; ratio ≤ 0.19.
  - Their expected contamination share is ≤ 2%, and their adjusted β is unchanged.

**(4) Does the transcriptome concordance hold without non-PT-enriched genes?** (`concordance_specificity_filtered.tsv`; discovery SN iPT dE FDR < 0.10 vs KPMP sn DKD-vs-HKD)
- **All genes.** 189 genes, 78.8% same sign, binomial p < 1e-4, Spearman 0.58.
- **Filtered on KPMP sn expression.**
  - Genes with immune or stromal CPM ≤ 1 × iPT CPM: 73 genes, 82.2% agree (p < 1e-4).
  - Stricter cut at ≤ 0.5 ×: 18 genes, 94.4% agree (p = 1e-4).
  - Only the genes higher in non-PT cells: 116 genes, 76.7% agree.
- **Filtered on discovery SN expression instead:** 87 genes, 83.9% agree; strict cut, 28 genes, 92.9% agree.
- **By gene class.**
  - Protein-coding only: 157 genes, 77.7% agree.
  - Non-coding RNA only: 30 genes, 90.0% agree.
  - Non-synaptic only: 171 genes, 77.2% agree.
  - Protein-coding, non-synaptic and not non-PT-higher: 51 genes, 78.4% agree.
- **Conclusion.** The concordance does not depend on genes that could come from leukocytes or stroma. It is, if anything, higher among iPT-enriched genes.
- **Caveat.** In KPMP sn, 70% of all tested genes have a higher participant-mean CPM in immune or stromal nuclei than in iPT. So the ratio ≤ 1 filter is permissive, and the strict filter is the more informative one.

**Gene class** (`gene_class_composition.tsv`, HGNC `locus_group`; synaptic = member of any GO BP term matching "SYNAP")
- Discovery FDR < 0.10 genes compared with all tested genes:
  - Protein-coding: 50% vs 81%.
  - HGNC non-coding RNA: 13% vs 6%.
  - Symbol not in HGNC: 37% vs 13%. Mostly clone-based lncRNA / unannotated loci; these were removed a priori from the locked programmes.
  - Synaptic: 5.6% vs 5.1%, so not enriched.
- In the locked programme, 4 of 20 genes are HGNC non-coding RNAs (LINC01090, ZFPM2-AS1 and PAX8-AS1 in `up`; LINC02055 in `down`). ZFPM2-AS1 and PAX8-AS1 are not in the KPMP matrices. Four genes carry the synaptic GO flag (RIMS2, NRCAM, RIMBP2, PDZRN3).
- **Gene length was not assessed:** there is no verified local gene-length annotation (no GTF under `data/raw`).

## Interpretation

- **Leukocyte or ambient RNA does not explain the programme-level DKD-vs-HKD shift.**
  - In KPMP sn, contamination indices are similar between groups, apart from a small leukocyte excess in DKD.
  - Adjusting for any index changes the DKD coefficient by ≤ 13%.
  - Transcriptome concordance is at least as high among iPT-enriched genes.
- **CD44 is the exception.** In KPMP sn it is predominantly expressed by immune and stromal nuclei, and its small iPT DKD excess can largely be accounted for by contamination (upper-bound estimate). CD44 should not be presented as an iPT-intrinsic marker on KPMP evidence. The KPMP-replicating genes are iPT-enriched: RIMS2, SYT14 and SULT1C4.
- **Limits.**
  - The ratio estimate assumes zero iPT expression of the strict lineage genes, so the shares are upper bounds.
  - Discovery has only 3 DKD libraries.
  - KPMP sc shows no programme effect to adjust.

## Files

- `contamination_indices_by_sample.tsv`: indices and programme scores per sample.
- `index_DKD_vs_HKD.tsv`: DKD vs HKD tests on each index.
- `score_models_adjusted.tsv`: programme scores adjusted for each index.
- `per_gene_models_adjusted.tsv`: the same models per programme gene.
- `marker_specificity.tsv`: programme genes and index genes, by cohort.
- `all_gene_specificity_KPMP_sn.tsv.gz`: specificity for every gene in KPMP sn.
- `expected_contamination_programme_genes.tsv`: expected versus observed iPT CPM for programme genes.
- `concordance_specificity_filtered.tsv`: concordance under each specificity filter.
- `gene_class_composition.tsv`: gene-class fractions.
- `PROVENANCE.json`
