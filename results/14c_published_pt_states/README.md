# 14c_published_pt_states — is the DKD-associated iPT shift the known failed-repair / adaptive PT programme?

Scripts: `scripts/stages/00c_download_published_signatures.sh` (download; URLs + md5 in `download_manifest.tsv`),
`scripts/lib/published_signatures.py` (parsers), `scripts/stages/14c_published_pt_states.py` (analysis).
Config key: `published_pt_states` in `config/run.yaml`. Exploratory analysis; the locked gene sets are unchanged.

## 1. Published signatures (`published_signatures.tsv`: signature, gene, rank, source, locator)

| signature | genes | source / locator |
|---|---|---|
| Kirita2020_FRPTC_DS2 | 200 | Kirita 2020 PNAS, Dataset S2 sheet MousePT_DEG, cluster 8 (the only cluster with the FR-PTC genes named in the text: Vcam1, Dcdc2a, Sema5a); avg_logFC>0, p_adj<0.05; MGI 1:1 human orthologs; top 200 by fold change |
| Kirita2020_FRPTC_text | 3 | Kirita 2020 Results text quote (Vcam1, Dcdc2a, Sema5a → VCAM1, DCDC2, SEMA5A) |
| Lake2023_aPT_snCv3 / _scCv3 | 398 / 444 | Lake 2023 Nature, Supplementary Table 13, "Full Conserved Adaptive State Gene Set", aPT rows |
| Lake2023_aPT_clinical | 26 | Lake 2023 Supplementary Table 29, column aPT |
| Lake2023_degenerative_clinical | 196 | Lake 2023 Supplementary Table 29, column "Degenerative" (pan-epithelial degenerative state) |
| Muto2021_PTVCAM1_vs_PT | 200 | Muto 2021 Nat Commun, Supplementary Data 5 (PT_VCAM1 vs PT); up, p_adj<0.05 (243 genes); top 200 |
| Muto2021_PTVCAM1_text | 5 | Muto 2021 Results text: VCAM1, HAVCR1, VIM, CD24, CD133 (=PROM1) |
| Abedini2024_iPT_ST8 | 200 | Abedini 2024 Nat Genet, Supplementary Table 8, cluster iPT; up, p_adj<0.05 (1,596 genes); top 200 |
| Abedini2024_iPT_VCAM1_text / _HAVCR1_text | 4 / 6 | Abedini 2024 Results text (iPT_VCAM1+: VCAM1, ACSL1, ASS1, ASPA; iPT_HAVCR1+: HAVCR1, NFKBIZ, IL18, ITGA3, PDGFB, ITGB1) |

The script only accepts a text-derived list if the quoted sentence appears word for word in the downloaded full text.
**Not obtained:** Gerhardt et al. 2023 JASN (doi:10.1681/ASN.0000000000000057). It is outside the PMC open-access subset, and both the PMC and publisher pages refuse scripted access. The only text reachable is the abstract, which names no genes, so no genes were included. That study is also mouse-only. Muto et al. is the 2021 Nat Commun paper (doi:10.1038/s41467-021-22368-w), which is the source of PT_VCAM1.
**Circularity caveats:** Abedini 2024 profiled the same donor cohort as the discovery set (GSE211785), and its iPT markers are iPT-vs-all-cells markers from that atlas. The Lake 2023 aPT sets were derived partly from KPMP biopsies, so they overlap the validation cohort.

## 2. Overlap (`overlap.tsv`; hypergeometric; background = 14,750 genes with a discovery SN-stratum iPT dE test)

- **12 locked markers:** 3 of 12 appear in any published signature: CD44 and RIMS2 (Kirita FR-PTC DS2; 2 of 12, p = 0.007) and SULT1C4 (Muto PT_VCAM1-vs-PT and Abedini iPT). For the union, 3 of 12 gives p = 0.047. None of the 12 are in any Lake aPT or degenerative set, and none are in the text-named markers (VCAM1, HAVCR1, etc.).
- **322 discovery dE genes (FDR < 0.10):** 14 fall in the union of 1,222 published genes, against 22.5 expected. That is mild depletion (p_depletion = 0.032) rather than enrichment. No single signature is enriched (all p ≥ 0.10). Lake aPT scCv3 has 0 of 439 against 8.1 expected (p_depletion = 2.5e-4).
- **Direction in discovery:** published genes lean slightly positive in the DKD − HKD iPT contrast. The union's median Wald statistic is 0.25 against 0.05 for other genes (MWU p = 4e-4). This comes mainly from Abedini iPT ST8 (p = 2e-5; same cohort) and Lake aPT snCv3 (p = 0.009). The published programmes are therefore at most weakly and diffusely shifted in DKD, not concentrated in the dE genes.

## 3. KPMP iPT participant pseudobulk (DKD n = 40 vs HKD n = 16 in sn; DKD n = 68 vs HKD n = 20 in sc)

`published_scores_DKD_vs_HKD.tsv`: scores are the mean z over measured genes, using the same `zmatrix` as 07c. The locked marker reproduces 07c: sn g = 0.65, Welch p = 0.067; sc g = 0.05.
- **sn:** the large failed-repair / adaptive programmes are not higher in DKD than in HKD: Kirita FR-PTC DS2 g = −0.06, Lake aPT snCv3 g = 0.08, Lake aPT scCv3 g = 0.38 (p = 0.18), Muto PT_VCAM1-vs-PT g = 0.25 (p = 0.46), Lake degenerative g = 0.21. Two sets are higher: Abedini iPT ST8 (g = 0.85, p = 0.033) and the 6-gene iPT_HAVCR1+ text set (g = 0.61, p = 0.038). The 3-gene Kirita text set is similar in size of effect but not significant (g = 0.64, p = 0.10).
- **sc:** no published signature differs (|g| ≤ 0.41, all p > 0.1). The marker does not differ in sc either.

`conditioned_models.tsv` fits marker ~ DKD + published score by OLS with HC3 errors, within DKD + HKD. Scores are in SD units. The unadjusted DKD β in sn is 0.63 (p = 0.061).
- **Conditioning on the large published programmes leaves the sn DKD effect intact:** Kirita FR-PTC DS2 retains 103% (β = 0.65, p = 0.047), Lake aPT snCv3 97% (p = 0.062), Lake aPT scCv3 87%, Muto PT_VCAM1-vs-PT 87%, and Lake degenerative 99%. With all four in one joint model, 64% is retained (β = 0.41, p = 0.16). Results are the same when the 12 markers are removed from the published sets.
- **Conditioning on the small VCAM1/HAVCR1 injury axes and on Abedini iPT markers removes much of the effect:** Abedini iPT ST8 retains 6% (marker–score ρ = 0.73), the iPT_HAVCR1+ text set 33% (ρ = 0.74), Kirita text (VCAM1/DCDC2/SEMA5A) 42%, and Lake aPT clinical 55%. So the sn DKD shift in the marker score moves together with a HAVCR1/VCAM1-type injury axis and with same-cohort iPT markers, but not with the broad FR-PTC/aPT transcriptional programmes.
- **sc:** there is no DKD effect to condition on (β = 0.05), so the sc models cannot be interpreted.

## 4. Transcriptome concordance after excluding published injury genes (`concordance_excluding_published.tsv`)

This uses the discovery dE (FDR < 0.10) and KPMP DESeq2 DKD vs HKD from 07c. The sign-agreement rate barely moves when published genes are removed:
- **sn, all dE genes:** 149/189 agree (78.8%, binomial p = 2.8e-16).
- **sn, excluding all 1,267 published genes:** 135/175 agree (77.1%, p = 1.5e-13; Spearman 0.57). The 14 dE genes inside published sets agree 14/14.
- **sn, all tested genes:** 59.4% agree, and 58.9% after the exclusion.
- **sc:** 60.9% agree (p = 0.002), and 58.6% after exclusion (p = 0.017).
- Removing any single signature changes the sn agreement by at most 0.7 percentage points.

## Bottom line

The DKD − HKD iPT shift is not simply the published failed-repair / adaptive PT programme:
- The shift shares few genes with it: 3 of 12 markers and 14 of 322 dE genes, which is fewer than expected by chance.
- The large FR-PTC and aPT signatures are not higher in DKD than in HKD in KPMP sn.
- Conditioning on those signatures leaves the sn DKD marker effect essentially unchanged.
- Excluding every published injury gene leaves transcriptome concordance at about 77%.

The marker score is, however, strongly correlated with a narrow HAVCR1/VCAM1 injury axis and with the same-cohort Abedini iPT markers, both of which are also modestly higher in DKD in sn. Conditioning on them removes most of the effect. Part of the marker shift can therefore be described as greater injury-marker intensity within iPT, while the transcriptome-wide DKD pattern extends beyond published injury genes. The sn primary effect itself is borderline (p = 0.061), and nothing replicates in sc.
