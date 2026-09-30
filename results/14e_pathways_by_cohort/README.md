# 14e_pathways_by_cohort

Script: `scripts/stages/14e_pathways_by_cohort.py`. Config: `gsea` (min 15, max 500 genes, 2000 permutations) and `pathways_by_cohort` in `config/run.yaml`. Seed 20260929.

## What was done

The pathway annotation in 09/12a ranked genes by a Stouffer combination of the two cohorts. Here each cohort is ranked on its own. We ran preranked GSEA (the `prerank` function of stage 09) on the DKD-vs-HKD Wald statistic of iPT pseudobulk (log2FC/SE), one run per cohort:

| ranking | source | genes ranked (full / technical genes removed) |
|---|---|---|
| discovery_SN | `03_pseudobulk_de/fine_SN_only/iPT.tsv` (3 DKD vs 8 HKD SN libraries) | 17,478 / 16,951 |
| discovery_SN_ribo_adj | `12a_technical_robustness/discovery_iPT_DKD_vs_HKD_adj_frac_ribo.tsv` | 17,478 / 16,951 |
| KPMP_sn | `07c_kpmp_concordance/kpmp_iPT_DKD_vs_HKD_sn.tsv` (40 DKD vs 16 HKD) | 14,224 / 13,693 |
| KPMP_sn_ribo_adj | `12a_technical_robustness/kpmp_sn_iPT_DKD_vs_HKD_adj_frac_ribo.tsv` | 14,224 / 13,693 |

We used two gene filters:
- `excl_technical` is the primary filter. It removes ribosomal-protein, mitochondrial-ribosomal and MT- genes, plus the members of the 12a translation/rRNA exclusion terms. It is run on Hallmark, Reactome and GO BP. GO BP is included because the pattern-specification terms are GO BP terms.
- `full` keeps all genes and is run on Hallmark and Reactome only. Its purpose is to keep the translation terms testable, because the primary filter removes most of their genes.

A pathway is called **consistent** when NES has the same sign in both cohorts and FDR < 0.25 in both. This is checked for four discovery × KPMP pairings. FDR 0.000 from gseapy means below the resolution of 2000 permutations and is written as < 0.001.

## Key pathways (primary filter `excl_technical`; NES / FDR)

| term | discovery_SN | discovery_SN_ribo_adj | KPMP_sn | KPMP_sn_ribo_adj | verdict |
|---|---|---|---|---|---|
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | 1.07 / 1.00 | 1.21 / 1.00 | 2.20 / <0.001 | 1.58 / 0.087 | KPMP only; not significant in discovery |
| HALLMARK_IL6_JAK_STAT3_SIGNALING | 1.44 / 0.19 | 1.14 / 1.00 | 1.97 / <0.001 | 1.45 / 0.074 | consistent (unadjusted discovery with either KPMP version); lost if discovery is ribo-adjusted |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | 1.22 / 0.70 | 1.09 / 0.97 | 1.85 / 0.001 | 1.42 / 0.070 | KPMP only |
| REACTOME_SIGNALING_BY_BMP | −1.90 / 0.29 | −1.96 / 0.34 | −1.67 / 0.29 | −1.50 / 0.25 | down in all four; FDR 0.25–0.34 (full filter: −1.96/0.18, −1.96/0.028, −1.69/0.24, −1.48/0.21 → consistent) |
| GOBP_ANTERIOR_POSTERIOR_PATTERN_SPECIFICATION | −2.30 / 0.014 | −2.19 / 0.060 | −1.49 / 0.33 | −1.72 / 0.20 | consistent with ribo-adjusted KPMP only |
| GOBP_PATTERN_SPECIFICATION_PROCESS | −2.05 / 0.072 | −1.86 / 0.19 | −1.54 / 0.32 | −1.62 / 0.27 | same direction; KPMP does not reach FDR 0.25 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 0.93 / 1.00 | −1.12 / 0.76 | 1.68 / 0.005 | −2.37 / <0.001 | sign depends on ribosome adjustment; not consistent |
| REACTOME_EUKARYOTIC_TRANSLATION_ELONGATION (full) | 1.87 / 0.093 | −2.54 / <0.001 | 3.15 / <0.001 | −2.75 / <0.001 | reverses with ribosome adjustment in both cohorts (technical) |
| REACTOME_TRANSLATION (full) | 1.58 / 0.40 | −1.49 / 0.39 | 2.33 / <0.001 | −2.37 / <0.001 | as above |
| HALLMARK_MYC_TARGETS_V1 | 1.55 / 0.083 | 1.18 / 1.00 | 1.78 / 0.002 | −1.53 / 0.071 | consistent only without adjustment; reverses in KPMP_adj |
| HALLMARK_ALLOGRAFT_REJECTION | 1.55 / 0.11 | 1.16 / 1.00 | 2.15 / <0.001 | 1.49 / 0.096 | consistent (unadjusted discovery) |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | 0.97 / 0.99 | 0.85 / 1.00 | 1.92 / 0.001 | 0.89 / 0.84 | KPMP unadjusted only |

## Which pathways are consistent in both cohorts

- **Up in DKD:** IL6-JAK-STAT3 and allograft rejection. Both have moderate discovery FDR (0.11–0.19) and hold against ribosome-adjusted KPMP. Neither survives ribosome adjustment of the discovery ranking. Both terms are rich in immune genes, so read them together with stage 14d. MYC targets V1/V2 are consistent only without adjustment.
- **Down in DKD:** BMP signalling (full-gene filter, FDR 0.03–0.24), anterior-posterior pattern specification (discovery FDR 0.014–0.060; KPMP 0.20 after ribosome adjustment), and HALLMARK_KRAS_SIGNALING_DN (discovery 0.001, KPMP 0.12). There is also a cluster of developmental/patterning GO BP terms: proximal-distal pattern formation, neuron fate specification, spinal cord / photoreceptor differentiation, and SMO activation / NOTCH2-3. This is the most reproducible direction; it appears in every pairing (`consistent_pathways.tsv`).
- **Not consistent:** TNFα/NFκB and EMT reach FDR < 0.1 only in KPMP. In the 09/12a combined ranking they looked like cross-cohort signals; that was driven by KPMP. OXPHOS, translation and MYC_V1 change sign with ribosome-fraction adjustment in both cohorts, so they are technical.
- **Pathway-level agreement** (Spearman r of NES across all terms, `pathway_level_agreement.tsv`, primary filter, discovery_SN vs KPMP_sn / vs KPMP_sn_ribo_adj): Hallmark 0.49 / 0.31, GO BP 0.32 / 0.14, Reactome 0.16 / −0.04.

## Files

- `gsea_by_cohort.tsv`: every GSEA result (ranking × gene filter × library).
- `ranking_sizes.tsv`: genes per ranking.
- `key_pathways_by_cohort.tsv`: key terms, wide format (NES and FDR per filter × ranking).
- `key_pathways_consistency.tsv`: key terms × discovery/KPMP pairing, with the consistency call.
- `pathway_consistency_all_terms.tsv`, `consistent_pathways.tsv`: the same call for every term.
- `pathway_level_agreement.tsv`: NES correlation and counts of consistent terms per library.
- `PROVENANCE.json`
