#!/usr/bin/env bash
# Master orchestrator. --resume (default) skips stages whose PROVENANCE.json exists;
# --force reruns everything; --skip=a,b skips named stages.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
eval "$(conda shell.bash hook)"
conda activate project016
MODE=resume; SKIP=""
for a in "$@"; do case $a in --force) MODE=force;; --resume) MODE=resume;; --skip=*) SKIP=${a#--skip=};; esac; done
mkdir -p logs
run(){ local st=$1; shift
  [[ ",$SKIP," == *",$st,"* ]] && { echo "skip $st"; return; }
  if [[ $MODE == resume && -f results/$st/PROVENANCE.json ]]; then echo "done $st"; return; fi
  echo ">> $st"; "$@" 2>&1 | tee logs/$st.log; }
run 00_download        bash scripts/stages/00_download.sh
run 00b_download_kpmp  bash -c "python scripts/stages/00b_download_kpmp.py expr && python scripts/stages/00b_download_kpmp.py visium_tif && python scripts/stages/00b_download_kpmp.py wsi"
run 01_cohort_audit    python scripts/stages/01_cohort_audit.py
run 02_discovery_pseudobulk python scripts/stages/02_discovery_pseudobulk.py
run 02b_prepare_gse195460   python scripts/stages/02b_prepare_gse195460.py
run 02c_library_qc     python scripts/stages/02c_library_qc.py
run 03_pseudobulk_de   python scripts/stages/03_pseudobulk_de.py
run 03b_technical_axis python scripts/stages/03b_technical_axis.py
run 03c_lodo_sn        python scripts/stages/03c_lodo_sn.py
run 03d_composition    python scripts/stages/03d_composition.py
run 04_scdisinfact     bash -c "python scripts/stages/04_scdisinfact.py primary && python scripts/stages/04_scdisinfact.py all_libraries"
run 05_lock_programme  python scripts/stages/05_lock_programme.py
run 06_validate_gse195460 python scripts/stages/06_validate_gse195460.py
run 07a_prepare_kpmp   python scripts/stages/07a_prepare_kpmp.py
run 07b_validate_kpmp  python scripts/stages/07b_validate_kpmp.py
run 07c_kpmp_concordance python scripts/stages/07c_kpmp_concordance.py
run 08_validate_ercb_arrays python scripts/stages/08_validate_ercb_arrays.py
run 09_signature_biology python scripts/stages/09_signature_biology.py
run 10a_visium_gse211785  python scripts/stages/10a_visium_gse211785.py
run 10b_visium_kpmp       python scripts/stages/10b_visium_kpmp.py
run 11a_histology_embed   bash -c "python scripts/stages/11a_histology_embed.py gse211785 && python scripts/stages/11a_histology_embed.py kpmp"
run 11c_histology_reliability python scripts/stages/11c_histology_reliability.py
run 11b_histology_association python scripts/stages/11b_histology_association.py all
run 12a_technical_robustness  python scripts/stages/12a_technical_robustness.py
run 12b_zenodo_gene_level     python scripts/stages/12b_zenodo_gene_level.py
run 12c_kpmp_clinical         python scripts/stages/12c_kpmp_clinical.py
run 13a_wsi_embed         python scripts/stages/13a_wsi_embed.py
run 13b_wsi_association   python scripts/stages/13b_wsi_association.py
run 00c_download_published_signatures bash scripts/stages/00c_download_published_signatures.sh
run 14a_concordance_permutation python scripts/stages/14a_concordance_permutation.py
run 14b_albuminuria_models   python scripts/stages/14b_albuminuria_models.py
run 14c_published_pt_states  python scripts/stages/14c_published_pt_states.py
run 14d_contamination_controls python scripts/stages/14d_contamination_controls.py
run 14e_pathways_by_cohort   python scripts/stages/14e_pathways_by_cohort.py
run 14f_site_sensitivity     python scripts/stages/14f_site_sensitivity.py
run 14g_exposure_trend       python scripts/stages/14g_exposure_trend.py
run 14h_injury_intensity_adjustment python scripts/stages/14h_injury_intensity_adjustment.py
run 15_test_summary          python scripts/stages/15_test_summary.py
