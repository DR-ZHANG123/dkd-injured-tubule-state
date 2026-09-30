#!/usr/bin/env bash
# 00_download: fetch all public inputs into data/raw and write MD5 manifest.
# Idempotent: wget -c resumes / skips completed files.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
R=data/raw
GEO=https://ftp.ncbi.nlm.nih.gov/geo/series
get(){ mkdir -p "$(dirname "$2")"; wget -q -c --tries=20 --waitretry=30 --retry-connrefused -O "$2" "$1" || { rm -f "$2"; echo "FAIL $1"; return 1; }; echo "ok $2"; }
# NCBI throttles >3 parallel connections; keep concurrency low
export -f get

# --- GSE211785 (discovery atlas: sc/snRNA + Visium with H&E) ---
B=$GEO/GSE211nnn/GSE211785/suppl
for f in GSE211785_scRNA-seq_snRNA-seq_snATAC-seq_metadata.txt.gz GSE211785_ST_metadata.txt.gz \
         GSE211785_EXPORT_ST_umap.txt.gz GSE211785_EXPORT_scRNA-seq_snRNA-seq_snATAC-seq_umap.txt.gz \
         GSE211785_EXPORT_ST_counts.rds.gz GSE211785_RAW.tar filelist.txt \
         GSE211785_V11Y24-076-A1-2.json.gz GSE211785_V11Y24-076-B1.json.gz GSE211785_V11Y24-076-C1.json.gz \
         GSE211785_V11Y24-076-D1.json.gz GSE211785_V11Y24-076_A1-HK3035.tif.gz GSE211785_V11Y24-076_B1-HK2852.tif.gz \
         GSE211785_V11Y24-076_C1-HK2871.tif.gz GSE211785_V11Y24-076_D1-HK2873.tif.gz GSE211785_V11Y24-079-A1_1.json.gz \
         GSE211785_V11Y24-079-B1_1.json.gz GSE211785_V11Y24-079-C1-1.json.gz GSE211785_V11Y24-79_A1_HK2770.tif.gz \
         GSE211785_V11Y24-79_B1_HK2844.tif.gz GSE211785_Susztak_SC_SN_ATAC_merged_PreSCVI_final.h5ad.gz \
         GSE211785_Susztak_KPMP_SC_SN_ATAC_merged_PostSCVI_final.h5ad.gz; do
  echo "$B/$f $R/GSE211785/$f"
done | xargs -P 3 -n 2 bash -c 'get "$0" "$1"' 

# --- GSE195460 (independent snRNA validation; RNA matrices only) ---
B=$GEO/GSE195nnn/GSE195460/suppl
for s in Control1 Control2 Control3 Control4 Control5 DN1 DN2 DN3; do
  get $B/GSE195460_${s}_filtered_feature_bc_matrix.h5 $R/GSE195460/GSE195460_${s}_filtered_feature_bc_matrix.h5 & done
# samples whose RNA matrix is stored only at GSM level
A=https://ftp.ncbi.nlm.nih.gov/geo/samples
get $A/GSM5837nnn/GSM5837792/suppl/GSM5837792_Control6_filtered_feature_bc_matrix.h5 $R/GSE195460/GSM5837792_Control6_filtered_feature_bc_matrix.h5 &
get $A/GSM5837nnn/GSM5837797/suppl/GSM5837797_DN4_filtered_feature_bc_matrix.h5 $R/GSE195460/GSM5837797_DN4_filtered_feature_bc_matrix.h5 &
get $A/GSM5837nnn/GSM5837799/suppl/GSM5837799_DN5_filtered_feature_bc_matrix.h5 $R/GSE195460/GSM5837799_DN5_filtered_feature_bc_matrix.h5 &
get $GEO/GSE195nnn/GSE195460/matrix/GSE195460_series_matrix.txt.gz $R/GSE195460/GSE195460_series_matrix.txt.gz &
wait

# --- GSE104954 (ERCB tubulointerstitial microarray, cross-disease) ---
B=$GEO/GSE104nnn/GSE104954
get $B/matrix/GSE104954-GPL22945_series_matrix.txt.gz $R/GSE104954/GSE104954-GPL22945_series_matrix.txt.gz
get $B/matrix/GSE104954-GPL24120_series_matrix.txt.gz $R/GSE104954/GSE104954-GPL24120_series_matrix.txt.gz

# --- Zenodo 19868428 (2026 DKD spatial atlas, CosMx + Xenium) ---
Z=https://zenodo.org/records/19868428/files
get "$Z/Diagnosis.xlsx?download=1" $R/zenodo_dkd_spatial/Diagnosis.xlsx
get "$Z/spatial_adata_xenium_cosmx_zenodo.h5ad?download=1" $R/zenodo_dkd_spatial/spatial_adata_xenium_cosmx_zenodo.h5ad

# --- HGNC gene table (sex-chromosome and approved-symbol filters) ---
get https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt $R/reference/hgnc_complete_set.txt
M=https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs
for f in h.all.v2024.1.Hs.symbols.gmt c2.cp.reactome.v2024.1.Hs.symbols.gmt c5.go.bp.v2024.1.Hs.symbols.gmt; do
  get $M/$f $R/reference/$f; done

( cd $R && find . -type f ! -name MD5SUMS -print0 | sort -z | xargs -0 md5sum > MD5SUMS )
echo "00_download done"
