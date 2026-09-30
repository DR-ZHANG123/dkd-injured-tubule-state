#!/usr/bin/env bash
# 00c_download_published_signatures: supplementary tables of published proximal-tubule injury-state
# studies, used by stage 14c (published PT states vs the DKD-associated iPT shift).
# Writes data/raw/published_signatures/{files} and MANIFEST.tsv (url, file, bytes, md5).
# Kirita 2020 PNAS supplementary datasets are served only through the Europe PMC bundle
# (pnas.org returns 403 to scripted clients); the member files are md5-recorded after extraction
# because the bundle itself is regenerated on each request.
# Not obtainable by script: Gerhardt et al. 2023 JASN (doi:10.1681/ASN.0000000000000057) is not in the
# PMC open-access subset (Europe PMC: "not open access"; PMC/LWW pages return captcha / 403).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$ROOT/data/raw/published_signatures"
mkdir -p "$OUT"; cd "$OUT"
UA="Mozilla/5.0 (X11; Linux x86_64)"
MAN=MANIFEST.tsv
printf "url\tfile\tbytes\tmd5\n" > "$MAN"

fetch () {  # url file
  [ -s "$2" ] || curl -sSfL --http1.1 -A "$UA" --retry 8 --retry-delay 5 --retry-all-errors -o "$2" "$1"
  printf "%s\t%s\t%s\t%s\n" "$1" "$2" "$(stat -c %s "$2")" "$(md5sum "$2" | cut -d' ' -f1)" >> "$MAN"
}

# Kirita et al. 2020 PNAS 117:15874 (PMC7355049) - Dataset S2 (mouse PT subcluster DEGs)
EPMC_K="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7355049/supplementaryFiles"
if [ ! -s kirita2020_sd02.xlsx ]; then
  curl -sSfL --http1.1 --retry 8 --retry-delay 5 --retry-all-errors -o kirita2020_supp.zip "$EPMC_K"
  unzip -oj kirita2020_supp.zip pnas.2005477117.sd02.xlsx pnas.2005477117.sapp.pdf -d .
  mv pnas.2005477117.sd02.xlsx kirita2020_sd02.xlsx; mv pnas.2005477117.sapp.pdf kirita2020_SIappendix.pdf
  rm -f kirita2020_supp.zip
fi
for f in kirita2020_sd02.xlsx kirita2020_SIappendix.pdf; do
  printf "%s\t%s\t%s\t%s\n" "$EPMC_K#${f}" "$f" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -d' ' -f1)" >> "$MAN"
done
# Kirita full text (for the marker genes named in the Results text)
fetch "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7355049/fullTextXML" kirita2020_fulltext.xml

# MGI mouse-human orthology (for Kirita mouse genes)
fetch "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt" MGI_HOM_MouseHumanSequence.rpt

# Lake et al. 2023 Nature 619:585 - Supplementary Tables 1-37 (zip)
fetch "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-05769-3/MediaObjects/41586_2023_5769_MOESM6_ESM.zip" lake2023_MOESM6_supp_tables.zip
for t in 10 12 13 29; do
  f="2021-07-12390D-s6/2021-07-12390D-Supplementary Table $t.xlsx"
  [ -s "lake2023_ST$t.xlsx" ] || { unzip -oj lake2023_MOESM6_supp_tables.zip "$f" -d . && mv "2021-07-12390D-Supplementary Table $t.xlsx" "lake2023_ST$t.xlsx"; }
  printf "%s\t%s\t%s\t%s\n" "lake2023_MOESM6_supp_tables.zip#$f" "lake2023_ST$t.xlsx" "$(stat -c %s "lake2023_ST$t.xlsx")" "$(md5sum "lake2023_ST$t.xlsx" | cut -d' ' -f1)" >> "$MAN"
done

# Muto et al. 2021 Nat Commun 12:2190 - ESM file 4 (snRNA cluster DEGs) and ESM file 8 (PT_VCAM1 vs PT DEGs)
fetch "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-22368-w/MediaObjects/41467_2021_22368_MOESM4_ESM.xlsx" muto2021_MOESM4.xlsx
fetch "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-22368-w/MediaObjects/41467_2021_22368_MOESM8_ESM.xlsx" muto2021_MOESM8.xlsx
fetch "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8044133/fullTextXML" muto2021_fulltext.xml

# Abedini et al. 2024 Nat Genet 56:1712 - Supplementary Tables 1-26 (single xlsx, ~365 MB)
fetch "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-024-01802-x/MediaObjects/41588_2024_1802_MOESM3_ESM.xlsx" abedini2024_MOESM3_supp_tables.xlsx
# Abedini author manuscript full text (PMC11592391, NIHMS2034232) for marker genes named in the Results text
fetch "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=11592391&rettype=xml" abedini2024_fulltext.xml

cat "$MAN"
