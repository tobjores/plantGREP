#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### get peaks from scATAC data (Marand et al., 2021, Cell; https://doi.org/10.1101/2020.09.27.315499) ###
# peaks called by Kerry with "simplepeaks"
if [ ! -f "${OUTDIR}/Maize_Seedling1.bed" ]; then
  unpigz -c ${EXTRADIR}/Maize_Seedling1.bed.gz \
    > ${OUTDIR}/Maize_Seedling1.bed
fi
if [ ! -f "${OUTDIR}/Maize_Seedling2.bed" ]; then
  unpigz -c ${EXTRADIR}/Maize_Seedling2.bed.gz \
    > ${OUTDIR}/Maize_Seedling2.bed
fi


### select peaks present in both (fresh and frozen) seedling samples ###
bedops -e 50% ${OUTDIR}/Maize_Seedling1.bed ${OUTDIR}/Maize_Seedling2.bed \
  > ${TMPDIR}/Maize_Seedling1_filtered.bed

bedops -e 50% ${OUTDIR}/Maize_Seedling2.bed ${OUTDIR}/Maize_Seedling1.bed \
  > ${TMPDIR}/Maize_Seedling2_filtered.bed

bedops -m ${TMPDIR}/Maize_Seedling1_filtered.bed ${TMPDIR}/Maize_Seedling2_filtered.bed \
  | sort-bed - \
  | awk -v OFS='\t' '{print $1, $2, $3, "Zm-"NR}' \
  > ${TMPDIR}/Maize_peaks.bed


### select 170 bp regions centered on peaks ###
awk -v OFS='\t' '{sub(/chr/, "", $1); mid = int(($2 + $3) / 2); print int($1), mid - 85, mid + 85, $4}' ${TMPDIR}/Maize_peaks.bed \
  | sort-bed - \
  > ${TMPDIR}/Maize_ACRs.bed


### add fragments tiling across the tb1 Hopscotch TE ###
MAXID=$(wc -l < ${TMPDIR}/Maize_ACRs.bed)

echo -e "1\t270489911\t270494841\ttb1" \
  | awk -v OFS='\t' -v maxid=${MAXID} '{i = 1; for (start = $2; start <= $3 - 170; start += 85) {print $1, start, start + 170, "Zm-"maxid - 1 + 2 * i"("$4"-"i")"; i += 0.5}}' \
  | sort-bed - ${TMPDIR}/Maize_ACRs.bed\
  > ${OUTDIR}/Maize_ACRs.bed


### get ACR sequences ###
. subscripts/get_sequences.sh Maize
