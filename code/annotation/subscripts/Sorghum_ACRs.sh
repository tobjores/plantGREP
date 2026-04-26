#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### get peaks from ATAC data (Parvathaneni et al., 2021, bioRxiv; https://doi.org/10.1101/2020.08.07.240580) ###
# peaks called by Kerry with MACS2
if [ ! -f "${OUTDIR}/Sorghum_peaks_1.bed" ]; then
  unpigz -c ${EXTRADIR}/Sorghum_peaks_1.bed.gz \
    > ${OUTDIR}/Sorghum_peaks_1.bed
fi
if [ ! -f "${OUTDIR}/Sorghum_peaks_2.bed" ]; then
  unpigz -c ${EXTRADIR}/Sorghum_peaks_2.bed.gz \
    > ${OUTDIR}/Sorghum_peaks_2.bed
fi
if [ ! -f "${OUTDIR}/Sorghum_peaks_3.bed" ]; then
  unpigz -c ${EXTRADIR}/Sorghum_peaks_3.bed.gz \
    > ${OUTDIR}/Sorghum_peaks_3.bed
fi


### create union peaks ###
bedops -m ${OUTDIR}/Sorghum_peaks_1.bed ${OUTDIR}/Sorghum_peaks_2.bed ${OUTDIR}/Sorghum_peaks_3.bed \
  > ${TMPDIR}/Sorghum_peaks_union.bed


### select union peaks detected in at least two samples ###
bedmap --echo --indicator --delim '\t' ${TMPDIR}/Sorghum_peaks_union.bed ${OUTDIR}/Sorghum_peaks_1.bed \
  | bedmap --echo --indicator --delim '\t' - ${OUTDIR}/Sorghum_peaks_2.bed \
  | bedmap --echo --indicator --delim '\t' - ${OUTDIR}/Sorghum_peaks_3.bed \
  | awk 'BEGIN{i = 1} ($4 + $5 + $6) >= 2 {print $1, $2, $3, "Sb-"i; i++}' \
  > ${TMPDIR}/Sorghum_peaks_union_2plus.bed


### select 170 bp regions centered on peaks ###
awk -v OFS='\t' '{sub(/chr/, "", $1); mid = int(($2 + $3) / 2); print int($1), mid - 85, mid + 85, $4}' ${TMPDIR}/Sorghum_peaks_union_2plus.bed \
  | sort-bed - \
  > ${OUTDIR}/Sorghum_ACRs.bed


### get ACR sequences ###
. subscripts/get_sequences.sh Sorghum
