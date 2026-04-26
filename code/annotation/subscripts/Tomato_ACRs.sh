#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### get peaks from ATAC data (Hendelman et al., 2021, Cell; https://doi.org/10.1016/j.cell.2021.02.001) ###
# peaks called by Sayeh with "simplepeaks"
if [ ! -f "${OUTDIR}/Tomato_peaks_1.bed" ]; then
  unpigz -c ${EXTRADIR}/Tomato_peaks_1.bed.gz \
    > ${OUTDIR}/Tomato_peaks_1.bed
fi
if [ ! -f "${OUTDIR}/Tomato_peaks_2.bed" ]; then
  unpigz -c ${EXTRADIR}/Tomato_peaks_2.bed.gz \
    > ${OUTDIR}/Tomato_peaks_2.bed
fi
if [ ! -f "${OUTDIR}/Tomato_peaks_3.bed" ]; then
  unpigz -c ${EXTRADIR}/Tomato_peaks_3.bed.gz \
    > ${OUTDIR}/Tomato_peaks_3.bed
fi
if [ ! -f "${OUTDIR}/Tomato_peaks_4.bed" ]; then
  unpigz -c ${EXTRADIR}/Tomato_peaks_4.bed.gz \
    > ${OUTDIR}/Tomato_peaks_4.bed
fi
if [ ! -f "${OUTDIR}/Tomato_peaks_5.bed" ]; then
  unpigz -c ${EXTRADIR}/Tomato_peaks_5.bed.gz \
    > ${OUTDIR}/Tomato_peaks_5.bed
fi
if [ ! -f "${OUTDIR}/Tomato_peaks_6.bed" ]; then
  unpigz -c ${EXTRADIR}/Tomato_peaks_6.bed.gz \
    > ${OUTDIR}/Tomato_peaks_6.bed
fi


### create union peaks ###
bedops -m ${OUTDIR}/Tomato_peaks_1.bed ${OUTDIR}/Tomato_peaks_2.bed ${OUTDIR}/Tomato_peaks_3.bed ${OUTDIR}/Tomato_peaks_4.bed ${OUTDIR}/Tomato_peaks_5.bed ${OUTDIR}/Tomato_peaks_6.bed\
  > ${TMPDIR}/Tomato_peaks_union.bed


### select union peaks detected in at least four samples ###
bedmap --echo --indicator --delim '\t' ${TMPDIR}/Tomato_peaks_union.bed ${OUTDIR}/Tomato_peaks_1.bed \
  | bedmap --echo --indicator --delim '\t' - ${OUTDIR}/Tomato_peaks_2.bed \
  | bedmap --echo --indicator --delim '\t' - ${OUTDIR}/Tomato_peaks_3.bed \
  | bedmap --echo --indicator --delim '\t' - ${OUTDIR}/Tomato_peaks_4.bed \
  | bedmap --echo --indicator --delim '\t' - ${OUTDIR}/Tomato_peaks_5.bed \
  | bedmap --echo --indicator --delim '\t' - ${OUTDIR}/Tomato_peaks_6.bed \
  | awk -v OFS='\t' 'BEGIN{i = 1} ($4 + $5 + $6 + $7 + $8 + $9) >= 4 {print $1, $2, $3, "Sl-"i; i++}' \
  > ${TMPDIR}/Tomato_peaks_union_4plus.bed


### select 170 bp regions centered on peaks ###
awk -v OFS='\t' '{sub(/chr/, "", $1); mid = int(($2 + $3) / 2); print int($1), mid - 85, mid + 85, $4}' ${TMPDIR}/Tomato_peaks_union_4plus.bed \
  | sort-bed - \
  > ${OUTDIR}/Tomato_ACRs.bed


### get ACR sequences ###
. subscripts/get_sequences.sh Tomato


### get tiles across promoter regions of genes of interest ###
. subscripts/Tomato_ROI-tiles.sh