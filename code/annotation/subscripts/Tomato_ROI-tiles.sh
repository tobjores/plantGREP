#!/bin/bash

### get coordinates of regions of interest ###
unpigz -c ${EXTRADIR}/Tomato_ROIs.fa.gz \
  > ${TMPDIR}/Tomato_ROIs.fa

awk -v FS='\t' -v OFS='\t' 'substr($1, 1, 1) == ">" {sub(/chr0?/, "", $3); if ($4 == "+") $5 - 1; print $3, $5 - 1, $6, toupper(substr($2, 1, 1)) substr($2, 2), ".", $4}' ${TMPDIR}/Tomato_ROIs.fa \
  | grep -v "2872981.*Solyc11g008650" \
  | sort-bed - \
  | uniq \
  > ${OUTDIR}/Tomato_ROIs.bed


### create tiles ###
awk -v OFS='\t' '
    {
      if ($6 == "+") {
        STOP = $2; TSTART = $3; i = 1
        while (TSTART > STOP) {
          print $1, TSTART - 170, TSTART, $4"-"i;
          TSTART -= 85; i += 0.5
        }
      } else {
        STOP = $3; TSTART = $2; i = 1
        while (TSTART < STOP) {
          print $1, TSTART, TSTART + 170, $4"-"i;
          TSTART += 85; i += 0.5
        }
      }
    }' ${OUTDIR}/Tomato_ROIs.bed \
  | sort-bed - \
  > ${TMPDIR}/Tomato_ROI-tiles.bed


### get sequences ###
bedtools getfasta -nameOnly -fi ${ANNOTATION}/Tomato_genome.fa -bed ${TMPDIR}/Tomato_ROI-tiles.bed \
  > ${TMPDIR}/Tomato_ROI-tiles.fa


### merge names of duplicate sequences and remove sequences with bases other than ACGT ###
awk '{if (NR % 2 == 1) name = $0;
      else if ($1 ~ /^[ACGTacgt]*$/) {
        if ($1 in seqs) seqs[$1] = seqs[$1]";"substr(name, 2); else seqs[$1] = name
      } else print substr(name, 2) > (tmpdir"/Tomato_ROI-tiles_nonACGT.txt")
    } END{
      for (seq in seqs) print seqs[seq]"\n"seq
    }' ${TMPDIR}/Tomato_ROI-tiles.fa \
  > ${OUTDIR}/Tomato_ROI-tiles_final.fa


### remove tiles with bases other than ACGT from the final bed file ###
if [ -f "${TMPDIR}/Tomato_ROI-tiles_nonACGT.txt" ]; then
  grep -v -f ${TMPDIR}/Tomato_ROI-tiles_nonACGT.txt ${TMPDIR}/Tomato_ROI-tiles.bed > ${OUTDIR}/Tomato_ROI-tiles_final.bed
else
  cp ${TMPDIR}/Tomato_ROI-tiles.bed ${OUTDIR}/Tomato_ROI-tiles_final.bed
fi
