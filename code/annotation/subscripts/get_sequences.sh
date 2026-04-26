#!/bin/bash

### define variables ###
SPECIES=${1}


### select regions not overlapping with core promoters ###
bedops -n 1 ${OUTDIR}/${SPECIES}_ACRs.bed ${ANNOTATION}/${SPECIES}_corePro.bed \
  > ${TMPDIR}/${SPECIES}_ACRs_nonPro.bed


### shift core promoter-overlapping ACRs upstream of the core promoter ###
bedops -e 1 ${OUTDIR}/${SPECIES}_ACRs.bed ${ANNOTATION}/${SPECIES}_corePro.bed \
  | bedmap --echo --echo-map --delim '\t' - ${ANNOTATION}/${SPECIES}_corePro.bed \
  | grep -v ";" \
  | awk -v species=${SPECIES} -v outdir=${OUTDIR} 'BEGIN{file = sprintf("%s/%s_ACRs_shifted.txt", outdir, species)}
    {
      print $4 > file
      if ($10 == "+") {if ($2 < $6 - 25) {$2 = $6 - 170; $3 = $6} else {$2 = $7; $3 = $7 + 170}}
      else {if ($3 > $7 + 25) {$2 = $7; $3 = $7 + 170} else {$2 = $6 - 170; $3 = $6}}
      print $1, $2, $3, $4
    }' \
  | bedops -n 1 - ${ANNOTATION}/${SPECIES}_corePro.bed \
  | bedops -n 50% - ${TMPDIR}/${SPECIES}_ACRs_nonPro.bed \
  | sort-bed - ${TMPDIR}/${SPECIES}_ACRs_nonPro.bed \
  > ${TMPDIR}/${SPECIES}_ACRs_fixed.bed


### shift ACRs overlapping the genome borders ###
awk -v OFS='\t' -v species=${SPECIES} -v annotation=${ANNOTATION} -v outdir=${OUTDIR} '
  BEGIN{
    file = sprintf("%s/%s_genome.fa.fai", annotation, species);
    while (getline < file) {
      split($0, line, "\t");
      genome[line[1]] = line[2]
    }
    file = sprintf("%s/%s_ACRs_shifted.txt", outdir, species)
  } {
    if ($2 < 0) {print $4 >> file; $3 = $3 - $2; $2 = $2 - $2; print}
    else if ($3 > genome[$1]) {print $4 >> file; d = $3 - genome[$1]; $2 = $2 - d; $3 = $3 - d; print}
    else {print}
  }' ${TMPDIR}/${SPECIES}_ACRs_fixed.bed \
  > ${TMPDIR}/${SPECIES}_ACRs_fixed2.bed


### get IDs of identical or similar (at least 80% overlap) ACRs ###
bedmap --echo --echo-map-id --fraction-either 1 --delim '\t' ${TMPDIR}/${SPECIES}_ACRs_fixed2.bed \
  | awk -v OFS='\t' -v species=${SPECIES} -v outdir=${TMPDIR} 'BEGIN{file = sprintf("%s/%s_ACRs_identical.txt", outdir, species)} {if ($5 ~ /;/) print $5 > file; else print}' \
  | bedmap --echo --echo-map-id --fraction-either 0.801 --delim '\t' - \
  | awk -v OFS='\t' '$6 ~ /;/ {print $6}' \
  | uniq \
  > ${TMPDIR}/${SPECIES}_ACRs_similar.txt

if [ -f "${TMPDIR}/${SPECIES}_ACRs_identical.txt" ]; then
  uniq ${TMPDIR}/${SPECIES}_ACRs_identical.txt > ${TMPDIR}/tmp
  mv ${TMPDIR}/tmp ${TMPDIR}/${SPECIES}_ACRs_identical.txt
else
  touch ${TMPDIR}/${SPECIES}_ACRs_identical.txt
fi


### select ACRs to remove (prefer the one that was shifted previously) ###
awk -v species=${SPECIES} -v outdir=${TMPDIR} '{shifted[$1]++} END{
    file = sprintf("%s/%s_ACRs_identical.txt", outdir, species);
    while (getline < file) {
      split($0, line, ";");
      if (line[1] in shifted) print line[2]; else print line[1]
    }
    file = sprintf("%s/%s_ACRs_similar.txt", outdir, species);
    while (getline < file) {
      split($0, line, ";");
      if (line[1] in shifted) print line[2]; else print line[1]
    }
  }' ${OUTDIR}/${SPECIES}_ACRs_shifted.txt \
  > ${TMPDIR}/${SPECIES}_ACRs_remove.txt


### remove ACRs ###
awk -v species=${SPECIES} -v outdir=${TMPDIR} 'BEGIN{
    file = sprintf("%s/%s_ACRs_remove.txt", outdir, species);
    while (getline < file) {
      remove[$0]++
    }
  } !($4 in remove) {print}' ${TMPDIR}/${SPECIES}_ACRs_fixed2.bed \
  > ${OUTDIR}/${SPECIES}_ACRs_final.bed


### get sequences (remove sequences with bases other than ACGT) ###
bedtools getfasta -nameOnly -fi ${ANNOTATION}/${SPECIES}_genome.fa -bed ${OUTDIR}/${SPECIES}_ACRs_final.bed \
  | awk -v species=${SPECIES} -v tmpdir=${TMPDIR} '{if (NR % 2 == 1) name = $0; else if ($1 ~ /^[ACGTacgt]*$/) print name"\n"$1; else print substr(name, 2) > (tmpdir"/"species"_nonACGT.txt")}' \
  > ${OUTDIR}/${SPECIES}_ACRs_final.fa


### remove ACRs with bases other than ACGT from the final bed file ###
if [ -f "${TMPDIR}/${SPECIES}_nonACGT.txt" ]; then
  grep -v -f ${TMPDIR}/${SPECIES}_nonACGT.txt ${OUTDIR}/${SPECIES}_ACRs_final.bed > ${OUTDIR}/tmp
  mv ${OUTDIR}/tmp ${OUTDIR}/${SPECIES}_ACRs_final.bed
fi