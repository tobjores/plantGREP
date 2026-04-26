#!/bin/bash

### define variables ###
SPECIES_A="${1}"
SPECIES_B="${2}"


### get CNSs in the second species that contain at least half an ACR ###
bedops --element-of 50% ${OUTDIR}/${SPECIES_B}_ACRs_final.bed ${OUTDIR}/${SPECIES_B}_CNSs.bed \
  | bedmap --echo --echo-map-id - ${OUTDIR}/${SPECIES_B}_CNSs.bed \
  | awk '{print $4}' \
  | sed 's/|/\t/' \
  > ${TMPDIR}/${SPECIES_B}_ACR2CNS.tsv


### get regions (CNS + 500bp on each side) in species 1 corresponding to CNSs in species 2 that contain at least half an ACR ###
awk -v OFS='\t' -v species="${SPECIES_A}" -v outdir=${OUTDIR} '{
    split($2, tmp, ";")
    for (i in tmp) {CNSs[tmp[i]]++}
  } END{
    file = sprintf("%s/%s_CNSs.bed", outdir, species)
    while (getline < file) {
      split($0, line, "\t")
      if (line[4] in CNSs) print line[1], line[2] - 500, line[3] + 500, line[4]
    }
  }' ${TMPDIR}/${SPECIES_B}_ACR2CNS.tsv \
  | sort-bed - \
  > ${TMPDIR}/${SPECIES_A}_regions.bed


### shift regions overlapping the genome borders ###
awk -v OFS='\t' -v species="${SPECIES_A}" -v annotation=${ANNOTATION} '
  BEGIN{
    file = sprintf("%s/%s_genome.fa.fai", annotation, species)
    while (getline < file) {
      split($0, line, "\t")
      genome[line[1]] = line[2]
    }
  } {
    if ($2 < 0) {$3 = $3 - $2; $2 = $2 - $2; print}
    else if ($3 > genome[$1]) {d = $3 - genome[$1]; $2 = $2 - d; $3 = $3 - d; print}
    else {print}
  }' ${TMPDIR}/${SPECIES_A}_regions.bed \
  | sort-bed - \
  > ${TMPDIR}/${SPECIES_A}_regions_fixed.bed


### get sequences ###
bedtools getfasta -name -fi ${ANNOTATION}/${SPECIES_A}_genome.fa -bed ${TMPDIR}/${SPECIES_A}_regions_fixed.bed \
  > ${TMPDIR}/${SPECIES_A}_regions.fa


### blast species 2 ACRs against species 1 regions ###
makeblastdb -in ${TMPDIR}/${SPECIES_A}_regions.fa -dbtype nucl -parse_seqids
blastn -db ${TMPDIR}/${SPECIES_A}_regions.fa -query ${OUTDIR}/${SPECIES_B}_ACRs_final.fa -outfmt 6 -out ${TMPDIR}/${SPECIES_B}_${SPECIES_A}.blast


### get coordinates for alignments longer than 100 bases with 80% to 98.5% identity ###
awk -v OFS='\t' -v speciesA="${SPECIES_A}" -v speciesB="${SPECIES_B}" -v tmpdir=${TMPDIR} 'BEGIN{
    file = sprintf("%s/%s_ACR2CNS.tsv", tmpdir, speciesB)
    while (getline < file) {
      split($0, ids, "\t")
      acr2cns[ids[1]] = ids[2]
    }
  } $3 < 98.5 && $3 >= 80 && $4 > 100 {
    SpACNS = substr($2, 1, index($2, ":") - 1)
    SpBCNS = acr2cns[$1]
    if (SpBCNS ~ SpACNS) {
      split($2, coords, ":")
      file1 = sprintf("%s/%s_%s_ACRs_1.bed", tmpdir, speciesA, speciesB)
      file2 = sprintf("%s/%s_%s_ACRs_2.bed", tmpdir, speciesA, speciesB)
      if ($10 > $9) {
        print coords[3], coords[4] + $9 - $7, coords[4] + $9 - $7 + 170, $1, ".", "+" > file1
        print coords[3], coords[4] + $10 - $8, coords[4] + $10 - $8 + 170, $1, ".", "+" > file2
      } else {
        print coords[3], coords[4] + $9 + $7 - 170, coords[4] + $9 + $7, $1, ".", "-" > file1
        print coords[3], coords[4] + $10 + $8 - 170, coords[4] + $10 + $8, $1, ".", "-" > file2
      }
    }
  }' ${TMPDIR}/${SPECIES_B}_${SPECIES_A}.blast


### keep only the first entry for each ACR ###
awk -v OFS='\t' '{if ($4 != lastID) print; lastID = $4}' ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_1.bed \
  | sort-bed - \
  > ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_1_fixed.bed
awk -v OFS='\t' '{if ($4 != lastID) print; lastID = $4}' ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_2.bed \
  | sort-bed - \
  > ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_2_fixed.bed


### get sequences ###
bedtools getfasta -s -name -fi ${ANNOTATION}/${SPECIES_A}_genome.fa -bed ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_1_fixed.bed > ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_1.fa
bedtools getfasta -s -name -fi ${ANNOTATION}/${SPECIES_A}_genome.fa -bed ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_2_fixed.bed > ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_2.fa


### find species 1 sequences with lowest Levenshtein distance to species 2 sequence ###
if [ "${SPECIES_A}" == "Sorghum" ]; then
  SPSHORT="Sb"
elif [ "${SPECIES_A}" == "Maize" ]; then
  SPSHORT="Zm"
else
  echo "unknown species: ${SPECIES_A}"
  exit 1
fi 

Rscript subscripts/find_best_seq.r ${SPSHORT} ${OUTDIR}/${SPECIES_B}_ACRs_final.fa ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_1.fa ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_2.fa \
  > ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs.tsv


### get cross-species ACRs which overlap the normal ACRs by less than 80% ###
awk -v FS='\t' -v OFS='\t' 'NR > 1 && $0 != "" {print $2, $3, $4, $1, $6, $5}' ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs.tsv \
  | sort-bed - \
  | bedops -n 80.1% - ${OUTDIR}/${SPECIES_A}_ACRs_final.bed \
  > ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_filtered.bed


### get IDs of identical or similar (at least 80% overlap) ACRs ###
bedmap --echo --echo-map-id --echo-map-score --fraction-either 1 --delim '\t' ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_filtered.bed \
  | awk -v OFS='\t' -v speciesA=${SPECIES_A} -v speciesB=${SPECIES_B} -v tmpdir=${TMPDIR} 'BEGIN{file = sprintf("%s/%s_%s_ACRs_identical.txt", tmpdir, speciesA, speciesB)} {if ($7 ~ /;/) print $7, $8 > file; else print}' \
  | bedmap --echo --echo-map-id --echo-map-score --fraction-either 0.801 --delim '\t' - \
  | awk -v OFS='\t' '$9 ~ /;/ {print $9, $10}' \
  | uniq \
  > ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_similar.txt

if [ -f "${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_identical.txt" ]; then
  uniq ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_identical.txt > ${TMPDIR}/tmp
  mv ${TMPDIR}/tmp ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_identical.txt
else
  touch ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_identical.txt
fi


### select ACRs to remove (keep the one with the lower Levenshtein distance) ###
cat ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_identical.txt ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_similar.txt \
  | awk '{
      n = split($1, IDs, ";")
      split($2, dist, ";")
      if (n > 2) print IDs[2]
      else if (dist[1] < dist [2]) print IDs[1]
      else print IDs[1]
    }' \
  > ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_remove.txt


### remove ACRs ###
awk -v OFS='\t' -v speciesA=${SPECIES_A} -v speciesB=${SPECIES_B} -v tmpdir=${TMPDIR} 'BEGIN{
    file = sprintf("%s/%s_%s_ACRs_remove.txt", tmpdir, speciesA, speciesB);
    while (getline < file) {
      remove[$0]++
    }
  } !($4 in remove) {print}' ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_filtered.bed \
  > ${OUTDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_final.bed

bedtools getfasta -s -nameOnly -fi ${ANNOTATION}/${SPECIES_A}_genome.fa -bed ${OUTDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_final.bed \
  | sed 's/([-+])//' \
  | awk -v speciesA=${SPECIES_A} -v speciesB=${SPECIES_B} -v tmpdir=${TMPDIR} '{if (NR % 2 == 1) name = $0; else if ($1 ~ /^[ACGTacgt]*$/) print name"\n"$1; else print substr(name, 2) > (tmpdir"/"speciesA"_"speciesB"_nonACGT.txt")}' \
  > ${OUTDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_final.fa


### remove ACRs with bases other than ACGT from the final bed file ###
if [ -f "${TMPDIR}/${SPECIES_A}_${SPECIES_B}_nonACGT.txt" ]; then
  grep -v -f ${TMPDIR}/${SPECIES_A}_${SPECIES_B}_nonACGT.txt ${OUTDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_final.bed > ${OUTDIR}/tmp
  mv ${OUTDIR}/tmp ${OUTDIR}/${SPECIES_A}_${SPECIES_B}_ACRs_final.bed
fi