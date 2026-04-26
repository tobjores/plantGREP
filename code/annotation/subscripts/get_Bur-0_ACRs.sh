#!/bin/bash

### map ID to Bur-0 coordinates (+500 bp on either side) ###
# initial mapping of Col-0 to Bur-0 coordinates by Kerry
if [ ! -f "${OUTDIR}/Arabidopsis_Col-0_Bur-0_peak_coords.tsv" ]; then
  unpigz -c ${EXTRADIR}/Arabidopsis_Col-0_Bur-0_peak_coords.tsv.gz \
    > ${OUTDIR}/Arabidopsis_Col-0_Bur-0_peak_coords.tsv
fi

# update the ID (has changed while selecting ACRs)
grep -v "(" ${OUTDIR}/Arabidopsis_ACRs.bed \
  | bedmap --echo --echo-map-id --fraction-either 1 --delim '\t' ${OUTDIR}/Arabidopsis_peaks_union.bed - \
  | awk -v OFS='\t' '{print $1, $2, $3, $5}' \
  > ${TMPDIR}/Arabidopsis_peaks_union_newID.bed

awk -v OFS='\t' 'NR > 1 && $6 != "NA" {print substr($1, 4, 4), $2, $3, substr($5, 4, 4), $6, $7}' ${OUTDIR}/Arabidopsis_Col-0_Bur-0_peak_coords.tsv \
  | bedmap --echo --echo-map-id --delim '\t' - ${TMPDIR}/Arabidopsis_peaks_union_newID.bed \
  | awk -v OFS='\t' '{print $4, $5 - 500, $6 + 500, $7}' \
  > ${TMPDIR}/Bur-0_peak_regions.bed


### shift peak regions overlapping the genome borders ###
awk -v OFS='\t' -v annotation=${ANNOTATION} '
  BEGIN{
    while (getline < (annotation"/Bur-0_genome.fa.fai")) {
      split($0, line, "\t");
      genome[line[1]] = line[2]
    }
  } {
    if ($2 < 0) {$3 = $3 - $2; $2 = $2 - $2; print}
    else if ($3 > genome[$1]) {d = $3 - genome[$1]; $2 = $2 - d; $3 = $3 - d; print}
    else {print}
  }' ${TMPDIR}/Bur-0_peak_regions.bed \
  | sort-bed - \
  > ${TMPDIR}/Bur-0_peak_regions_fixed.bed


### merge peak regions linked to the same ID ###
awk -v OFS='\t' 'BEGIN{getline; lastID = $4; laststart = $2; lastline = $0} {if ($4 == lastID) {$2 = laststart; lastline = $0} else {print lastline} lastline = $0; lastID = $4; laststart = $2} END{print lastline}' ${TMPDIR}/Bur-0_peak_regions_fixed.bed \
  > ${TMPDIR}/Bur-0_peak_regions_fixed2.bed


### get Bur-0 sequences ###
bedtools getfasta -name -fi ${ANNOTATION}/Bur-0_genome.fa -bed ${TMPDIR}/Bur-0_peak_regions_fixed2.bed \
  > ${TMPDIR}/Bur-0_peak_regions.fa


### blast Col-0 ACRs against Bur-0 peak regions ###
makeblastdb -in ${TMPDIR}/Bur-0_peak_regions.fa -dbtype nucl -parse_seqids
blastn -db ${TMPDIR}/Bur-0_peak_regions.fa -query ${OUTDIR}/Arabidopsis_ACRs_final.fa -outfmt 6 -out ${TMPDIR}/Col_Bur.blast


### get Bur-0 coordinates for alignments longer than 100 bases with 80% to 98.5% identity ###
awk -v OFS='\t' -v tmpdir=${TMPDIR} '
  $1 == substr($2, 1, index($2, ":") - 1) && $3 < 98.5 && $3 >= 80 && $4 > 100 {
    split($2, coords, ":")
    file1 = (tmpdir"/Bur-0_ACRs_1.bed")
    file2 = (tmpdir"/Bur-0_ACRs_2.bed")
    if ($10 > $9) {
      print coords[3], coords[4] + $9 - $7, coords[4] + $9 - $7 + 170, $1, ".", "+" > file1
      print coords[3], coords[4] + $10 - $8, coords[4] + $10 - $8 + 170, $1, ".", "+" > file2
    } else {
      print coords[3], coords[4] + $9 + $7 - 170, coords[4] + $9 + $7, $1, ".", "-" > file1
      print coords[3], coords[4] + $10 + $8 - 170, coords[4] + $10 + $8, $1, ".", "-" > file2
    }
  }' ${TMPDIR}/Col_Bur.blast


### keep only the first entry for each ACR ###
awk -v OFS='\t' '{if ($4 != lastID) print; lastID = $4}' ${TMPDIR}/Bur-0_ACRs_1.bed \
  | sort-bed - \
  > ${TMPDIR}/Bur-0_ACRs_1_fixed.bed
awk -v OFS='\t' '{if ($4 != lastID) print; lastID = $4}' ${TMPDIR}/Bur-0_ACRs_2.bed \
  | sort-bed - \
  > ${TMPDIR}/Bur-0_ACRs_2_fixed.bed


### get sequences ###
bedtools getfasta -s -name -fi ${ANNOTATION}/Bur-0_genome.fa -bed ${TMPDIR}/Bur-0_ACRs_1_fixed.bed > ${TMPDIR}/Bur-0_ACRs_1.fa
bedtools getfasta -s -name -fi ${ANNOTATION}/Bur-0_genome.fa -bed ${TMPDIR}/Bur-0_ACRs_2_fixed.bed > ${TMPDIR}/Bur-0_ACRs_2.fa


### find Bur-0 sequences with lowest Levenshtein distance to Col-0 sequence ###
Rscript subscripts/find_best_seq.r Bur-0 ${OUTDIR}/Arabidopsis_ACRs_final.fa ${TMPDIR}/Bur-0_ACRs_1.fa ${TMPDIR}/Bur-0_ACRs_2.fa \
  > ${TMPDIR}/Bur-0_ACRs.tsv


### convert data to .bed file ###
awk -v OFS='\t' 'NR > 1 && $0 != "" {print $2, $3, $4, $1, $6, $5}' ${TMPDIR}/Bur-0_ACRs.tsv \
  | sort-bed - \
  > ${TMPDIR}/Bur-0_ACRs_best.bed


### get IDs of identical or similar (at least 80% overlap) ACRs ###
bedmap --echo --echo-map-id --echo-map-score --fraction-either 1 --delim '\t' ${TMPDIR}/Bur-0_ACRs_best.bed \
  | awk -v OFS='\t' -v tmpdir=${TMPDIR} 'BEGIN{file = tmpdir"/Bur-0_ACRs_identical.txt"} {if ($7 ~ /;/) print $7, $8 > file; else print}' \
  | bedmap --echo --echo-map-id --echo-map-score --fraction-either 0.801 --delim '\t' - \
  | awk -v OFS='\t' '$9 ~ /;/ {print $9, $10}' \
  | uniq \
  > ${TMPDIR}/Bur-0_ACRs_similar.txt

if [ -f "${TMPDIR}/Bur-0_ACRs_identical.txt" ]; then
  uniq ${TMPDIR}/Bur-0_ACRs_identical.txt > ${TMPDIR}/tmp
  mv ${TMPDIR}/tmp ${TMPDIR}/Bur-0_ACRs_identical.txt
else
  touch ${TMPDIR}/Bur-0_ACRs_identical.txt
fi


### select ACRs to remove (keep the one with the lower Levenshtein distance) ###
cat ${TMPDIR}/Bur-0_ACRs_identical.txt ${TMPDIR}/Bur-0_ACRs_similar.txt \
  | awk '{
      n = split($1, IDs, ";")
      split($2, dist, ";")
      if (n > 2) print IDs[2]
      else if (dist[1] < dist [2]) print IDs[1]
      else print IDs[1]
    }' \
  > ${TMPDIR}/Bur-0_ACRs_remove.txt


### remove ACRs ###
awk -v OFS='\t' -v tmpdir=${TMPDIR} 'BEGIN{
    file = (tmpdir"/Bur-0_ACRs_remove.txt")
    while (getline < file) {
      remove[$0]++
    }
  } !($4 in remove) {print}' ${TMPDIR}/Bur-0_ACRs_best.bed \
  > ${OUTDIR}/Bur-0_ACRs_final.bed

bedtools getfasta -s -nameOnly -fi ${ANNOTATION}/Bur-0_genome.fa -bed ${OUTDIR}/Bur-0_ACRs_final.bed \
  | sed 's/([-+])//' \
  | awk -v tmpdir=${TMPDIR} '{if (NR % 2 == 1) name = $0; else if ($1 ~ /^[ACGTacgt]*$/) print name"\n"$1; else print substr(name, 2) > (tmpdir"/Bur-0_nonACGT.txt")}' \
  > ${OUTDIR}/Bur-0_ACRs_final.fa


### remove ACRs with bases other than ACGT from the final bed file ###
if [ -f "${TMPDIR}/Bur-0_nonACGT.txt" ]; then
  grep -v -f ${TMPDIR}/Bur-0_nonACGT.txt ${OUTDIR}/Bur-0_ACRs_final.bed > ${OUTDIR}/tmp
  mv ${OUTDIR}/tmp ${OUTDIR}/Bur-0_ACRs_final.bed
fi