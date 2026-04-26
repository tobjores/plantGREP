#!/bin/bash

### set up output file ###
echo -e "species\tGC" > ${OUTDIR}/genome_GC.tsv

### calculate GC content ###
for SPECIES in "Arabidopsis" "Tomato" "Maize" "Sorghum"
do
  awk -v OFS='\t' '$1 ~ "^[0-9]+$" {print $1, 0, $2}' ${ANNOTATION}/${SPECIES}_genome.fa.fai  | bedtools nuc -fi ${ANNOTATION}/${SPECIES}_genome.fa -bed - | awk -v OFS='\t' -v species=${SPECIES} 'NR > 1 {AT += $6 + $9; GC += $7 + $8} END {print species, GC / (GC + AT)}' >> ${OUTDIR}/genome_GC.tsv
done
