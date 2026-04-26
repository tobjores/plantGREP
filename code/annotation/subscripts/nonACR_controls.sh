#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### pick non-ACR controls for each species ###
for SPECIES in "Arabidopsis" "Maize" "Sorghum" "Tomato"
do

  # get ACRs
  cat ${OUTDIR}/${SPECIES}*final.bed | cut -f 1-4 | sort-bed - > ${TMPDIR}/${SPECIES}_exclude.bed

  # get species abbreviation and number of controls
  if [ "${SPECIES}" == "Arabidopsis" ]; then
    SPSHORT="At"
    NCONTROLS="1673"
    NLINES="2000"
  elif [ "${SPECIES}" == "Maize" ]; then
    SPSHORT="Zm"
    NCONTROLS="1561"
    NLINES="2000"
  elif [ "${SPECIES}" == "Sorghum" ]; then
    SPSHORT="Sb"
    NCONTROLS="1066"
    NLINES="2000"
  elif [ "${SPECIES}" == "Tomato" ]; then
    SPSHORT="Sl"
    NCONTROLS="2164"
    NLINES="3000"
  else
    echo "unknown species: ${SPECIES}"
    exit 1
  fi
  
  # create a genome file
  if [ "${SPECIES}" == "Tomato" ]; then
    awk -v OFS='\t' '$1 ~ /^[[:digit:]]+$/ {print $1, $2}' ${ANNOTATION}/${SPECIES}_genome.fa.fai > ${TMPDIR}/${SPECIES}.genome
  else
    awk -v OFS='\t' '$1 ~ /^[[:digit:]]$/ {print $1, $2}' ${ANNOTATION}/${SPECIES}_genome.fa.fai > ${TMPDIR}/${SPECIES}.genome
  fi
  
  # shuffle ACRs
  head -n ${NLINES} ${TMPDIR}/${SPECIES}_exclude.bed \
    | bedtools shuffle -seed 928 -noOverlapping -excl ${TMPDIR}/${SPECIES}_exclude.bed -f 0 -g ${TMPDIR}/${SPECIES}.genome -i - \
    | sort-bed - \
    > ${TMPDIR}/${SPECIES}_shuffledACRs.bed
  
  # get sequences (remove sequences with bases other than ACGT)
  bedtools getfasta -fi ${ANNOTATION}/${SPECIES}_genome.fa -bed ${TMPDIR}/${SPECIES}_shuffledACRs.bed \
    | awk -v OFS='\t' -v species=${SPECIES} -v spshort=${SPSHORT} -v outdir=${OUTDIR} -v ncontrols=${NCONTROLS} 'BEGIN{i = 1} {
          if (i > ncontrols) exit;
          if (NR % 2 == 1) split($1, coords, "[>:-]");
          else if ($0 ~ /^[ACGT]*$/) {print ">"spshort"-sh"i"\n"$0; print coords[2], coords[3], coords[4], spshort"-sh"i > (outdir"/"species"_shuffledACRS_final.bed"); i++}
        }' \
    > ${OUTDIR}/${SPECIES}_shuffledACRs_final.fa

done
