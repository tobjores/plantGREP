#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### get annotations for synthesized sequences  from Arabidopsis (Col-0 and Bur-0), maize, sorghum, and tomato ###
for SPECIES in "Arabidopsis" "Bur-0" "Maize" "Sorghum" "Tomato"
do
  # combine all sequences of a given species
  cat ${OUTDIR}/${SPECIES}_*_final.bed \
    | sort-bed - \
    > ${TMPDIR}/${SPECIES}_sequences.bed

  # get perbase cutcount files
  NFILES=1

  if [ "${SPECIES}" = "Arabidopsis" ]
  then
    unpigz -c ${EXTRADIR}/Arabidopsis_Col-0_perbase.bed.gz \
      | sed 's/chr0\?//g' \
      > ${TMPDIR}/${SPECIES}_perbase_${NFILES}.bed
  elif [ "${SPECIES}" = "Bur-0" ]
  then
    unpigz -c ${EXTRADIR}/Arabidopsis_Bur-0_perbase.bed.gz \
      | sed 's/chr0\?//g' \
      > ${TMPDIR}/${SPECIES}_perbase_${NFILES}.bed
  else
    for FILE in ${EXTRADIR}/${SPECIES}*perbase*
    do
      unpigz -c ${FILE} \
        | sed 's/chr0\?//g' \
        > ${TMPDIR}/${SPECIES}_perbase_${NFILES}.bed

      let "NFILES=NFILES+1"
    done
    let "NFILES=NFILES-1"
  fi

  # get sum of cutcounts for synthesized sequences
  for ID in $(seq 1 ${NFILES})
  do
    bedmap --sum ${TMPDIR}/${SPECIES}_sequences.bed ${TMPDIR}/${SPECIES}_perbase_${NFILES}.bed \
      > ${TMPDIR}/${SPECIES}_cutcount_${ID}.tsv
  done

  paste ${TMPDIR}/${SPECIES}_sequences.bed ${TMPDIR}/${SPECIES}_cutcount_*.tsv \
    | awk -v OFS='\t' 'function colsum(first_col, last_col) {
        if (first_col == 0) first_col = 1
        if (last_col == 0) last_col = NF
        sum = 0
        for (i = first_col; i <= last_col; i++) sum += $i
        return sum
      } BEGIN{
        print "chromosome", "start", "end", "id", "cutcount"
      } {print $1, $2, $3, $4, colsum(5)}' \
    > ${TMPDIR}/${SPECIES}_wCC.tsv

  # get midpoint coordinates
  awk 'NR > 1 {print $1, $2 + 85, $3 - 84}' ${TMPDIR}/${SPECIES}_wCC.tsv \
    > ${TMPDIR}/${SPECIES}_center.bed
  
  # get region in which the midpoint resides
  bedmap --echo-map-id ${TMPDIR}/${SPECIES}_center.bed ${ANNOTATION}/${SPECIES}_features.bed \
    | awk 'BEGIN{
        print "region"
      } {
        if ($1 == "") region = "intergenic"
        else if ($1 ~ /CDS/) region = "CDS"
        else if ($1 ~ /five_prime_UTR/) region = "5'"'"'UTR"
        else if ($1 ~ /three_prime_UTR/) region = "3'"'"'UTR"
        else if ($1 ~ /ncRNA/) region = "ncRNA"
        else region = "intron"
        print region
      }' \
    > ${TMPDIR}/${SPECIES}_regions.tsv

  # get closest TSS
  closest-features --delim '\t' --closest --dist ${TMPDIR}/${SPECIES}_center.bed ${ANNOTATION}/${SPECIES}_TSS.bed \
    | awk -v OFS='\t' 'BEGIN{print "closest_TSS", "TSS_dist"} {
        if ($9 == "+") {$10 = $10 * -1}
        print $7"("$9")", $10
      }' \
    > ${TMPDIR}/${SPECIES}_TSS.tsv

  # get species abbreviation
  if [ "${SPECIES}" == "Arabidopsis" ] || [ "${SPECIES}" == "Bur-0" ]; then
    SPSHORT="At"
  elif [ "${SPECIES}" == "Maize" ]; then
    SPSHORT="Zm"
  elif [ "${SPECIES}" == "Sorghum" ]; then
    SPSHORT="Sb"
  elif [ "${SPECIES}" == "Tomato" ]; then
    SPSHORT="Sl"
  else
    echo "unknown species: ${SPECIES}"
    exit 1
  fi
  
  # get overlap with histone ChIP peaks (Arabidopsis only)
  if [ "${SPECIES}" == "Arabidopsis" ]; then
    zgrep "Col-0" ${EXTRADIR}/remap2022_histone_nr_macs2_TAIR10_v1_0.bed.gz \
      | awk -F '	|:' -v OFS='\t' '{print $1, $2, $3, $4}' \
      | bedmap --range 500 --echo-map-id-uniq ${TMPDIR}/${SPECIES}_center.bed - \
      | awk 'BEGIN{print "histones"} {print}' \
      > ${TMPDIR}/${SPECIES}_histones.tsv
  else
    echo "histones" > ${TMPDIR}/${SPECIES}_histones.tsv
  fi

  # combine annotations
  paste ${TMPDIR}/${SPECIES}_wCC.tsv ${TMPDIR}/${SPECIES}_regions.tsv ${TMPDIR}/${SPECIES}_TSS.tsv ${TMPDIR}/${SPECIES}_histones.tsv \
    | awk -v OFS='\t' -v spshort=${SPSHORT} 'NR == 1 {print $0, "species"} NR > 1 {print $0, spshort}' \
    > ${TMPDIR}/${SPECIES}_annotation.tsv
done

### combine species annotations ###
awk -v OFS='\t' 'NR == 1 {print} FNR > 1 {print}' ${TMPDIR}/*_annotation.tsv \
  | pigz \
  > ${OUTDIR}/tamsACR_annotation.tsv.gz
