#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### get Maize-Sorghum CNS alignments (https://doi.org/10.1101/gr.266528.120) ###
if [ ! -f "${OUTDIR}/CNS_Maize_Sorghum.sam" ]; then
  wget -nv -O ${OUTDIR}/CNS_Maize_Sorghum.sam https://figshare.com/ndownloader/files/26641469
fi

### get coordinates from .sam file (ignore chromosomes other than 1-10) ###
# uses .bed file-like coordinates (i.e. start is 0-based and end is 1-based)
awk -v OFS='\t' '
  function ref_len(cigar) {
    reflen = 0
    split(cigar, cigar_number, "[HDMI]", cigar_type)
    for (i in cigar_type) {
      if (cigar_type[i] ~ /[DM]/) reflen += cigar_number[i]
    }
    return reflen
  }
  BEGIN{print "id", "Sb_chr", "Sb_start", "Sb_end", "Zm_chr", "Zm_start", "Zm_end"; id = 0}
  $1 ~ /^[[:digit:]]+$/ && $3 ~ /^[[:digit:]]+$/ {
    id += 1
    split($6, cigar, "H")
    print "CNS-"id, $1, cigar[1] - 1, cigar[1] + length($10) - 1, $3, $4 - 1, $4 + ref_len(cigar[2]) - 1
  }' ${OUTDIR}/CNS_Maize_Sorghum.sam \
  > ${OUTDIR}/CNS_Maize_Sorghum.tsv


### create .bed files for each species ###
awk -v OFS='\t' 'NR > 1 {print $5, $6, $7, $1}' ${OUTDIR}/CNS_Maize_Sorghum.tsv | sort-bed - > ${OUTDIR}/Maize_CNSs.bed
awk -v OFS='\t' 'NR > 1 {print $2, $3, $4, $1}' ${OUTDIR}/CNS_Maize_Sorghum.tsv | sort-bed - > ${OUTDIR}/Sorghum_CNSs.bed


### get cross-species ACRs ###
. subscripts/get_cross-species_ACRs.sh Maize Sorghum
. subscripts/get_cross-species_ACRs.sh Sorghum Maize
