#!/bin/bash

# usage:
# extract_barcode.sh <basename of reads>

# requires:
# - PANDAseq (version 2.11)

### define variables ###
if [[ ${1} == *"/"* ]]; then
  SUBDIR=`dirname ${1}`
  SUBDIR="/"${SUBDIR}
  SAMPLE=`basename ${1}`
else
  SUBDIR=""
  SAMPLE=${1}
fi

READDIR="../../reads/barcodes${SUBDIR}"
OUTDIR="../../data/barcode_counts${SUBDIR}"
mkdir -p ${OUTDIR}

if [[ -z "$TMPDIR" ]]; then
  TMPDIR=${OUTDIR}/tmp${RANDOM}
  mkdir -p ${TMPDIR}
fi

BCLENGTH=18


### print temporary directory name ###
echo "reads are at: ${READDIR}"
echo "curent working directory is at: ${OUTDIR}"
echo "temporary directory is at: ${TMPDIR}"


### extracting barcodes ###
R1LEN="$(zcat "${READDIR}/${SAMPLE}_1.fastq.gz" | head -2 | awk 'NR % 2 == 0 {print length($1)}')"
R2LEN="$(zcat "${READDIR}/${SAMPLE}_2.fastq.gz" | head -2 | awk 'NR % 2 == 0 {print length($1)}')"

if [[ ${R1LEN} -gt ${BCLENGTH} ]]; then
  echo "trimming read 1 ($(date +"%T"))"

  zcat ${READDIR}/${SAMPLE}_1.fastq.gz | awk -v BCLEN=${BCLENGTH} '{if (NR % 2 == 0) print substr($1, 1, BCLEN); else print}' | pigz > ${TMPDIR}/${SAMPLE}_1_trimmed.fastq.gz

  READ1=${TMPDIR}/${SAMPLE}_1_trimmed.fastq.gz
else
  READ1=${READDIR}/${SAMPLE}_1.fastq.gz
fi

if [[ ${R2LEN} -gt ${BCLENGTH} ]]; then
  echo "trimming read 2 ($(date +"%T"))"

  zcat ${READDIR}/${SAMPLE}_2.fastq.gz | awk -v BCLEN=${BCLENGTH} '{if (NR % 2 == 0) print substr($1, 1, BCLEN); else print}' | pigz > ${TMPDIR}/${SAMPLE}_2_trimmed.fastq.gz

  READ2=${TMPDIR}/${SAMPLE}_2_trimmed.fastq.gz
else
  READ2=${READDIR}/${SAMPLE}_2.fastq.gz
fi


### assemble barcodes ###
echo "assembling barcodes ($(date +"%T"))"

mkdir -p ${TMPDIR}/panda_logs

pandaseq -N -f ${READ1} -r ${READ2} -d bfsrkm -G ${TMPDIR}/panda_logs/log_${SAMPLE}_barcodes.txt.bz2 -w ${TMPDIR}/${SAMPLE}_barcodes.fa


### count assembled barcodes of correct length ###
awk -v BCLEN=${BCLENGTH} 'NR % 2 == 0 {if (length($1) == BCLEN) barcodes[$1]++} END{for (bc in barcodes) print barcodes[bc], bc}' ${TMPDIR}/${SAMPLE}_barcodes.fa | pigz -f > ${OUTDIR}/barcode_${SAMPLE}.count.gz
