#!/bin/bash

# requires:
# - PANDAseq (version 2.11)
# - bioawk (version 1.0)
# - bowtie (version 2.5.3)
# - R (version 4.4.1)

### define variables ###
NSLOTS=`nproc`

READBASENAME="tamsACR_subassembly"

REFSEQ="../../data/refseq/tamsACR_sequences.fa.gz"

READDIR="../../reads/subassembly/tamsACR"
OUTDIR="../../data/subassembly"
mkdir -p ${READDIR} ${OUTDIR}

if [[ -z "$TMPDIR" ]]; then
  TMPDIR=${OUTDIR}/tmp${RANDOM}
  mkdir -p ${TMPDIR}
fi

CONTROLDIR="../../data/extra_files"


### print temporary directory name ###
# echo "reads will be deposited at: ${READDIR}"
echo "read files are at: ${READDIR}"
echo "final files will be deposited at: ${OUTDIR}"
echo "temporary directory is at: ${TMPDIR}"


### assemble inserts and barcodes ###
echo "assembling reads ($(date +"%T"))"

mkdir -p ${TMPDIR}/panda_logs

pandaseq -N -f ${READDIR}/${READBASENAME}_insert_1.fastq.gz -r ${READDIR}/${READBASENAME}_insert_2.fastq.gz -d bfsrkm -G ${TMPDIR}/panda_logs/log_tamsACR_insert.txt.bz2 -w ${TMPDIR}/tamsACR_insert.fa
pandaseq -N -f ${READDIR}/${READBASENAME}_barcode_1.fastq.gz -r ${READDIR}/${READBASENAME}_barcode_2.fastq.gz -d bfsrkm -G ${TMPDIR}/panda_logs/log_tamsACR_barcode.txt.bz2 -w ${TMPDIR}/tamsACR_barcode.fa


### extract barcodes and split them by library ###
echo "spliting barcodes ($(date +"%T"))"

awk -v tmpdir=${TMPDIR} -v OFS='\t' '{
    sub(/:[ACGTN0+]*;.*$/, "", $1);
    NAME=substr($1, 2);
    getline;
    if (length($1) == 18) {
      if (substr($1, 7, 1)substr($1, 13, 1) == "AG") print NAME, $0 > (tmpdir "/barcodes_pPSntF_tamsACR_35SprAG.tsv");
      else if (substr($1, 7, 1)substr($1, 13, 1) == "GA") print NAME, $0 > (tmpdir "/barcodes_pPSntR_tamsACR_35SprGA.tsv");
    }
  }' ${TMPDIR}/tamsACR_barcode.fa

rm ${TMPDIR}/tamsACR_barcode.fa


### join barcodes and inserts ###
echo "joining barcodes and inserts ($(date +"%T"))"

for LIB in "pPSntF_tamsACR_35SprAG" "pPSntR_tamsACR_35SprGA"
do
  bioawk -c fastx -t -v library=${LIB} -v tmpdir=${TMPDIR} '
      BEGIN{
        while (getline < (tmpdir "/barcodes_" library ".tsv")) {
          split($0, line, "\t")
          barcode[line[1]] = line[2]
        }
        print "barcode", "sequence", "assembly_count"
      } {
        sub(/:[ACGTN0+]*;.*$/, "", $name);
        if ($name in barcode) {counts[barcode[$name]$seq]++; delete barcode[$name]}
      } END{
        for (combo in counts) if (counts[combo] >= 5) print substr(combo, 1,  18), substr(combo, 19), counts[combo]
      }
    ' ${TMPDIR}/tamsACR_insert.fa \
    > ${TMPDIR}/subassembly_${LIB}.tsv
  
  rm ${TMPDIR}/barcodes_${LIB}.tsv

  echo "finished ${LIB} ($(date +"%T"))"
done

rm ${TMPDIR}/tamsACR_insert.fa


### align ACR sequences ###
echo "aligning ACRs ($(date +"%T"))"

# build index #
mkdir -p ${TMPDIR}/bowtie2-index
bowtie2-build ${REFSEQ} ${TMPDIR}/bowtie2-index/index 1>&2

for LIB in "pPSntF_tamsACR_35SprAG" "pPSntR_tamsACR_35SprGA"
do
  # convert subassembly to fasta #
  if [[ ${LIB} == "pPSntR"* ]]; then
    bioawk 'NR > 1 {print ">" $1 ":" $3 "\n" revcomp($2)}' ${TMPDIR}/subassembly_${LIB}.tsv \
     > ${TMPDIR}/${LIB}.fa
  else
    bioawk 'NR > 1 {print ">" $1 ":" $3 "\n" $2}' ${TMPDIR}/subassembly_${LIB}.tsv \
     > ${TMPDIR}/${LIB}.fa
  fi

  # align #
  bowtie2 -p ${NSLOTS} --xeq --no-unal -f --norc -x ${TMPDIR}/bowtie2-index/index -U ${TMPDIR}/${LIB}.fa 2> ${OUTDIR}/alignment_${LIB}.txt > ${TMPDIR}/${LIB}.sam


  ### extract alignment info from sam file ###
  bioawk -c sam -t '
      BEGIN{print "barcode", "assembly_count", "id", "start", "cigar", "sequence", "GC"}
      {split($qname, query, ":"); tmp = $seq; gsub(/A|T/, "", tmp); print query[1], query[2], $rname, $pos, $cigar, $seq, length(tmp)/length($seq)}
    ' ${TMPDIR}/${LIB}.sam \
    > ${TMPDIR}/${LIB}.tsv


  ### convert CIGAR string to variant info ###
  Rscript cigar_to_variant.r ${TMPDIR}/${LIB}.tsv ${REFSEQ}
  mv ${TMPDIR}/${LIB}_variant.tsv ${TMPDIR}/subassembly_${LIB}.tsv

  echo "finished ${LIB} ($(date +"%T"))"
done


### filter subassemblies ###
Rscript filter_subassembly_tamsACR.r ${TMPDIR}/subassembly_pPSntF_tamsACR_35SprAG.tsv ${TMPDIR}/subassembly_pPSntR_tamsACR_35SprGA.tsv ${CONTROLDIR}/barcodes_pPSntF_controls_35SprCC.txt.gz \
  > ${TMPDIR}/subassembly_tamsACR_filtered.tsv


### compress final files and copy to main directory ###
echo "copying files ($(date +"%T"))"

for FILE in ${TMPDIR}/subassembly_*.tsv
do
  FILENAME=`basename ${FILE}`
  pigz -c ${FILE} > ${OUTDIR}/${FILENAME}.gz
done

echo "all done ($(date +"%T"))"
