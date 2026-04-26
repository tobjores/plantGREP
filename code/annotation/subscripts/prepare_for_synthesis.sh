#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### concatenate final sequence files ###
# Arabidopsis, maize and sorghum sequences
cat $(ls ${OUTDIR}/*final.fa | grep -v "Tomato") > ${TMPDIR}/amsACR_array.fa

# Tomato sequences
cat $(ls ${OUTDIR}/*final.fa | grep "Tomato") > ${TMPDIR}/tACR_array.fa


### blast sequences against each other to find very similar ones (amsACR only) ###
makeblastdb -in ${TMPDIR}/amsACR_array.fa -dbtype nucl -parse_seqids
blastn -db ${TMPDIR}/amsACR_array.fa -query ${TMPDIR}/amsACR_array.fa -ungapped -outfmt 6 -out ${TMPDIR}/amsACR_array.blast


### for sequence pairs with less than three mismatches, remove one (amsACR only) ###
awk '$1 != $2 && ($4-$5) > 167' ${TMPDIR}/amsACR_array.blast > ${TMPDIR}/amsACR_array_critical.blast

awk '{if (!($2 in remove)) {remove[$1]; print ">"$1}}' ${TMPDIR}/amsACR_array_critical.blast > ${TMPDIR}/amsACR_array_remove.txt

awk -v OFS='\t' -v tmpdir=${TMPDIR} 'BEGIN{while (getline < (tmpdir"/amsACR_array_remove.txt")) {remove[$0]++}} {if ($0 in remove) getline; else print}' ${TMPDIR}/amsACR_array.fa > ${TMPDIR}/amsACR_array_filtered.fa


### remove duplicated sequences (tACR only) ###
awk '{name = $1; getline; if ($1 in seqs) {if (substr(name, 1, 6) == ">Solyc") print seqs[$1]; else print name}; seqs[$1] = name}' ${TMPDIR}/tACR_array.fa > ${TMPDIR}/tACR_array_remove.txt

awk -v OFS='\t' -v tmpdir=${TMPDIR} 'BEGIN{while (getline < (tmpdir"/tACR_array_remove.txt")) {remove[$0]++}} {if ($0 in remove) getline; else print}' ${TMPDIR}/tACR_array.fa > ${TMPDIR}/tACR_array_filtered.fa


### final processing ###
for ARRAY in "amsACR" "tACR"
do

  # mutate BsaI recognition sites
  awk -v OFS='\t' -v outdir=${OUTDIR} -v array=${ARRAY} '{
      NAME=$1
      getline
      SEQ=sprintf("%sC", $1)
      MUT=""
      if(SEQ ~ /(GGTCTC)|(GAGACC)/) {
        while(index(SEQ, "GGTCTC") > 0) {MUT=MUT "T" index(SEQ, "GGTCTC") + 2 ">A;"; sub(/GGTCTC/, "GGACTC", SEQ)};
        while(index(SEQ, "GAGACC") > 0) {MUT=MUT "A" index(SEQ, "GAGACC") + 3 ">T;"; sub(/GAGACC/, "GAGTCC", SEQ)};
        print substr(NAME, 2), substr(MUT, 1, length(MUT) - 1) > (outdir"/"array"_array_mutations.tsv")
      }
      print NAME"\n"substr(SEQ, 1, 170)
    }' ${TMPDIR}/${ARRAY}_array_filtered.fa \
    > ${OUTDIR}/${ARRAY}_array_sequences.fa

  # add cloning adapters
  awk '{if (NR % 2 == 0) $1 = "GAGCGGTCTCCACTC"$1"CTGTAGAGACCGGGC"; print}' ${OUTDIR}/${ARRAY}_array_sequences.fa > ${OUTDIR}/${ARRAY}_array_sequences_with_adapter.fa

done
