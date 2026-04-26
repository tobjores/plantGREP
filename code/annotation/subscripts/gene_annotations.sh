#!/bin/bash

### define variables ###
OUTDIR="${1}"
mkdir -p ${OUTDIR}


### download genomes and annotations ###
if [ ! -f "${OUTDIR}/Arabidopsis_genome.fa" ]; then
  curl https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas \
    > ${OUTDIR}/Arabidopsis_genome.fa
  samtools faidx ${OUTDIR}/Arabidopsis_genome.fa
fi
if [ ! -f "${OUTDIR}/Arabidopsis_annotation.gff3" ]; then
  curl https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/archived/Araport11_GFF3_genes_transposons.Mar92021.gff.gz \
    | unpigz \
    | grep -v "^#" \
    | sed 's/^Chr//g' \
    > ${OUTDIR}/Arabidopsis_annotation.gff3
fi

if [ ! -f "${OUTDIR}/Bur-0_genome.fa" ]; then
  curl http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/fasta/UNMASKED/bur_0.v7.fas \
    | sed 's/Chr//g' \
    > ${OUTDIR}/Bur-0_genome.fa
  samtools faidx ${OUTDIR}/Bur-0_genome.fa
fi
if [ ! -f "${OUTDIR}/Bur-0_annotation.gff3" ]; then
  curl http://mtweb.cs.ucl.ac.uk/mus/www/19genomes/annotations/consolidated_annotation_9.4.2011/gene_models/consolidated_annotation.Bur_0.gff3.bz2 \
    | bunzip2 \
    | grep -v "^#" \
    | sed 's/^Chr//g' \
    > ${OUTDIR}/Bur-0_annotation.gff3
fi

if [ ! -f "${OUTDIR}/Maize_genome.fa" ]; then
  curl ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/zea_mays/dna/Zea_mays.B73_RefGen_v4.dna.toplevel.fa.gz \
    | unpigz \
    > ${OUTDIR}/Maize_genome.fa
  samtools faidx ${OUTDIR}/Maize_genome.fa
fi
if [ ! -f "${OUTDIR}/Maize_annotation.gff3" ]; then
  curl ftp://ftp.ensemblgenomes.org/pub/plants/release-50/gff3/zea_mays/Zea_mays.B73_RefGen_v4.50.chr.gff3.gz \
    | unpigz \
    | grep -v "^#" \
    > ${OUTDIR}/Maize_annotation.gff3
fi

if [ ! -f "${OUTDIR}/Maize_dTSSs_shoot.gff3" ]; then
  unpigz -c ${EXTRADIR}/Maize_dTSSs_shoot.gff3.gz \
    > ${OUTDIR}/Maize_dTSSs_shoot.gff3
fi
if [ ! -f "${OUTDIR}/Maize_dTSSs_root.gff3" ]; then
  unpigz -c ${EXTRADIR}/Maize_dTSSs_root.gff3.gz \
    > ${OUTDIR}/Maize_dTSSs_root.gff3
fi

if [ ! -f "${OUTDIR}/Sorghum_genome.fa" ]; then
  curl ftp://ftp.ensemblgenomes.org/pub/plants/release-50/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz \
    | unpigz \
    > ${OUTDIR}/Sorghum_genome.fa
  samtools faidx ${OUTDIR}/Sorghum_genome.fa
fi
if [ ! -f "${OUTDIR}/Sorghum_annotation.gff3" ]; then
  curl ftp://ftp.ensemblgenomes.org/pub/plants/release-50/gff3/sorghum_bicolor/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.50.chr.gff3.gz \
    | unpigz \
    | grep -v "^#" \
    > ${OUTDIR}/Sorghum_annotation.gff3
fi

# the tomato genome and annotations were supplied by the Lippman lab
if [ ! -f "${OUTDIR}/Tomato_genome.fa" ]; then
  unpigz -c ${EXTRADIR}/Tomato_genome.fa.gz \
    > ${OUTDIR}/Tomato_genome.fa
  samtools faidx ${OUTDIR}/Tomato_genome.fa
fi
if [ ! -f "${OUTDIR}/Tomato_annotation.gff3" ]; then
  unpigz -c ${EXTRADIR}/Tomato_annotation.gff3.gz \
    > ${OUTDIR}/Tomato_annotation.gff3
fi


### get gene coordinates for Arabidopsis (Col-0 and Bur-0), maize, sorghum, and tomato ###
for SPECIES in "Arabidopsis" "Bur-0" "Maize" "Sorghum" "Tomato"
do

  if [ "${SPECIES}" = "Arabidopsis" ]
  then
    awk '{if($3 == "gene" && $1 ~ /[0-9]+/) print}' ${OUTDIR}/${SPECIES}_annotation.gff3 \
      | grep "type=mirna" \
      | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
      > ${OUTDIR}/${SPECIES}_miRNA_genes.gff3
  elif [ "${SPECIES}" = "Bur-0" ]
  then
    awk -v OFS='\t' '{if($3 == "miRNA" && $1 ~ /[0-9]+/) print $1, $2, $3, $4, $5, $6, $7, $8, $9}' ${OUTDIR}/${SPECIES}_annotation.gff3 \
      > ${OUTDIR}/${SPECIES}_miRNA_genes.gff3   
  elif [ "${SPECIES}" = "Maize" ] || [ "${SPECIES}" = "Sorghum" ]
  then
    awk '{if($3 == "ncRNA_gene" && $1 ~ /[0-9]+/) print}' ${OUTDIR}/${SPECIES}_annotation.gff3 \
      | grep "type=pre_miRNA" \
      | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
      > ${OUTDIR}/${SPECIES}_miRNA_genes.gff3
  else
    touch ${OUTDIR}/${SPECIES}_miRNA_genes.gff3
  fi

  if [ "${SPECIES}" = "Tomato" ] || [ "${SPECIES}" = "Bur-0" ]
  then
    awk -v FS='\t' -v OFS='\t' '{if($3 == "gene" && $1 ~ /[0-9]+/) print $1, $2, $3, $4, $5, $6, $7, $8, $9}' ${OUTDIR}/${SPECIES}_annotation.gff3 \
      > ${OUTDIR}/${SPECIES}_protein_coding_genes.gff3
  else
    awk '{if($3 == "gene" && $1 ~ /[0-9]+/) print}' ${OUTDIR}/${SPECIES}_annotation.gff3 \
      | grep "type=protein_coding" \
      | awk -v OFS='\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' \
      > ${OUTDIR}/${SPECIES}_protein_coding_genes.gff3
  fi
    
  if [ ${SPECIES} = 'Maize' ]
  then
    Rscript subscripts/correct_Maize_TSSs.r ${OUTDIR}/Maize_dTSSs_shoot.gff3 ${OUTDIR}/Maize_dTSSs_root.gff3 ${OUTDIR}/Maize_protein_coding_genes.gff3
  fi

  cat ${OUTDIR}/${SPECIES}_protein_coding_genes.gff3 ${OUTDIR}/${SPECIES}_miRNA_genes.gff3 \
    | awk -v FS='[\t;]' -v OFS='\t' '$4 < $5 {sub("ID=(Gene:|gene:)?", "", $9); sub("\\.[0-9]+", "", $9); print $1, $4 - 1, $5, $9, $6, $7}' \
    | sort-bed - \
    > ${OUTDIR}/${SPECIES}_genes.bed

  awk -v OFS='\t' '{if ($6 == "+") {$3 = $2 + 1} else {$2 = $3 - 1}; print}' ${OUTDIR}/${SPECIES}_genes.bed \
    | sort-bed - \
    > ${OUTDIR}/${SPECIES}_TSS.bed

  awk -v OFS='\t' '{if ($6 == "+") {$3 = $2 + 5; $2 = $2 - 60} else {$2 = $3 - 5; $3 = $3 + 60}; print}' ${OUTDIR}/${SPECIES}_genes.bed \
    | sort-bed - \
    > ${OUTDIR}/${SPECIES}_corePro.bed

  awk -v FS='\t' -v OFS='\t' '$3 ~ /^gene|CDS|UTR|[^m]RNA/ {print $1, $4 -1, $5, $3, ".", $7}' ${OUTDIR}/${SPECIES}_annotation.gff3 \
    | sort-bed - \
    > ${OUTDIR}/${SPECIES}_features.bed

done
