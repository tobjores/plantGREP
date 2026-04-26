#!/bin/bash

# requires:
#   - bedops (version 2.4.35)
#   - bedtools (version 2.29.2)
#   - samtools (version 1.10)
#   - blast (version 2.11.0)
#   - R (version 4.0.0) with libraries:
#     - dplyr (version 1.0.2)
#     - readr (version 1.4.0)
#     - stringr (version 1.4.0)
#     - tidyr (version 1.1.2)
#     - tibble (version 3.0.4)
#     - readxl (version 1.3.1)
#     - Biostrings (version 2.56.0)

### define variables ###
# main output directory
OUTDIR="../../data/annotation"
mkdir -p ${OUTDIR}

if [[ -z "$TMPDIR" ]]; then
  TMPDIR=${OUTDIR}/tmp${RANDOM}
  mkdir -p ${TMPDIR}
fi

# directory with required extra files
EXTRADIR="../../data/extra_files"

# directory for gene annotations (and genomes)
ANNOTATION="${OUTDIR}/gene_annotations"

# directory for ACR sequences
ACRDIR="${OUTDIR}/ACRs"

# directory for reference sequences
REFSEQ="../../data/refseq"
mkdir -p ${REFSEQ}


### gene annotations ###
. subscripts/gene_annotations.sh ${ANNOTATION}


### ACR sequences ###
. subscripts/Arabidopsis_ACRs.sh ${ACRDIR}
. subscripts/Maize_ACRs.sh ${ACRDIR}
. subscripts/Sorghum_ACRs.sh ${ACRDIR}
. subscripts/CNS_ACRs_Maize_Sorghum.sh ${ACRDIR}
. subscripts/Tomato_ACRs.sh ${ACRDIR}


### control and extra sequences ###
. subscripts/nonACR_controls.sh ${ACRDIR}
. subscripts/Enhancers+Insulators.sh ${ACRDIR}


### prepare sequences for array synthesis ###
. subscripts/prepare_for_synthesis.sh ${ACRDIR}


### merge and copy final sequence files to reference sequence directory ###
cat ${ACRDIR}/*_array_sequences.fa | pigz > ${REFSEQ}/tamsACR_sequences.fa.gz


### annotate sequences ###
. subscripts/ACR_annotations.sh ${ACRDIR}


### calculate genomic GC content ###
. subscripts/genome_GC.sh
