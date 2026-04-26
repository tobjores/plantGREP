### process arguments ###
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop('Three files are required: 1. experimental TSSs from shoots, 2. experimental TSSs from roots, 3. the file to correct.')
} else if (length(args) > 3) {
  warning(paste0('More than three arguments passed. Ignoring arguments 4 - ', length(args), '.'))
  shootfile <- args[1]
  rootfile <- args[2]
  correctfile <- args[3]
} else {
  shootfile <- args[1]
  rootfile <- args[2]
  correctfile <- args[3]
}


### test if data file exists ###
if (! file.exists(shootfile)) {
  stop(paste0('The file "', shootfile, '" does not exist.'))
}
if (! file.exists(rootfile)) {
  stop(paste0('The file "', rootfile, '" does not exist.'))
}
if (! file.exists(correctfile)) {
  stop(paste0('The file "', correctfile, '" does not exist.'))
}


### load libraries ###
library(dplyr)
library(tidyr)
library(readr)
library(tibble)


### load files ###
# the 'eglab_*' files were obtained from the Grotewold lab; Mejía-Guerra et al., 2015, Plant Cell, doi: 10.1105/tpc.15.00630
shoot <- read_tsv(shootfile, comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = substr(annotation, 19, 32)
  )

root <- read_tsv(rootfile, comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = substr(annotation, 19, 32)
  )

reference <- read_tsv(correctfile, comment = '#', col_names = c('chromosome', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'annotation')) %>%
  mutate(
    gene = substr(annotation, 9, 22)
  )


### merge shoot and root dTSSs
correctedTSS <- union(shoot, root)


### get reference TSSs ###
referenceTSS <- reference %>%
  mutate(
    TSS = if_else(strand == '+', start, stop)
  ) %>%
  select(gene, TSS) %>%
  deframe()


### extract dTSSs linked to 2+ genes and link them to the closest reference TSS ###
multipleTSS <- correctedTSS %>%
  filter(grepl(',', annotation, fixed = TRUE)) %>%
  rename('gene1' = gene) %>%
  mutate(
    gene2 = substr(annotation, 34, 47),
    gene3 = substr(annotation, 49, 62),
    gene4 = substr(annotation, 64, 77)
  ) %>%
  pivot_longer(
    cols = starts_with('gene'),
    values_to = 'gene',
    names_to = 'gene_id',
    names_prefix = 'gene'
  ) %>%
  filter(grepl('Zm00001d', gene, fixed = TRUE)) %>%
  mutate(
    refTSS = referenceTSS[gene]
  ) %>%
  group_by(across(c(-gene, -gene_id, -refTSS))) %>%
  summarise(
    gene = first(gene[abs(stop - refTSS) == min(abs(stop - refTSS))])
  ) %>%
  ungroup() %>%
  select(annotation, gene) %>%
  deframe()

correctedTSS <- correctedTSS %>%
  mutate(
    gene = if_else(annotation %in% names(multipleTSS), multipleTSS[annotation], gene)
  ) %>%
  select(gene, stop) %>%
  deframe()


### export names of genes with a dTSS in the Grotewold data ###
outdir <- dirname(correctfile)
write_lines(names(correctedTSS), paste0(outdir, '/Maize_corrected_TSS.txt'))


### correct TSS annotation in B73 reference annotation ###
TSS <- reference %>%
  mutate(
    start = if_else(strand == '+' & gene %in% names(correctedTSS), correctedTSS[gene], start),
    stop = if_else(strand == '-' & gene %in% names(correctedTSS), correctedTSS[gene], stop)
  ) %>% 
  select(-gene)


### save corrected annotation to file ###
write_tsv(TSS, correctfile, col_names = FALSE)
