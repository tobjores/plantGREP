library(tidyverse)
library(universalmotif)
library(RcppRoll)
library(readxl)
library(openxlsx)
library(gprofiler2)

# data loading and preprocessing ------------------------------------------

### load functions to export the plot data ##
source('code/analysis/pggfplot.R')

pggf_config(out_dir = 'figures/rawData/')


### do not print code or preview plots when sourced ###
if (sys.nframe() != 0) {
  pggf_config(print_code = FALSE, preview_plot = FALSE)
  options(readr.show_col_types = FALSE)
}


### define helpers ###
condition_order <- c('light', 'dark', 'warm', 'cold', 'maize')
species_short <- c('At', 'Sl', 'Zm', 'Sb')
species_to_short <- setNames(species_short, c('Arabidopsis', 'Tomato', 'Maize', 'Sorghum'))


### load tamsACR replicate data ###
tamsACR_reps <- readRDS('data/RData/tamsACR_reps.rds') |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )


### load main tamsACR data and add annotation ###
tamsACR_annotation <- read_tsv('data/annotation/ACRs/tamsACR_annotation.tsv.gz')

tamsACR_data <- readRDS('data/RData/tamsACR_main.rds') |>
  left_join(
    tamsACR_annotation,
    by = 'id'
  ) |>
  mutate(
    species = ordered(species, levels = species_short),
    condition = ordered(condition, levels = condition_order)
  )


### load eVal replicate data ###
eVal_reps <- readRDS('data/RData/eVal_reps.rds') |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )


### load main eVal data ###
eVal_data <- readRDS('data/RData/eVal_main.rds') |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )

### split eVal data by experiment ###
# TFBS shuffle
eVal_TFBS_shuffle <- eVal_data |>
  filter(str_detect(id, '_(WT)|(shuffle)')) |>
  separate_wider_delim(
    id,
    delim = '_',
    names = c('id', NA, 'TF', 'TFBS')
  ) |>
  select(condition, id, TF, TFBS, GC, enrichment) |>
  pivot_wider(
    names_from = TFBS,
    values_from = enrichment
  ) |>
  drop_na(WT, shuffle) |>
  mutate(
    diff = shuffle - WT,
    TF = paste0('TF_', TF)
  )

# TFBS insertion
eVal_TFBS_insertion <- eVal_data |>
  filter(str_detect(id, 'rnd-[dm][0-9]+_.*_.*_.*_.*') | str_detect(id, 'rnd-[dm][0-9]+_evo0_start_0$')) |>
  mutate(
    id = str_replace(id, 'evo0_start_0', 'WT_NA_NA_NA')
  ) |>
  separate_wider_delim(
    id,
    delim = '_',
    names = c('id', 'set', 'TF1', 'TF2', 'TF3')
  ) |>
  mutate(
    across(starts_with('TF'), ~ if_else(.x == 'NA', NA, .x)),
    across(starts_with('TF'), as.integer),
    across(starts_with('TF'), ~ if_else(is.na(.x), NA, paste0('TF_', .x))),
    n_TFs = 3 - (is.na(TF1) + is.na(TF2) + is.na(TF3)),
    GC_bin = if_else(str_detect(id, '-d'), 'low GC', 'high GC'),
    GC_bin = ordered(GC_bin, levels = c('low GC', 'high GC'))
  )

# evolution
eVal_evolution <- eVal_data |>
  filter(str_detect(id, 'evo')) |>
  separate_wider_delim(
    id,
    delim = '_',
    names = c('id', 'orientation', 'muts_per_round', 'objective', 'round'),
    too_few = 'align_end'
  ) |>
  unite(
    'id',
    id,
    orientation,
    na.rm = TRUE
  ) |>
  mutate(
    muts_per_round = str_replace(muts_per_round, 'evo', '')
  ) |>
  separate_longer_delim(
    c(round, objective, muts_per_round),
    delim = '|'
  ) |>
  mutate(
    across(c(round, muts_per_round), as.numeric)
  )

# validation sequences
eVal_validation <- eVal_data |>
  mutate(
    id = str_replace(id, '_evo0_start_0', ''),
    id = str_replace(id, '_TF_[0-9]+_WT', '_fwd')
  ) |>
  filter(id == 'noEnh' | str_detect(id, '_(fwd|rev)$')) |>
  separate_wider_delim(
    id,
    delim = '_',
    names = c('id', 'orientation'),
    too_few = 'align_start'
  )


### load dual-luciferase data
dualLuc_data <- readRDS('data/RData/DL_data.rds')


### load sequences ###
tamsACR_seqs <- Biostrings::readDNAStringSet('data/refseq/tamsACR_sequences.fa.gz')
validation_seqs <- Biostrings::readDNAStringSet('data/refseq/validation_sequences.fa')


### load transcription factor motifs ###
TF_motifs <- read_meme('data/extra_files/TF-clusters.meme')


### helper functions ###
# keep track of figure numbers
cur_fig <- 'Fig. 0a'
cur_ext_fig <- 'Extended Data Fig. 0a'
cur_supp_fig <- 'Supplementary Fig. 0a'
cur_supp_table <- 0
cur_supp_data <- 0

nextfigure <- function(cur_fig, increment = 1) {
  fig_parts <- str_match(cur_fig, '(.*Fig\\. )([0-9]+)(.)')[-1]
  fig_parts[2] <- as.integer(fig_parts[2]) + increment
  fig_parts[3] <- 'a'
  return(paste0(fig_parts, collapse = ''))
}

nextsubfig <- function(cur_fig, increment = 1) {
  fig_parts <- str_match(cur_fig, '(.*Fig\\. )([0-9]+)(.)')[-1]
  fig_parts[3] <- letters[which(letters == fig_parts[3]) + increment]
  return(paste0(fig_parts, collapse = ''))
}

nexttable <- function(cur_table, increment = 1) {
  return(cur_table + increment)
}

# reverse-complement a sequence
rev_comp <- function(x) {
  Biostrings::DNAStringSet(x) |>
    Biostrings::reverseComplement() |>
    as.character()
}


### define styles for exported excel tables ###
xlsx_bold <- createStyle(textDecoration = 'bold')
xlsx_wrap <- createStyle(wrapText = TRUE)
xlsx_2digit <- createStyle(numFmt = '0.00')
xlsx_mixed <- createStyle(numFmt = '[=1]0;[<0.001]0.00E+0;0.0000')
xlsx_seq_font <- createStyle(fontName = 'Courier New')
xlsx_right <- createStyle(halign = 'right')
xlsx_center <- createStyle(halign = 'center')
xlsx_green <- createStyle(bgFill = '#92D050')

supp_data_dir <- file.path('data', 'supplementary_data')


# Figure 1 ----------------------------------------------------------------
cur_fig <- nextfigure(cur_fig)
cur_fig <- nextsubfig(cur_fig, 2)# skip assay scheme & conditions

### PCA of tamsACR experiments ###
# reshape replicate data
tamsACR_matrix <- tamsACR_reps |>
  select(experiment, condition, id, orientation, enrichment) |>
  pivot_wider(
    names_from = c(condition, experiment),
    values_from = enrichment
  ) |>
  select(-id, -orientation) |>
  drop_na() |>
  mutate(
    across(everything(), ~ scale(.x)[,1])
  ) |>
  t()

# perform PCA
tamsACR_pca <- tamsACR_matrix |>
  prcomp(
    center = TRUE,
    scale. = TRUE
  )

# export data
tamsACR_pca$x |>
  as_tibble() |>
  bind_cols(
    sample = rownames(tamsACR_matrix)
  ) |>
  separate_wider_delim(
    sample,
    delim = '_',
    names = c('condition', 'experiment')
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  ) |>
  pggfplot(
    filename = 'tamsACR_PCA',
    x = PC2,
    y = PC1,
    axis_annotation = list(
      PC1 = summary(tamsACR_pca)$importance[2, 'PC1'] * 100,
      PC2 = summary(tamsACR_pca)$importance[2, 'PC2'] * 100
    )
  ) |>
  pggf_scatter(color = condition)


### enhancer strength by species (only light shown in Fig. 1; rest in Extended Data Fig. 1) ###
cur_fig <- nextsubfig(cur_fig)

tamsACR_data |>
  drop_na(species) |>
  pggfplot(
    filename = 'tamsACR_species',
    x = species,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin()


### correlation of enhancer strength in fwd and rev orientation (only light shown in Fig. 1; rest in Extended Data Fig. 1) ###
cur_fig <- nextsubfig(cur_fig)

tamsACR_data |>
  drop_na(orientation) |>
  select(condition, id, orientation, enrichment) |>
  pivot_wider(
    names_from = orientation,
    values_from = enrichment
  ) |>
  drop_na(fwd, rev) |>
  pggfplot(
    filename = 'tamsACR_cor_orientation',
    x = fwd,
    y = rev,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_hexbin() |>
  pggf_stats(!! pggf_stat_fns$correlation)


### fold-change in enhancer strength in fwd or rev orientation ###
cur_fig <- nextsubfig(cur_fig)

tamsACR_data |>
  group_by(condition, id) |>
  filter(n() == 2) |>
  summarise(
    diff = abs(enrichment[orientation == 'fwd'] - enrichment[orientation == 'rev'])
  ) |>
  ungroup() |>
  pggfplot(
    filename = 'tamsACR_orientation',
    x = condition,
    y = diff
  ) |>
  pggf_violin()


### enhancer strength in stable Arabidopsis lines ###
cur_fig <- nextsubfig(cur_fig)

# order enhancers by strength
dl_enh_order <- dualLuc_data |>
  filter(! enhancer %in% c('noEnh', '35S(GG)')) |>
  group_by(enhancer) |>
  summarise(
    l2ratio = median(l2ratio)
  ) |>
  arrange(l2ratio) |>
  pull(enhancer)

dl_enh_order <- c('noEnh', dl_enh_order)

# export enhancer strength data
dualLuc_data |>
  mutate(
    enhancer = ordered(enhancer, levels = dl_enh_order, labels = c('none', dl_enh_order[-1]))
  ) |>
  drop_na(enhancer) |>
  pggfplot(
    filename = 'DL_strength',
    x = enhancer,
    y = l2ratio
  ) |>
  pggf_boxplot()


### compare dual-luciferase and Plant STARR-seq data ###
cur_fig <- nextsubfig(cur_fig)

# combine data
dualLuc_cor <- dualLuc_data |>
  group_by(enhancer, orientation) |>
  summarise(
    CI = 0.5 * diff(t.test(l2ratio)$conf.int),
    l2ratio = mean(l2ratio),
    .groups = 'drop'
  ) |>
  inner_join(
    tamsACR_data |>
      filter(condition == 'light') |>
      select('enhancer' = id, orientation, enrichment),
    by = c('enhancer', 'orientation')
  ) |>
  mutate(
    enhancer = ordered(enhancer, levels = dl_enh_order, labels = c('none', dl_enh_order[-1]))
  ) |>
  drop_na(enhancer)

# export correlation plot
dualLuc_cor |>
  pggfplot(
    filename = 'DL_cor',
    x = enrichment,
    y = l2ratio
  ) |>
  pggf_trendline() |>
  pggf_scatter(
    color = enhancer,
    y_error = CI
  ) |>
  pggf_stats(!! pggf_stat_fns$correlation)


### correlation of enhancer strength and gene expression data (only light shown in Fig. 1; rest in Extended Data Fig. 3) ###
cur_fig <- nextsubfig(cur_fig)

# map gene names to the ones in the most recent annotation
Zm_gene_map <- read_tsv('data/extra_files/B73v4_to_B73v5.tsv', col_names = c('V4', 'V5')) |>
  deframe()

Sl_gene_map <- read_lines('data/extra_files/Sl_new_gene_IDs.txt')
names(Sl_gene_map) <- str_sub(Sl_gene_map, 1, 14)

tamsACR_genes <- tamsACR_data |>
  filter(str_sub(id, 1, 5) != 'Solyc' & ! str_detect(id, 'sh')) |>
  mutate(
    closest_TSS = str_replace(closest_TSS, '\\([+-]\\)', ''),
    closest_TSS = case_match(
      species,
      'Sl' ~ Sl_gene_map[closest_TSS],
      'Zm' ~ Zm_gene_map[closest_TSS],
      .default = closest_TSS
    )
  ) |>
  drop_na(closest_TSS) |>
  separate_longer_delim(
    closest_TSS,
    delim = ','
  )

# calculate mean enrichment from sequences in forward and reverse orientation
tamsACR_genes <- tamsACR_genes |>
  group_by(condition, species, id, GC, closest_TSS) |>
  summarise(
    enrichment = mean(enrichment)
  ) |>
  ungroup()

# download expression data
expression_data_IDs <- c(
  'At' = 'E-MTAB-7978',
  'Sl' = 'E-MTAB-4812',
  'Zm' = 'E-MTAB-4342',
  'Sb' = 'E-MTAB-5956'
)

for (species in names(expression_data_IDs)) {
  if (! file.exists(paste0('data/expression_data/', species, '_FPKM.tsv'))) {
    download.file(
      paste0('https://www.ebi.ac.uk/gxa/experiments-content/', expression_data_IDs[species], '/resources/ExperimentDownloadSupplier.RnaSeqBaseline/fpkms.tsv'),
      paste0('data/expression_data/', species, '_FPKM.tsv')
    )
  }
}

# function to load and preprocess expression data
load_expression_data <- function(species) {
  expression_data <- read_tsv(paste0('data/expression_data/', species,'_FPKM.tsv'), comment = '#') |>
    mutate(
      across(where(is.numeric), ~ if_else(.x <= 10, NA, log10(.x)))
    ) |>
    select(-`Gene Name`) |>
    rename('id' = `Gene ID`) |>
    pivot_longer(
      where(is.numeric),
      names_to = 'sample',
      values_to = 'expression'
    ) |>
    drop_na(expression)
  
  return(expression_data)
}

# load expression data
expression_data <- tibble(species = ordered(species_short, levels = species_short)) |>
  group_by(species) |>
  reframe(
    load_expression_data(species)
  )

# combine expression and Plant STARR-seq data and export plot
tamsACR_genes |>
  group_by(condition, species, closest_TSS) |>
  summarise(
    enrichment = max(enrichment)
  ) |>
  ungroup() |>
  inner_join(
    expression_data,
    join_by(closest_TSS == id, species == species),
    relationship = 'many-to-many'
  ) |>
  group_by(condition, species, sample) |>
  summarise(
    correlation = cor(enrichment, expression),
    n = n()
  ) |>
  ungroup() |>
  mutate(
    tissue = if_else(str_detect(sample, 'leaf'), 'leaves', 'other'),
  ) |>
  pggfplot(
    filename = 'tamsACR_expression',
    x = species,
    y = correlation,
    facet_col = condition
  ) |>
  pggf_boxplot() |>
  pggf_scatter(
    color = tissue
  )


# Supplementary Figure 1 --------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### replicate correlation plots ###
rep_combis <- tamsACR_reps |>
  distinct(condition, experiment) |>
  group_by(condition) |>
  reframe(
    exps = t(combn(experiment, 2))
  ) |>
  group_by(condition) |>
  mutate(
    title = paste(
      'replicate',
      as.integer(ordered(exps[, 1], levels = unique(c(exps)))),
      'vs.',
      as.integer(ordered(exps[, 2], levels = unique(c(exps))))
    )
  ) |>
  ungroup()

tamsACR_reps_wide <- tamsACR_reps |>
  select(experiment, condition, id, orientation, enrichment) |>
  pivot_wider(
    names_from = c(condition, experiment),
    values_from = enrichment
  )

for (cond in unique(rep_combis$condition)) {
  rep_combis |>
    filter(condition == cond) |>
    group_by(title) |>
    reframe(
      tamsACR_reps_wide |>
        select(
          'rep1' = paste(cond, exps[,1], sep = '_'),
          'rep2' = paste(cond, exps[,2], sep = '_')
        )
    ) |>
    drop_na() |>
    pggfplot(
      filename = paste0('tamsACR_cor_reps_', cond),
      x = rep1,
      y = rep2,
      facet_col = title,
      scales = 'square'
    ) |>
    pggf_hexbin() |>
    pggf_stats(!! pggf_stat_fns$correlation) |>
    print()
  
  cur_supp_fig <- nextsubfig(cur_supp_fig)
}


# Supplementary Figure 2 --------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### correlation of enhancer strength across conditions ###
condition_combis <- tamsACR_data |>
  distinct(condition) |>
  arrange(condition) |>
  pull() |>
  combn(2) |>
  t()

tamsACR_data_wide <- tamsACR_data |>
  select(condition, id, orientation, enrichment) |>
  pivot_wider(
    names_from = condition,
    values_from = enrichment,
    names_prefix = 'enrichment_'
  )

tibble(
  condition_1 = condition_combis[,1],
  condition_2 = condition_combis[,2]
) |>
  group_by(pick(everything())) |>
  reframe(
    tamsACR_data_wide |>
      select(
        enrichment_1 = paste0('enrichment_', condition_1),
        enrichment_2 = paste0('enrichment_', condition_2)
      )
  ) |>
  drop_na() |>
  pggfplot(
    filename = 'tamsACR_cor_conds',
    x = enrichment_1,
    y = enrichment_2,
    facet_col = condition_1,
    facet_row = condition_2,
    scales = 'square',
    axis_annotation = expand_grid(
      condY = unique(condition_combis[,2]),
      condX = unique(condition_combis[,1])
    )
  ) |>
  pggf_hexbin() |>
  pggf_stats(!! pggf_stat_fns$correlation)


# Extended Data Figure 1 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### extra conditions/species not shown in Fig. 1d,e ###
# nothing to do here


# Supplementary Figure 3 --------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### scheme of dual-luciferase assay; LaTeX only ###
# nothing to do here


# Supplementary Figure 4 --------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### ACR sequences vs. non-ACR controls ###
tamsACR_data |>
  drop_na(species) |>
  filter(str_sub(id, 1, 2) != 'So') |>
  mutate(
    ACR = ! str_detect(id, 'sh'),
    ACR = ordered(ACR, levels = c('TRUE', 'FALSE'), labels = c('ACR', 'non-ACR'))
  ) |>
  pggfplot(
    filename = 'tamsACR_non-ACR_ctrls',
    x = species,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin(
    half = ACR,
    signif = list(test = 'Wilcox', p_adjust = 'bonferroni', save_test_results = cur_supp_fig)
  )


# Extended Data Figure 2 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### enhancer strength by chromatin accessibility ###
tamsACR_data |>
  drop_na(cutcount) |>
  filter(str_sub(id, 1, 2) != 'So' & ! str_detect(id, 'sh')) |>
  group_by(species) |>
  mutate(
    accessibility = cut_number(cutcount, n = 3, labels = c('low', 'medium', 'high'))
  ) |>
  ungroup() |>
  pggfplot(
    filename = 'tamsACR_cutcount',
    x = accessibility,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin(
    signif = list(test = 'Tukey', save_test_results = cur_ext_fig, order_by = 'sample')
  )


### histone modifications ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

tamsACR_histones <- tamsACR_data |>
  filter(species == 'At' & ! str_detect(id, 'Bur-0')) |>
  mutate(
    histone_mod = case_when(
      str_detect(histones, 'ac') & str_detect(histones, 'me') ~ 'ac+me',
      str_detect(histones, 'ac') ~ 'ac',
      str_detect(histones, 'me') ~ 'me',
      .default = 'none'
    ),
    histone_mod = ordered(histone_mod, levels = c('none', 'ac', 'me', 'ac+me'))
  )

tamsACR_histones |>
  pggfplot(
    filename = 'tamsACR_histones',
    x = histone_mod,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin(
    signif = list(test = 'Tukey', save_test_results = cur_ext_fig, order_by = 'sample')
  )


### enhancer strength by genomic region ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

tamsACR_data |>
  drop_na(region) |>
  mutate(
    region = ordered(
      region,
      levels = c("intergenic", "5'UTR", "CDS", "intron", "3'UTR", "ncRNA")
    )
  ) |>
  pggfplot(
    filename = 'tamsACR_region',
    x = region,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin(
    signif = list(test = 'Tukey', save_test_results = cur_ext_fig, order_by = 'sample')
  )


### enhancer strength by distance to TSS ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

tamsACR_data |> 
  drop_na(TSS_dist) |>
  mutate(
    TSS_dist = cut(
      TSS_dist,
      breaks = c(min(TSS_dist), -5000, -1000, 0, 1000, 5000, max(TSS_dist)),
      labels = c('u5+', 'u5', 'u1', 'd1', 'd5', 'd5+'),
      include.lowest = TRUE
    )
  ) |>
  pggfplot(
    filename = 'tamsACR_TSS_dist',
    x = TSS_dist,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin(
    signif = list(test = 'Tukey', save_test_results = cur_ext_fig, order_by = 'sample')
  )


# Extended Data Figure 3 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### differentially expressed genes in light vs. dark in Arabidopsis ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

# load data
light_DEGs <- read_excel(
  path = 'data/extra_files/41598_2017_4524_MOESM5_ESM.xls',
  range = 'A3:J4262',
  .name_repair = ~ if_else(.x == '', paste0('V', seq_along(.x)), .x),
  na = c('', 'inf', 'inf(-)')
) |>
  select('gene' = `Locus ID`, 'light' = `Control-FPKM`, 'dark' = `Dark-FPKM`, 'expression_fc' = `log2(fold_change)`, 'category' = `up-or-down`) |>
  mutate(
    category = ordered(category, levels = c('down-regulated under darkness', 'up-regulated under darkness'), labels = c('dark_down', 'dark_up'))
  )

# combine DEGs and Plant STARR-seq data and export plot
tamsACR_genes |>
  filter(condition %in% c('light', 'dark') & species == 'At') |>
  pivot_wider(
    names_from = condition,
    values_from = enrichment
  ) |>
  drop_na(light, dark) |>
  mutate(
    enrichment_DvL = dark - light,
    enrichment_DvL = enrichment_DvL - median(enrichment_DvL),
    threshold = log2(1.5)
  ) |>
  inner_join(
    light_DEGs,
    join_by(closest_TSS == gene)
  ) |>
  mutate(
    specificity = case_when(
      enrichment_DvL > threshold ~ 'dark',
      enrichment_DvL < -threshold ~ 'light',
      .default = 'none'
    ),
    specificity = ordered(specificity, levels = c('light', 'none', 'dark'))
  ) |>
  count(category, specificity) |>
  group_by(category) |>
  mutate(
    percent = n / sum(n) * 100
  ) |>
  ungroup() |>
  pggfplot(
    filename = 'tamsACR_dark_light_specificity',
    x = category,
    y = percent
  ) |>
  pggf_bar(
    group = specificity
  )


### condition-specificity of enhancer strength vs. tissue-specificity of gene expression ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

# calculate condition-/tissue-specificity using the tau index
enhancer_tau <- tamsACR_genes |>
  group_by(condition) |>
  mutate(
    enrichment = enrichment - median(enrichment)
  ) |>
  ungroup() |>
  mutate(
    enrichment = (enrichment - min(enrichment)) / (max(enrichment) - min(enrichment))
  ) |>
  group_by(species, closest_TSS, id, GC) |>
  summarise(
    n = n(),
    tau = sum(1 - (enrichment / max(enrichment))) / (n() - 1)
  ) |>
  ungroup() |>
  drop_na(tau)

expression_tau <- expression_data |>
  group_by(species, id) |>
  summarise(
    n = n(),
    tau = sum(1 - (expression / max(expression))) / (n() - 1)
  ) |>
  ungroup() |>
  drop_na(tau)

# export plot data
enhancer_tau |>
  group_by(species, closest_TSS) |>
  summarise(
    tau = max(tau)
  ) |>
  ungroup() |>
  inner_join(
    expression_tau,
    join_by(closest_TSS == id, species),
    suffix = c('_enrichment', '_expression')
  ) |>
  group_by(species) |>
  mutate(
    tau_exp_bin = cut_number(tau_expression, n = 3, labels = c('low', 'medium', 'high'))
  ) |>
  ungroup() |>
  pggfplot(
    filename = 'tamsACR_specificity',
    x = tau_exp_bin,
    y = tau_enrichment
  ) |>
  pggf_violin(
    signif = list(test = 'Tukey', save_test_results = cur_ext_fig, order_by = 'sample')
  )
  

### GO terms ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

# creating GMT files with ACR IDs instead of gene names
for (species in species_short) {
  if (! file.exists(paste0('data/GO-terms/GOslim_ACR_', species, '.gmt'))) {
    read_delim(paste0('data/GO-terms/GOslim_', species, '.gmt'), delim = ';', col_names = 'data') |>
      separate_wider_delim(
        data,
        delim = '\t',
        names = c('GO_term', 'description', 'gene'),
        too_many = 'merge'
      ) |>
      separate_longer_delim(
        gene,
        delim = '\t'
      ) |>
      inner_join(
        distinct(tamsACR_genes, id, closest_TSS),
        join_by(gene == closest_TSS),
        relationship = 'many-to-many'
      ) |>
      group_by(GO_term, description) |>
      summarise(
        ACRs = paste(id, collapse = '\t')
      ) |>
      ungroup() |>
      write_tsv(
        paste0('data/GO-terms/GOslim_ACR_', species, '.gmt'),
        col_names = FALSE
      )
  }
}

# upload GO slim annotations
if (file.exists('data/GO-terms/GOslim_IDs_ACR.tsv')) {
  # load custom GMT IDs
  GOslim_IDs_ACR <- read_tsv('data/GO-terms/GOslim_IDs_ACR.tsv') |>
    deframe()
} else {
  # upload custom GMTs
  GOslim_IDs_ACR <- tibble(species = species_short) |>
    rowwise() |>
    mutate(
      ID = upload_GMT_file(paste0('data/GO-terms/GOslim_ACR_', species, '.gmt'))
    ) |>
    ungroup()
  
  write_tsv(GOslim_IDs_ACR, 'data/GO-terms/GOslim_IDs_ACR.tsv')
  
  GOslim_IDs_ACR <- GOslim_IDs_ACR |>
    deframe()
}

# find enriched GO terms
GOterms <- tamsACR_genes |>
  group_by(condition, species, id) |>
  summarise(
    enrichment = mean(enrichment)
  ) |>
  group_by(condition, species) |> 
  arrange(-enrichment, .by_group = TRUE) |>
  reframe(
    gost(
      query = id[seq_len(0.05 * n())],
      organism = GOslim_IDs_ACR[unique(species)],
      significant = FALSE,
      domain_scope = 'custom_annotated',
      custom_bg = id
    ) |>
      _$result
  ) |>
  ungroup()

# export plot data (only GO terms with significance in at least 3 species in any condition)
GOterms |>
  group_by(condition, term_name) |>
  mutate(
    n_signif = sum(significant),
    significant = as.numeric(significant)
  ) |>
  group_by(term_name) |>
  filter(any(n_signif >= 3)) |>
  ungroup() |>
  mutate(
    term_name = ordered(term_name, levels = c('photosynthesis', 'thylakoid', 'chloroplast', 'plastid', 'response to biotic stimulus', 'response to external stimulus', 'response to chemical', 'response to stress'))
  ) |>
  pggfplot(
    filename = 'tamsACR_GO-terms',
    x = species,
    y = term_name,
    facet_col = condition
  ) |>
  pggf_heatmap(
    score = significant
  )


# Extended Data Figure 4 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### enhancer strength by GC ###
GC_breaks <- tamsACR_data |>
  drop_na(GC) |>
  reframe(
    GC = quantile(GC * 100, seq(0, 1, .2)),
    GC = case_match(
      names(GC),
      '0%' ~ floor(GC),
      '100%' ~ ceiling(GC),
      .default = round(GC)
    )
  ) |>
  pull(GC)

tamsACR_data |>
  drop_na(GC) |>
  mutate(
    GC_bin = cut(
      x = GC * 100,
      breaks = GC_breaks,
      include.lowest = TRUE,
      ordered_result = TRUE
    )
  ) |>
  pggfplot(
    filename = 'tamsACR_GC',
    x = GC_bin,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin(
    signif = list(test = 'Tukey', save_test_results = cur_ext_fig, order_by = 'sample')
  )


### GC distribution of test sequences ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

tamsACR_data |>
  distinct(id, species, GC) |>
  drop_na(species, GC) |>
  mutate(
    GC = GC * 100
  ) |>
  pggfplot(
    'GC_distribution',
    x = GC,
    y = species
  ) |>
  pggf_ridgeline()


### compare average GC content of genomes and ACR sequences ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

tamsACR_data |>
  filter(str_sub(id, 1, 5) != 'Solyc' & ! str_detect(id, 'sh')) |>
  distinct(id, species, GC) |>
  drop_na(species, GC) |>
  group_by(species) |>
  summarise(
    GC = mean(GC)
  ) |>
  inner_join(
    read_tsv('data/annotation/genome_GC.tsv') |>
      mutate(
        species = species_to_short[species]
      ),
    by = 'species',
    suffix = c('_ACR', '_genome')
  ) |>
  pivot_longer(
    starts_with('GC_'),
    names_to = c('.value', 'region'),
    names_pattern = '(GC)_(.*)'
  ) |>
  mutate(
    species = ordered(species, species_short),
    GC = GC * 100
  ) |>
  group_by(species) |>
  mutate(
    diff = GC[region == 'ACR'] - GC,
  ) |>
  ungroup() |>
  unite(
    'group',
    species,
    region,
    sep = ':',
    remove = FALSE
  ) |>
  pggfplot(
    filename = 'GC_ACR_genome',
    x = GC,
    y = species
  ) |>
  pggf_scatter(
    color = group
  ) |>
  pggf_quiver(
    u = diff
  )


# Figure 2 ----------------------------------------------------------------
cur_fig <- nextfigure(cur_fig)

### effect of transcription factor binding sites ###
# load transcription factor binding site counts
TF_hits <- read_tsv('data/extra_files/TF_hits.tsv.gz')

# GC normalize enhancer strength
GC_coeffs <- tamsACR_data |>
  drop_na(GC) |>
  group_by(condition) |>
  reframe(
    lm(enrichment ~ GC) |>
      coefficients() |>
      enframe()
  ) |>
  mutate(
    name = if_else(name != 'GC', 'intercept', name)
  ) |>
  unite(
    'name',
    condition,
    name
  ) |>
  deframe()

tamsACR_data_GCnorm <- tamsACR_data |>
  drop_na(GC) |>
  mutate(
    enrichment = enrichment - (GC_coeffs[paste0(condition, '_GC')] * GC + GC_coeffs[paste0(condition, '_intercept')])
  )

# combine enhancer strength and TF hits
tamsACR_TFs <- tamsACR_data_GCnorm |>
  inner_join(
    TF_hits,
    by = 'id'
  )

# determine effect of presence/absence of a single transcription factor binding site
TF_effects <- tamsACR_TFs |>
  group_by(condition) |>
  summarise(
    across(
      starts_with('TF_'),
      ~ mean(enrichment[.x == 1]) - mean(enrichment[.x == 0])
    )
  ) |>
  ungroup() |>
  pivot_longer(
    starts_with('TF_'),
    names_to = 'TF',
    values_to = 'effect'
  )

# sort results by effect in "light" condition
TF_order <- TF_effects |>
  filter(condition == 'light') |>
  arrange(effect) |>
  pull(TF)

TF_effects <- TF_effects |>
  mutate(
    TF = ordered(TF, levels = TF_order)
  )

# export plot data
TF_effects |>
  pggfplot(
    filename = 'TF_effects',
    x = TF,
    y = condition,
  ) |>
  pggf_heatmap(
    score = effect,
    dendrogram = 'y'
  )
  
# export TF family names
TF_motifs |>
  sapply(function(x) x['altname']) |>
  str_replace('_', '-') |>
  as_tibble() |>
  mutate(
    id = paste0('TF_', seq_len(n())),
    family = paste0(str_replace_all(value, '/', '\\\\slash '), ' (', seq_len(n()), ')')
  ) |>
  select(id, family) |>
  write_tsv(paste0(pggf_defaults$out_dir, 'TF_families.tsv'))

# export positions of TFs for validation
validation_TFs <- read_tsv('data/extra_files/TFs_for_validation.txt') |>
  mutate(
    TF = paste0('TF_', TF)
  )
  
tibble('pos'= which(TF_order %in% validation_TFs$TF)) |>
  write_tsv(paste0(pggf_defaults$out_dir, 'TF_position.tsv'))


### compare shuffle effect with TF effect in tamsACR library ###
cur_fig <- nextsubfig(cur_fig, increment = 2) # skip scheme of TFBS shuffling experiment

eVal_TFBS_shuffle |>
  group_by(condition, TF) |>
  summarise(
    effect = -mean(diff),
  ) |>
  ungroup() |>
  inner_join(
    TF_effects,
    by = c('condition', 'TF'),
    suffix = c('_shuffle', '_tamsACR')
  ) |>
  pggfplot(
    filename = 'TFBS_shuffle_tamsACR',
    x = effect_shuffle,
    y = effect_tamsACR,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_trendline() |>
  pggf_scatter() |>
  pggf_stats(!! pggf_stat_fns$correlation)


### nucleotide content of random sequences ###
cur_fig <- nextsubfig(cur_fig)

eVal_TFBS_insertion |>
  distinct(GC_bin, tmp = str_sub(id, 5, 5)) |>
  group_by(GC_bin) |>
  reframe(
    validation_seqs[str_detect(names(validation_seqs), paste0('rnd-', tmp, '[0-9]+_evo0_start_0$'))] |>
      Biostrings::alphabetFrequency(as.prob = TRUE, baseOnly = TRUE) |>
      colMeans() |>
      enframe(
        name = 'base',
        value = 'freq'
      )
  ) |>
  filter(base != 'other') |>
  mutate(
    angle = freq * 360
  ) |>
  group_by(GC_bin) |>
  summarise(
    tmp = paste(base, cumsum(angle), sep = '/', collapse = ',')
  ) |>
  ungroup() |>
  pivot_wider(
    names_from = GC_bin,
    values_from = tmp
  ) |>
  write_tsv(paste0(pggf_defaults$out_dir, 'synEnh_nucFreq.tsv'))


### export selected TFBSs ###
synEnh_TFs <- eVal_TFBS_insertion |>
  distinct(TF1) |>
  drop_na() |>
  pull()

synEnh_TFs[order(as.numeric(str_replace(synEnh_TFs, 'TF_', '')))] |>
  matrix(ncol = 2) |>
  as_tibble(.name_repair = ~ paste0('col', seq_along(.x))) |>
  write_tsv(paste0(pggf_defaults$out_dir, 'synEnh_TFs.tsv'))


### strength of synthetic enhancers created by inserting TFBSs into random sequences ###
cur_fig <- nextsubfig(cur_fig)

eVal_TFBS_insertion |>
  mutate(
    n_TFs = ordered(n_TFs)
  ) |>
  pggfplot(
    filename = 'synEnh_strength',
    x = n_TFs,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_boxplot(
    group = GC_bin,
    signif = list(test = 'Tukey', save_test_results = cur_fig, order_by = 'sample')
  ) |>
  pggf_hline(
    fn = median(y[x == 1]),
    group = GC_bin
  )


# Supplementary Figure 5 --------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### examples of TF effects ###
# select validation TFs
TF_examples <- tamsACR_TFs |>
  select(condition, orientation, id, any_of(validation_TFs$TF), enrichment) |>
  pivot_longer(
    any_of(validation_TFs$TF),
    names_to = 'TF',
    values_to = 'hits'
  ) |>
  filter(hits <= 1) |>
  mutate(
    hits = if_else(hits == 0, 'absent', 'present'),
    hits = ordered(hits, levels = c('absent', 'present')),
    condition_name = condition
  ) |>
  mutate(
    TF = ordered(TF, levels = TF_order),
    TF = droplevels(TF)
  )

# plot data
TF_examples |>
  pggfplot(
    filename = 'TF_effect_examples',
    x = TF,
    y = enrichment,
    facet_row = condition,
    scales = 'free_y',
    preview_plot = FALSE
  ) |>
  pggf_violin(
    half = hits,
    signif = list(test = 'Wilcox', p_adjust = 'bonferroni', save_test_results = cur_supp_fig),
    extra_cols = condition_name
  )


# Supplementary Figure 6 --------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### correlation of TF effects across conditions ###
condition_combis <- TF_effects |>
  distinct(condition) |>
  arrange(condition) |>
  pull() |>
  combn(2) |>
  t()

TF_effects_wide <- TF_effects |>
  pivot_wider(
    names_from = condition,
    values_from = effect
  )

tibble(
  condition_1 = condition_combis[,1],
  condition_2 = condition_combis[,2]
) |>
  group_by(pick(everything())) |>
  reframe(
    TF_effects_wide |>
      select(
        effect_1 = all_of(condition_1),
        effect_2 = all_of(condition_2)
      )
  ) |>
  drop_na() |>
  pggfplot(
    filename = 'TF_effect_cor_conds',
    x = effect_1,
    y = effect_2,
    facet_col = condition_1,
    facet_row = condition_2,
    scales = 'square',
    axis_annotation = expand_grid(
      condY = unique(condition_combis[,2]),
      condX = unique(condition_combis[,1])
    )
  ) |>
  pggf_scatter() |>
  pggf_stats(!! pggf_stat_fns$correlation)


# Extended Data Figure 5 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### explain sequence categories ###
# since we use GC-normalized enhancer strength, we can't categorize sequences based on the noEnh control
# instead, we use the non-ACR controls as a baseline
tamsACR_categories <- tamsACR_TFs |>
  group_by(condition) |>
  mutate(
    mean = mean(enrichment[str_detect(id, 'sh')]),
    sd = sd(enrichment[str_detect(id, 'sh')]),
    category = case_when(
      enrichment < mean - 1 * sd ~ 'repressive',
      enrichment > mean + 3 * sd ~ 'strong_enh',
      enrichment > mean + 2 * sd ~ 'medium_enh',
      enrichment > mean + 1 * sd ~ 'weak_enh',
      .default = 'inactive'
    ),
    category = ordered(category, levels = rev(c('repressive', 'inactive', 'weak_enh', 'medium_enh', 'strong_enh')))
  ) |>
  ungroup()

# export plot data
tamsACR_categories |>
  filter(condition == 'light') |>
  mutate(
    sample = 'all',
    pos = 5
  ) |>
  bind_rows(
    tamsACR_categories |>
      filter(condition == 'light' & str_detect(id, 'sh')) |>
      mutate(
        sample = 'non-ACR',
        pos = 1
      )
  ) |>
  pggfplot(
    filename = 'category_def',
    x = pos,
    y = enrichment
  ) |>
  pggf_violin() |>
  pggf_hline(
    list(c(
      'mean' = mean[1],
      'mean-1' = max(y[category == 'repressive']),
      'mean+1' = min(y[category == 'weak_enh']),
      'mean+2' = min(y[category == 'medium_enh']),
      'mean+3' = min(y[category == 'strong_enh']),
      'max' = max(y),
      'min' = min(y)
    ))
  )


### number of sequences in each category ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

tamsACR_categories |>
  count(condition, category) |>
  mutate(
    category_label = category
  ) |>
  pggfplot(
    filename = 'category_numbers',
    x = n,
    y = category,
    facet_col = condition
  ) |>
  pggf_bar(
    orientation = 'y',
    group = category_label,
    stack = TRUE
  )


### enrichment/depletion of transcription factor binding sites in sequences grouped by their strength ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

# for each category, we calculate the average counts for all TFBSs and normalize them to the "inactive" sequences
categories_TFs <- tamsACR_categories |>
  group_by(condition, category) |>
  summarise(
    across(starts_with('TF_'), mean)
  ) |>
  group_by(condition) |>
  mutate(
    across(starts_with('TF_'), ~ log2(.x / .x[category == 'inactive']))
  ) |>
  ungroup() |>
  filter(category != 'inactive') |>
  pivot_longer(
    starts_with('TF_'),
    names_to = 'TF',
    values_to = 'norm_count'
  ) |>
  inner_join(
    TF_effects,
    by = c('condition', 'TF')
  ) |>
  filter(is.finite(norm_count))

# export plot data
categories_TFs |>
  pggfplot(
    filename = 'categories_TFs',
    x = effect,
    y = norm_count,
    facet_col = condition,
    facet_row = category
  ) |>
  pggf_scatter() |>
  pggf_trendline() |>
  pggf_stats(!! pggf_stat_fns$trendline)


# Supplementary Figure 7 --------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### replicate correlation of eVal library ###
rep_combis <- eVal_reps |>
  distinct(condition, experiment) |>
  group_by(condition) |>
  reframe(
    exps = t(combn(experiment, 2))
  ) |>
  group_by(condition) |>
  mutate(
    title = paste(
      'replicate',
      as.integer(ordered(exps[, 1], levels = unique(c(exps)))),
      'vs.',
      as.integer(ordered(exps[, 2], levels = unique(c(exps))))
    )
  ) |>
  ungroup()

eVal_reps_wide <- eVal_reps |>
  select(experiment, condition, id, enrichment) |>
  pivot_wider(
    names_from = c(condition, experiment),
    values_from = enrichment
  )

for (cond in unique(rep_combis$condition)) {
  rep_combis |>
    filter(condition == cond) |>
    group_by(title) |>
    reframe(
      eVal_reps_wide |>
        select(
          'rep1' = paste(cond, exps[,1], sep = '_'),
          'rep2' = paste(cond, exps[,2], sep = '_')
        )
    ) |>
    drop_na() |>
    pggfplot(
      filename = paste0('eVal_cor_reps_', cond),
      x = rep1,
      y = rep2,
      facet_col = title,
      scales = 'square'
    ) |>
    pggf_hexbin() |>
    pggf_stats(!! pggf_stat_fns$correlation) |>
    print()
  
  cur_supp_fig <- nextsubfig(cur_supp_fig)
}


### PCA of eVal experiments ###
# reshape replicate data
eVal_matrix <- eVal_reps |>
  select(experiment, condition, id, enrichment) |>
  pivot_wider(
    names_from = c(condition, experiment),
    values_from = enrichment
  ) |>
  select(-id) |>
  drop_na() |>
  mutate(
    across(everything(), ~ scale(.x)[,1])
  ) |>
  t()

# perform PCA
eVal_pca <- eVal_matrix |>
  prcomp(
    center = TRUE,
    scale. = TRUE
  )

# export plot data
eVal_pca$x |>
  as_tibble() |>
  bind_cols(
    sample = rownames(eVal_matrix)
  ) |>
  separate_wider_delim(
    sample,
    delim = '_',
    names = c('condition', 'experiment')
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  ) |>
  pggfplot(
    filename = 'eVal_PCA',
    x = PC2,
    y = PC1,
    axis_annotation = list(
      PC1 = summary(eVal_pca)$importance[2, 'PC1'] * 100,
      PC2 = summary(eVal_pca)$importance[2, 'PC2'] * 100
    )
  ) |>
  pggf_scatter(color = condition)


# Supplementary Figure 8 ----------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### correlation between tamsACR and eVal libraries ###
tamsACR_data |>
  select(condition, id, orientation, enrichment) |>
  inner_join(
    eVal_validation |>
      select(condition, id, orientation, enrichment),
    by = c('condition', 'id', 'orientation'),
    suffix = c('_tamsACR', '_eVal')
  ) |>
  pggfplot(
    filename = 'tamsACR_eVal_cor',
    x = enrichment_tamsACR,
    y = enrichment_eVal,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_hexbin() |>
  pggf_stats(!! pggf_stat_fns$correlation)


# Supplementary Figure 9 ---------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### effect of TFBS shuffling ###
# sort by effect in light
shuffle_order <- eVal_TFBS_shuffle |>
  filter(condition == 'light') |>
  group_by(TF) |>
  summarise(
    diff = median(diff)
  ) |>
  arrange(diff) |>
  pull(TF)

# export plot data
eVal_TFBS_shuffle |>
  mutate(
    TF = ordered(TF, levels = shuffle_order)
  ) |>
  pggfplot(
    filename = 'TFBS_shuffle_effect',
    x = TF,
    y = diff,
    facet_col = condition
  ) |>
  pggf_boxplot(
    signif = list(test = 'Wilcox', comparisons = 'one-sample tests', p_adjust = 'bonferroni', save_test_results = cur_supp_fig)
  )


# Extended Data Figure 6 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### effect of TFBS insertion ###
cur_ext_fig <- nextsubfig(cur_ext_fig)# skip experiment scheme

# select synthetic enhancers with single transcription factor binding site
eVal_TFBS_insertion_single <- eVal_TFBS_insertion |>
  filter(n_TFs %in% c(0, 1)) |>
  group_by(condition, id) |>
  filter(any(set == 'WT')) |>
  mutate(
    diff = enrichment - enrichment[set == 'WT']
  ) |>
  ungroup() |>
  filter(set != 'WT') |>
  select(condition, id, GC_bin, set, starts_with('TF'), diff) |>
  mutate(
    pos = case_when(
      ! is.na(TF1) ~ 'TF1',
      ! is.na(TF2) ~ 'TF2',
      ! is.na(TF3) ~ 'TF3'
    )
  ) |>
  unite(
    'TF',
    starts_with('TF'),
    na.rm = TRUE
  )

# sort by effect in light
insertion_order <- eVal_TFBS_insertion_single |>
  filter(condition == 'light') |>
  group_by(TF) |>
  summarise(
    diff = median(diff)
  ) |>
  arrange(-diff) |>
  pull(TF)

# export plot data
eVal_TFBS_insertion_single |>
  mutate(
    TF = ordered(TF, levels = insertion_order)
  ) |>
  pggfplot(
    filename = 'TFBS_insertion_effect',
    x = TF,
    y = diff,
    facet_col = condition
  ) |>
  pggf_boxplot(
    signif = list(test = 'Wilcox', comparisons = 'one-sample tests', p_adjust = 'bonferroni', save_test_results = cur_ext_fig)
  )


### compare effect of TFBS insertion with TF effect in tamsACR library ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

eVal_TFBS_insertion_single |>
  group_by(condition, TF) |>
  summarise(
    effect = mean(diff),
  ) |>
  ungroup() |>
  inner_join(
    TF_effects,
    by = c('condition', 'TF'),
    suffix = c('_insertion', '_tamsACR')
  ) |>
  pggfplot(
    filename = 'TFBS_insertion_tamsACR',
    x = effect_insertion,
    y = effect_tamsACR,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_trendline() |>
  pggf_scatter() |>
  pggf_stats(!! pggf_stat_fns$correlation)


# Extended Data Figure 7 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### predict synEnh strength with different inputs ###
cur_ext_fig <- nextsubfig(cur_ext_fig)# skip explanation of methods for TF effects

# get strength of random sequences
base_strength <- eVal_TFBS_insertion |>
  filter(set == 'WT') |>
  unite(
    'tmp',
    condition,
    id
  ) |>
  select(tmp, enrichment) |>
  deframe()

# get effect of TFBS insertion in a given position and random sequence
single_insertion_effects <- eVal_TFBS_insertion_single |>
  group_by(condition, id, pos, TF) |>
  summarise(
    effect = mean(diff)
  ) |>
  ungroup() |>
  unite(
    'tmp',
    condition,
    id,
    pos,
    TF
  ) |>
  select(tmp, effect) |>
  deframe()

# combine data and predict enrichment
synEnh_prediction <- eVal_TFBS_insertion |>
  filter(n_TFs >= 2) |>
  mutate(
    across(starts_with('TF'), ~ single_insertion_effects[paste(condition, id, cur_column(), .x, sep = '_')], .names = 'effect_{.col}'),
    across(starts_with('effect_TF'), ~ replace_na(.x, 0)),
    base_strength = base_strength[paste(condition, id, sep = '_')],
    prediction = base_strength + effect_TF1 + effect_TF2 + effect_TF3
  ) |>
  drop_na(prediction) |>
  select(condition, id, enrichment, base_strength, starts_with(c('TF', 'effect_TF')), prediction)

# randomize TF effects
set.seed(928)

synEnh_prediction_random <- synEnh_prediction |>
  group_by(condition) |>
  mutate(
    across(starts_with('effect_TF'), ~ sample(.x, n())),
    prediction = base_strength + effect_TF1 + effect_TF2 + effect_TF3,
    input = 'random'
  ) |>
  ungroup()

# use only sequence id but not position
single_insertion_effects_noPos <- eVal_TFBS_insertion_single |>
  group_by(condition, id, TF) |>
  summarise(
    effect = mean(diff)
  ) |>
  ungroup() |>
  unite(
    'tmp',
    condition,
    id,
    TF
  ) |>
  select(tmp, effect) |>
  deframe()

synEnh_prediction_noPos <- synEnh_prediction |>
  mutate(
    across(starts_with('TF'), ~ single_insertion_effects_noPos[paste(condition, id, .x, sep = '_')], .names = 'effect_{.col}'),
    across(starts_with('effect_TF'), ~ replace_na(.x, 0)),
    prediction = base_strength + effect_TF1 + effect_TF2 + effect_TF3
  ) |>
  drop_na(prediction)

# use only position but not sequence id
single_insertion_effects_noID <- eVal_TFBS_insertion_single |>
  group_by(condition, pos, TF) |>
  summarise(
    effect = mean(diff)
  ) |>
  ungroup() |>
  unite(
    'tmp',
    condition,
    pos,
    TF
  ) |>
  select(tmp, effect) |>
  deframe()

synEnh_prediction_noID <- synEnh_prediction |>
  mutate(
    across(starts_with('TF'), ~ single_insertion_effects_noID[paste(condition, cur_column(), .x, sep = '_')], .names = 'effect_{.col}'),
    across(starts_with('effect_TF'), ~ replace_na(.x, 0)),
    prediction = base_strength + effect_TF1 + effect_TF2 + effect_TF3
  ) |>
  drop_na(prediction)

# use mean across sequences and position
single_insertion_effects_mean <- eVal_TFBS_insertion_single |>
  group_by(condition, TF) |>
  summarise(
    effect = mean(diff)
  ) |>
  ungroup() |>
  unite(
    'tmp',
    condition,
    TF
  ) |>
  select(tmp, effect) |>
  deframe()

synEnh_prediction_mean <- synEnh_prediction |>
  mutate(
    across(starts_with('TF'), ~ single_insertion_effects_mean[paste(condition, .x, sep = '_')], .names = 'effect_{.col}'),
    across(starts_with('effect_TF'), ~ replace_na(.x, 0)),
    prediction = base_strength + effect_TF1 + effect_TF2 + effect_TF3
  ) |>
  drop_na(prediction)

# use TF effect from tamsACR library
single_insertion_effects_tamsACR <- TF_effects |>
  unite(
    'tmp',
    condition,
    TF
  ) |>
  select(tmp, effect) |>
  deframe()

synEnh_prediction_tamsACR <- synEnh_prediction |>
  mutate(
    across(starts_with('TF'), ~ single_insertion_effects_tamsACR[paste(condition, .x, sep = '_')], .names = 'effect_{.col}'),
    across(starts_with('effect_TF'), ~ replace_na(.x, 0)),
    prediction = base_strength + effect_TF1 + effect_TF2 + effect_TF3
  ) |>
  drop_na(prediction)

# combine all and export plot data
bind_rows(
  'id+pos' = synEnh_prediction,
  'no_pos' = synEnh_prediction_noPos,
  'no_id' = synEnh_prediction_noID,
  'mean' = synEnh_prediction_mean,
  'tamsACR' = synEnh_prediction_tamsACR,
  'random' = synEnh_prediction_random,
  .id = 'input'
) |>
  mutate(
    input = ordered(input, levels = c('id+pos', 'no_pos', 'no_id', 'mean', 'tamsACR', 'random'))
  ) |>
  group_by(condition, input) |>
  summarise(
    rsquare = cor(enrichment, prediction)^2
  ) |>
  ungroup() |>
  pggfplot(
    filename = 'synEnh_pred_accuracy',
    x = input,
    y = rsquare,
    facet_col = condition
  ) |>
  pggf_bar()


### correlation plot for synEnh predictions ###
cur_fig <- nextsubfig(cur_fig)

synEnh_prediction |>
  pggfplot(
    filename = 'synEnh_prediction',
    x = prediction,
    y = enrichment,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_hexbin() |>
  pggf_stats(!! pggf_stat_fns$correlation)


# Supplementary Figure 10 -------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### variability of TF effect sizes in synthetic enhancers ###
cur_supp_fig <- nextsubfig(cur_supp_fig) # skip scheme

# calculate mean and standard deviation of TF effects in synthetic enhancers across sequences or positions
synEnh_effects_noPos <- eVal_TFBS_insertion_single |>
  group_by(condition, TF, id) |>
  summarise(
    mean = mean(diff),
    sd = sd(diff)
  ) |>
  ungroup()

synEnh_effects_noId <- eVal_TFBS_insertion_single |>
  group_by(condition, TF, pos) |>
  summarise(
    mean = mean(diff),
    sd = sd(diff)
  ) |>
  ungroup()

# combine and export plot data
bind_rows(
  'acrossPos' = synEnh_effects_noPos,
  'acrossSeq' = synEnh_effects_noId,
  .id = 'sample'
) |>
  drop_na(sd) |>
  pggfplot(
    filename = 'synEnh_TF_sd',
    x = sample,
    y = sd,
    facet_col = condition
  ) |>
  pggf_violin(
    signif = list(test = 'Wilcox', comparisons = list(c(1, 2)), p_adjust = 'none', save_test_results = cur_supp_fig)
  )


# Extended Data Figure 8 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### strength of synthetic enhancers with species-/condition-specific TFs ###
cur_ext_fig <- nextsubfig(cur_ext_fig)# skip table with selected TFBS

# get species-/condition-specificity of TFs
TF_specificities <- TF_effects |>
  filter(TF %in% synEnh_TFs) |>
  pivot_wider(
    names_from = condition,
    values_from = effect
  ) |>
  mutate(
    TF = as.character(TF),
    `tobacco-specific` = (light + dark + warm + cold) / 4 - maize,
    `maize-specific` = -`tobacco-specific`,
    `light-specific` = light - dark,
    `dark-specific` = -`light-specific`,
    `warm-specific` = warm - cold,
    `cold-specific` = -`warm-specific`
  ) |>
  select(TF, ends_with('-specific')) |>
  pivot_longer(
    ends_with('-specific'),
    names_to = 'specificity',
    values_to = 'effect'
  ) |>
  mutate(
    specificity = ordered(specificity, levels = c('tobacco-specific', 'maize-specific', 'light-specific', 'dark-specific', 'warm-specific', 'cold-specific'))
  ) |>
  group_by(specificity) |>
  slice_max(effect, n = 1) |>
  ungroup()

# export selected TFs
TF_specificities |>
  arrange(specificity) |>
  write_tsv(paste0(pggf_defaults$out_dir, 'synEnh_TF_specificities.tsv'))

TF_specificities <- TF_specificities |>
  select(TF, specificity) |>
  deframe()

# select putative species-/condition-specific synthetic enhancers
synEnh_specific <- eVal_TFBS_insertion |>
  filter(set == 'WT' | (TF1 == TF2 & TF1 == TF3)) |>
  mutate(
    set = if_else(set == 'WT', 'noTFBS', TF_specificities[TF1])
  ) |>
  drop_na(set)

# export plot data
synEnh_specific |>
  mutate(
    set = ordered(set, levels = c('noTFBS', paste0(c('tobacco', 'maize', 'light', 'dark', 'warm', 'cold'), '-specific'))),
    target = case_when(
      condition == 'light' & set %in% c('tobacco-specific', 'light-specific') ~ 'high',
      condition == 'light' & set %in% c('maize-specific', 'dark-specific') ~ 'low',
      condition == 'dark' & set %in% c('tobacco-specific', 'dark-specific') ~ 'high',
      condition == 'dark' & set %in% c('maize-specific', 'light-specific') ~ 'low',
      condition == 'warm' & set %in% c('tobacco-specific', 'warm-specific') ~ 'high',
      condition == 'warm' & set %in% c('maize-specific', 'cold-specific') ~ 'low',
      condition == 'cold' & set %in% c('tobacco-specific', 'cold-specific') ~ 'high',
      condition == 'cold' & set %in% c('maize-specific', 'warm-specific') ~ 'low',
      condition == 'maize' & set == 'maize-specific' ~ 'high',
      condition == 'maize' & set == 'tobacco-specific' ~ 'low',
      .default = 'any'
    )
  ) |>
  pggfplot(
    filename = 'synEnh_specific_strength',
    x = set,
    y = enrichment,
    facet_col = condition,
    # facet_row = GC_bin
  ) |>
  pggf_boxplot(
    extra_cols = target,
    signif = list(test = 'Tukey', save_test_results = cur_ext_fig, order_by = 'sample')
  ) |>
  pggf_hline(
    fn = median(y[x == 1])
  )


### species-/condition-specificity of synthetic enhancers ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

# species-specificity
species_spec <- synEnh_specific |>
  filter(set %in% c('noTFBS', 'tobacco-specific', 'maize-specific')) |>
  group_by(set, GC_bin, id, pick(starts_with('TF'))) |>
  filter(n() > 1 & any(condition == 'maize')) |>
  summarise(
    specificity = mean(enrichment[condition != 'maize']) - enrichment[condition == 'maize'],
    category = 'species-specificity'
  ) |>
  ungroup()

# light-specificity
light_spec <- synEnh_specific |>
  filter(set %in% c('noTFBS', 'light-specific', 'dark-specific') & condition %in% c('light', 'dark')) |>
  group_by(set, GC_bin, id, pick(starts_with('TF'))) |>
  filter(n() == 2) |>
  summarise(
    specificity = enrichment[condition == 'light'] - enrichment[condition == 'dark'],
    category = 'light-specificity'
  ) |>
  ungroup()

# temperature-specificity
temp_spec <- synEnh_specific |>
  filter(set %in% c('noTFBS', 'warm-specific', 'cold-specific') & condition %in% c('warm', 'cold')) |>
  group_by(set, GC_bin, id, pick(starts_with('TF'))) |>
  filter(n() == 2) |>
  summarise(
    specificity = enrichment[condition == 'warm'] - enrichment[condition == 'cold'],
    category = 'temperature-specificity'
  ) |>
  ungroup()

# combine specificities, normalize data, and export plot
bind_rows(species_spec, light_spec, temp_spec) |>
  group_by(category) |>
  mutate(
    specificity = specificity - median(specificity[set == 'noTFBS']),
    set = ordered(set, levels = c('tobacco-specific', 'light-specific', 'warm-specific', 'noTFBS', 'cold-specific', 'dark-specific', 'maize-specific'))
  ) |>
  ungroup() |>
  mutate(
    category = ordered(category, levels = paste0(c('species', 'light', 'temperature'), '-specificity')),
    target_species = if_else(set == 'noTFBS', 'none', str_replace(set, '-specific', ''))
  ) |>
  pggfplot(
    filename = 'synEnh_specificity',
    x = set,
    y = specificity,
    facet_col = category,
    scales = 'free_x'
  ) |>
  pggf_boxplot(
    extra_cols = target_species,
    signif = list(test = 'Tukey', save_test_results = cur_ext_fig, order_by = 'sample')
  )


# Figure 3 ----------------------------------------------------------------
cur_fig <- nextfigure(cur_fig)

### compare plantGREP predictions to measured data ###
cur_fig <- nextsubfig(cur_fig)

# function to make column names of prediction dataframe more verbose
fix_col_names <- function(col_names) {
  model_outputs <- c('light', 'dark', 'warm', 'cold', 'maize')
  
  name_parts <- col_names |>
    str_match('(.*)_(.*)')
  
  new_names <- paste(
    if_else(name_parts[,2] == 'target', 'enrichment', 'prediction'),
    model_outputs[as.numeric(name_parts[,3]) + 1],
    sep = '_'
  )
  
  return(new_names)
}

# load plantGREP predictions
plantGREP_predictions <- read_tsv('data/modelling_data/plantGREP_predictions.tsv') |>
  rename_with(fix_col_names) |>
  pivot_longer(
    everything(),
    names_to = c('.value', 'condition'),
    names_pattern = '(.*)_(.*)'
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )

# export plot data
plantGREP_predictions |>
  pggfplot(
    filename = 'plantGREP_predictions',
    x = enrichment,
    y = prediction,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_hexbin() |>
  pggf_stats(!! pggf_stat_fns$correlation)


### example DeepLIFT output ###
cur_fig <- nextsubfig(cur_fig)

# load and export data for plot of example sequence (At-12806_fwd)
deeplift_example <- read_tsv('data/modelling_data/DeepLIFT/DeepLIFT_example.tsv')

deeplift_example |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  ) |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'contrib'
  ) |>
  pggfplot(
    filename = 'DeepLIFT_example',
    x = pos,
    y = contrib,
    facet_row = condition,
    preview_plot = FALSE
  ) |>
  pggf_logo(
    letters = base
  )

# export plot for TCP binding motif
convert_type(TF_motifs[[22]], 'ICM')['motif'] |>
  t() |>
  as_tibble() |>
  mutate(
    pos = seq_len(n())
  ) |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'IC'
  ) |>
  pggfplot(
    filename = 'TF_22',
    x = pos,
    y = IC,
    preview_plot = FALSE
  ) |>
  pggf_logo(
    letters = base
  )

# export plot for MYB binding motif
convert_type(TF_motifs[[13]], 'ICM')['motif'] |>
  t() |>
  as_tibble() |>
  mutate(
    pos = n() - seq_len(n()) + 1
  ) |>
  rename_with(rev_comp, c(A, C, G, T)) |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'IC'
  ) |>
  pggfplot(
    filename = 'TF_13_RC',
    x = pos,
    y = IC,
    preview_plot = FALSE
  ) |>
  pggf_logo(
    letters = base
  )

# export plot for bZIP binding motif
convert_type(TF_motifs[[7]], 'ICM')['motif'] |>
  t() |>
  as_tibble() |>
  mutate(
    pos = n() - seq_len(n()) + 1
  ) |>
  rename_with(rev_comp, c(A, C, G, T)) |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'IC'
  ) |>
  pggfplot(
    filename = 'TF_7_RC',
    x = pos,
    y = IC,
    preview_plot = FALSE
  ) |>
  pggf_logo(
    letters = base
  )


### TF-modisco analysis ###
cur_fig <- nextsubfig(cur_fig)

# load data
modisco_data <- read_tsv('data/modelling_data/TF-MoDISco/TF-MoDISco_patterns.tsv.gz') |>
  separate_wider_regex(
    pattern,
    patterns = c('pattern_type' = '.*', '_pattern_', 'pattern_id' = '[0-9]+')
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order),
    pattern_type = ordered(pattern_type, levels = c('pos', 'neg')),
    pattern_id = as.integer(pattern_id)
  ) |>
  arrange(condition, pattern_type, pattern_id)

# trim TF-MoDISco patterns
trim_pattern <- function(pattern, threshold = 0.3, extend = 0) {
  trimmed <- pattern |>
    arrange(pos) |>
    mutate(
      .tp.score = abs(cwm_A) + abs(cwm_C) + abs(cwm_G) + abs(cwm_T),
      .tp.keep = .tp.score >= threshold * max(.tp.score),
      .tp.extended = between(pos, min(which(.tp.keep)) - extend, max(which(.tp.keep)) + extend)
    ) |>
    filter(.tp.extended) |>
    select(-starts_with(fixed('.tp.')))
  
  return(trimmed)
}

modisco_trimmed <- modisco_data |>
  nest_by(condition, pattern_type, pattern_id) |>
  reframe(
    trim_pattern(data, extend = 1)
  ) |>
  pivot_longer(
    starts_with(c('ppm_', 'cwm_')),
    names_to = c('motif_type', '.value'),
    names_pattern = '(.*)_(.)'
  )

# add reverse-complemented patterns
modisco_trimmed <- modisco_trimmed |>
  group_by(condition, pattern_type, pattern_id) |>
  mutate(
    orientation = 'rev',
    pos = max(pos) - pos + 1
  ) |>
  arrange(pos, .by_group = TRUE) |>
  ungroup() |>
  rename_with(rev_comp, c(A, C, G, T)) |>
  bind_rows(modisco_trimmed) |>
  mutate(
    orientation = replace_na(orientation, 'fwd'),
    orientation = ordered(orientation, levels = c('fwd', 'rev'))
  ) |>
  arrange(condition, pattern_type, pattern_id, orientation, motif_type, pos)

# generate motifs from modisco ppms
modisco_motifs <- modisco_trimmed |>
  filter(motif_type == 'ppm') |>
  nest_by(condition, pattern_type, pattern_id, orientation) |>
  summarise(
    motif = list(
      data[c('A', 'C', 'G', 'T')] |>
        as.matrix() |>
        t() |>
        create_motif(name = paste(condition, pattern_type, pattern_id, orientation, sep = '_'))
    )
  ) |>
  ungroup() |>
  pull(motif)

# find matches to known transcription factor binding motifs
modisco_TF_match <- compare_motifs(c(TF_motifs, modisco_motifs), seq_along(TF_motifs), tryRC = FALSE) |>
  as_tibble() |>
  filter(! str_detect(target, 'TF-cluster')) |>
  separate_wider_delim(
    target,
    delim = '_',
    names = c('condition', 'pattern_type', 'pattern_id', 'orientation')
  ) |>
  select(condition, pattern_type, pattern_id, orientation, 'TF' = subject, 'TF_id' = subject.i, score, 'p_value' = Pval) |>
  group_by(condition, pattern_type, pattern_id, TF) |>
  slice_min(p_value, n = 1) |>
  ungroup() |>
  mutate(
    condition = ordered(condition, levels = condition_order),
    pattern_type = ordered(pattern_type, levels = c('pos', 'neg')),
    pattern_id = as.integer(pattern_id),
    orientation = ordered(orientation, levels = c('fwd', 'rev')),
    TF = str_replace(TF, '-cluster', '')
  ) |>
  arrange(condition, pattern_type, pattern_id, p_value)

# select transcription factors to show in this figure
TFs_to_show <- c(21, 15, 13, 6, 36, 2, 40, 35)

# select corresponding modisco patterns
modisco_patterns_to_show <- modisco_TF_match |>
  filter(TF_id %in% TFs_to_show) |>
  group_by(condition, pattern_type, pattern_id) |>
  slice_min(p_value) |>
  group_by(condition, TF_id) |>
  slice_max(score) |>
  ungroup() |>
  unite(
    'modisco_motif',
    condition,
    pattern_type,
    pattern_id,
    orientation
  ) |>
  mutate(
    TF = ordered(TF, levels = paste0('TF_', TFs_to_show))
  )

# align motifs and patterns and extract alignment info
modisco_patterns_aligned <- modisco_patterns_to_show |>
  group_by(TF) |>
  reframe(
    view_motifs(c(TF_motifs[[unique(TF_id)]], filter_motifs(modisco_motifs, name = modisco_motif)), tryRC = FALSE, return.raw = TRUE) |>
      sapply(function(x) min(which(colSums(x) > 0))) |>
      enframe(name = 'motif', value = 'offset'),
    offset = offset - offset[str_detect(motif, 'TF-cluster_')],
    TFBS_len = nchar(TF_motifs[[unique(TF_id)]]['consensus'])
  )

# add the modisco patterns and trim them to only the part matching the TF motif
modisco_patterns_aligned <- modisco_patterns_aligned |>
  filter(! str_detect(motif, 'TF-cluster_')) |>
  separate_wider_delim(
    motif,
    delim = '_',
    names = c('condition', 'pattern_type', 'pattern_id', 'orientation')
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order),
    pattern_id = as.integer(pattern_id)
  ) |>
  left_join(
    modisco_trimmed,
    by = c('condition', 'pattern_type', 'pattern_id', 'orientation')
  ) |>
  filter(motif_type == 'cwm') |>
  group_by(condition, pattern_type, pattern_id, orientation) |>
  mutate(
    pos = pos - min(pos) + 1 + offset
  ) |>
  ungroup() |>
  filter(between(pos, 1, TFBS_len))
  
# export patterns for plotting
modisco_patterns_aligned |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'contrib'
  ) |>
  right_join(
    modisco_patterns_aligned |>
      expand(condition, TF),
    by = c('condition', 'TF')
  ) |>
  pggfplot(
    filename = 'modisco_motifs_examples',
    x = pos,
    y = contrib,
    facet_col = TF,
    facet_row = condition,
    scales = 'free_x',
    space = 'free_x',
    preview_plot = FALSE
  ) |>
  pggf_logo(
    letters = base
  )

# export TF motifs for plotting
modisco_patterns_to_show |>
  distinct(TF, TF_id) |>
  group_by(TF) |>
  reframe(
    convert_type(TF_motifs[[TF_id]], 'ICM')['motif'] |>
      t() |>
      as_tibble() |>
      mutate(
        pos = seq_len(n())
      )
  ) |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'IC'
  ) |>
  pggfplot(
    filename = 'modisco_TFs_examples',
    x = pos,
    y = IC,
    facet_col = TF,
    scales = 'free_x',
    space = 'free_x',
    axis_annotation = modisco_patterns_to_show |> distinct(TF) |> arrange(TF) |> select('TF name' = TF),
    preview_plot = FALSE
  ) |>
  pggf_logo(
    letters = base
  )


# Extended Data Figure 9 --------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### linear model predictions for TF motif matches ###
cur_ext_fig <- nextsubfig(cur_ext_fig)# skip model scheme

# load linear model predictions for TF hits (generated by file `code/modelling/linear_models.py`)
TF_hits_predictions <- read_tsv('data/modelling_data/lm_TF_hits_predictions.tsv.gz') |>
  pivot_longer(
    -c(id, orientation),
    names_to = c('.value', 'condition'),
    names_pattern = '(.*)_(.*)'
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )

# export plot data
TF_hits_predictions |>
  pggfplot(
    filename = 'lm_TF_hits_predictions',
    x = enrichment,
    y = prediction,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_hexbin() |>
  pggf_stats(!! pggf_stat_fns$correlation)
  

### linear model predictions for TF motif scores ###
cur_ext_fig <- nextsubfig(cur_ext_fig, increment = 2)# skip model scheme

# load TF scores
TF_motif_scores <- read_tsv('data/extra_files/TF_motif_scores.tsv.gz')

# load linear model predictions for TF scores (generated by file `code/modelling/linear_models.py`)
TF_scores_predictions <- read_tsv('data/modelling_data/lm_TF_scores_predictions.tsv.gz') |>
  pivot_longer(
    -c(id, orientation),
    names_to = c('.value', 'condition'),
    names_pattern = '(.*)_(.*)'
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )

# export plot data
TF_scores_predictions |>
  pggfplot(
    filename = 'lm_TF_scores_predictions',
    x = enrichment,
    y = prediction,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_hexbin() |>
  pggf_stats(!! pggf_stat_fns$correlation)


### compare model coefficients ###
cur_ext_fig <- nextsubfig(cur_ext_fig)

# combine coefficients
TF_coefficients <- bind_rows(
  'hits' = read_tsv('data/modelling_data/lm_TF_hits_coefficients.tsv.gz'),
  'scores' = read_tsv('data/modelling_data/lm_TF_scores_coefficients.tsv.gz'),
  .id = 'metric'
) |>
  pivot_longer(
    starts_with('coefficient'),
    names_to = c('.value', 'condition'),
    names_pattern = '(.*)_(.*)'
  ) |>
  pivot_wider(
    names_from = metric,
    values_from = coefficient,
    names_prefix = 'coefficient_'
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )

# add TF labels
TF_coefficients_labelled <- TF_coefficients |>
  mutate(
    TF = as.numeric(str_replace(TF, 'TF_', '')),
    type = case_when(
      condition == 'light' & coefficient_scores > 0.25 ~ 'high',
      condition == 'light' & coefficient_scores < -0.2 ~ 'low',
      condition == 'light' & coefficient_hits <  -0.35 ~ 'low',
      condition == 'dark' & coefficient_scores > 0.25 ~ 'high',
      condition == 'dark' & coefficient_scores < -0.22 ~ 'low',
      condition == 'dark' & coefficient_hits < -0.4 ~ 'low',
      condition == 'warm' & coefficient_scores > 0.25 ~ 'high',
      condition == 'warm' & coefficient_scores < -0.15 ~ 'low',
      condition == 'warm' & coefficient_hits < -0.32 ~ 'low',
      condition == 'cold' & coefficient_scores > 0.2 ~ 'high',
      condition == 'cold' & coefficient_scores < -0.15 ~ 'low',
      condition == 'cold' & coefficient_hits < -0.3 ~ 'low',
      condition == 'maize' & coefficient_scores > 0.15 ~ 'high',
      condition == 'maize' & coefficient_scores < -0.08 ~ 'low'
    ),
    pin = case_when(
      is.na(type) ~ NA,
      TF == 15 & condition == 'warm' ~ 45,
      TF == 22 & condition == 'maize' ~ -45,
      TF == 2 & condition == 'dark' ~ 135,
      TF == 2 & condition == 'warm' ~ 0,
      TF == 1 ~ 45,
      TF == 58 ~ 90,
      TF == 33 ~ 180,
      TF == 71 ~ 0,
      TF %in% c(7, 35, 42) ~ -45,
      TF %in% c(15, 22) ~ 135,
      type == 'high' ~ 90,
      type == 'low' ~ -90
    )
  )

# export plot data
TF_coefficients_labelled |>
  pggfplot(
    filename = 'lm_TF_coefficients',
    x = coefficient_hits,
    y = coefficient_scores,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_scatter(
    extra_cols = c(TF, pin)
  ) |>
  pggf_stats(!! pggf_stat_fns$correlation)


# Supplementary Figure 11 -------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### correlation between TF effects and TF model coefficients ###
TF_coeff_cor_effect <- TF_coefficients |>
  inner_join(TF_effects, by = c('condition', 'TF')) |>
  pivot_longer(
    starts_with('coefficient_'),
    names_to = c('.value', 'model'),
    names_pattern = '(.*)_(.*)'
  ) |>
  mutate(
    model = ordered(model, levels = c('hits', 'scores'))
  )

TF_coeff_cor_effect |>
  pggfplot(
    filename = 'lm_TF_cor_effect',
    x = effect,
    y = coefficient,
    facet_col = condition,
    facet_row = model,
    scales = 'square'
  ) |>
  pggf_scatter() |>
  pggf_trendline() |>
  pggf_stats(!! pggf_stat_fns$correlation)


# Extended Data Figure 10 -------------------------------------------------
cur_ext_fig <- nextfigure(cur_ext_fig)

### linear model predictions for kmer counts ###
cur_ext_fig <- nextsubfig(cur_ext_fig)# skip model scheme

# load linear model predictions for kmers (generated by file `code/modelling/linear_models.py`)
kmer_predictions <- read_tsv('data/modelling_data/lm_kmers_predictions.tsv.gz') |>
  pivot_longer(
    -c(id, orientation),
    names_to = c('.value', 'condition'),
    names_pattern = '(.*)_(.*)'
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )

# export plot data
kmer_predictions |>
  pggfplot(
    filename = 'lm_kmers_predictions',
    x = enrichment,
    y = prediction,
    facet_col = condition,
    scales = 'square'
  ) |>
  pggf_hexbin() |>
  pggf_stats(!! pggf_stat_fns$correlation)


# Supplementary Figure 12-21 ----------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### TF-MoDISco motifs and TF matches ###
## select modisco motifs to show
# since the orientation from modisco is essentially arbitrary, we use the orientation best matching a TF motif (or fwd if there is no match)
# for each modisco motif we will show only the best matching TF motif
modisco_best_TF_match <- modisco_TF_match |>
  group_by(condition, pattern_type, pattern_id) |>
  slice_min(p_value) |>
  ungroup() |>
  right_join(
    modisco_trimmed |>
      filter(motif_type == 'cwm'),
    by = c('condition', 'pattern_type', 'pattern_id', 'orientation')
  ) |>
  group_by(condition, pattern_type, pattern_id) |>
  filter(! is.na(TF) | (all(is.na(TF)) & orientation == 'fwd')) |>
  ungroup()

# function to determine y tick positions for the modiso motif plots
get_ytick <- function(pos, score, group = NULL) {
  ytick <- tibble(pos = pos, score = score, group = group) |>
    group_by(pick(any_of('group'), pos)) |>
    summarise(
      max = max(sum(score[score > 0]), abs(sum(score[score < 0]))),
      .groups = 'drop'
    ) |>
    group_by(pick(any_of('group'))) |>
    summarise(
      max = max(max),
      ytick = case_when(
        max >= 0.1 ~ paste0(seq(-round(max, 1) - 0.1, round(max, 1) + 0.1, 0.1), collapse = ','),
        max >= 0.05 ~ paste0(seq(-round(max, 1) - 0.1, round(max, 1) + 0.1, 0.05), collapse = ','),
        .default = paste(-floor(max * 100) / 100, 0, floor(max * 100) / 100, sep = ',')
      ),
      .groups = 'drop'
    ) |>
    arrange(pick(any_of('group'))) |>
    select(ytick)
  
  return(ytick)
}

# export modisco motifs for plotting
modisco_best_TF_match |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'contrib'
  ) |>
  nest_by(condition, pattern_type) |>
  summarise(
    pggfplot = list(
      pggfplot(
        data = data,
        filename = paste('modisco_motifs', pattern_type, condition, sep = '_'),
        x = pos,
        y = contrib,
        facet_row = pattern_id,
        scales = 'free_y',
        axis_annotation = get_ytick(data$pos, data$contrib, data$pattern_id),
        preview_plot = FALSE
      ) |>
        pggf_logo(
          letters = base,
          align = 'center'
        )
    )
  ) |>
  pull(pggfplot)

# export TF motifs for plotting
modisco_best_TF_match |>
  distinct(condition, pattern_type, pattern_id, TF, TF_id) |>
  group_by(condition, pattern_type, pattern_id, TF) |>
  reframe(
    motif = ifelse(
      is.na(TF_id),
      list(NA),
      list(
        convert_type(TF_motifs[[TF_id]], 'ICM')['motif'] |>
          t() |>
          as_tibble() |>
          mutate(
            pos = seq_len(n())
          )
      )
    )
  ) |>
  unnest(motif) |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'IC'
  ) |>
  nest_by(condition, pattern_type) |>
  reframe(
    pggfplot = list(
      pggfplot(
        data = data,
        filename = paste('modisco_TFs', pattern_type, condition, sep = '_'),
        x = pos,
        y = IC,
        facet_row = pattern_id,
        axis_annotation = data |> distinct(pattern_id, TF) |> arrange(pattern_id) |> select('TF name' = TF),
        preview_plot = FALSE
      ) |>
        pggf_logo(
          letters = base,
          align = 'center'
        )
    )
  ) |>
  pull(pggfplot)

# adjust number of Supplementary Figures
cur_supp_fig <- nextfigure(cur_supp_fig, increment = nrow(distinct(modisco_best_TF_match, condition, pattern_type)) - 1)


# Supplementary Figure 22 -------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### predictions from plantGREP ensemble ###
# load predictions
ensemble_predictions <- read_tsv('data/modelling_data/plantGREP_ensemble/ensemble_predictions.tsv.gz') |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  )

# randomly select models to combine (pick 50 different combinations, making sure that each model appears once in the first position)
set.seed(928)

model_list <- lapply(seq_len(50), function(x) c(x, sample(seq_len(50)[-x], 49)))

# calculate rsquares for using n models (0 < n < 51) based on the model_list; do this for every model_list element
ensemble_rsquares <- tibble(n_models = seq_len(50)) |>
  group_by(n_models) |>
  reframe(
    seed = seq_len(50),
    models = sapply(model_list, function(x) list(x[1:n_models]))
  ) |>
  group_by(n_models, seed) |>
  reframe(
    ensemble_predictions |>
      nest_by(condition) |>
      reframe(
        enrichment = data$enrichment,
        prediction = rowMeans(data[paste0('prediction_', models[[1]])])
      ) |>
      group_by(condition) |>
      summarise(
        rsquare = cor(enrichment, prediction)^2
      )
  )

# export raw plot data
ensemble_rsquares |>
  filter(n_models %in% c(1, 5, 10, 20, 50)) |>
  mutate(
    n_models = ordered(n_models)
  ) |>
  pggfplot(
    filename = 'plantGREP_ensemble_raw',
    x = n_models,
    y = rsquare,
    facet_col = condition
  ) |>
  pggf_boxplot(
    signif = list(test = 'Tukey', save_test_results = cur_supp_fig, order_by = 'sample')
  )

# export summary plot data
cur_supp_fig <- nextsubfig(cur_supp_fig)

ensemble_rsquares |>
  group_by(condition, n_models) |>
  summarise(
    sd = sd(rsquare),
    rsquare = mean(rsquare)
  ) |>
  ungroup() |>
  pggfplot(
    filename = 'plantGREP_ensemble_summary',
    x = n_models,
    y = rsquare
  ) |>
  pggf_line(
    group = condition,
    error = sd
  )


# Figure 4 ----------------------------------------------------------------
cur_fig <- nextfigure(cur_fig)

### in silico evolution by number of mutations ###
cur_fig <- nextsubfig(cur_fig, increment = 1)# skip in silico evolution scheme

eVal_evolution |>
  filter(objective %in% c('start', 'constitutive') & muts_per_round %in% c(0, 1)) |>
  mutate(
    round = ordered(round)
  ) |>
  pggfplot(
    filename = 'evolution_rounds',
    x = round,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin(
    signif = list(test = 'Tukey', save_test_results = cur_fig, order_by = 'sample')
  )


### in silico evolution specificity ###
cur_fig <- nextsubfig(cur_fig)

# species-specificity
species_spec <- eVal_evolution |>
  filter(objective %in% c('start', 'tobacco-specific', 'maize-specific') & muts_per_round <= 1) |>
  group_by(id, objective, round) |>
  filter(n() > 1 & any(condition == 'maize')) |>
  summarise(
    specificity = mean(enrichment[condition != 'maize']) - enrichment[condition == 'maize'],
    facet = 'species-specificity'
  ) |>
  ungroup() |>
  mutate(
    round = if_else(objective == 'tobacco-specific', -round, round)
  )

# light-specificity
light_spec <- eVal_evolution |>
  filter(objective %in% c('start', 'light-specific', 'dark-specific') & condition %in% c('light', 'dark') & muts_per_round <= 1) |>
  group_by(id, objective, round) |>
  filter(n() == 2) |>
  summarise(
    specificity = enrichment[condition == 'light'] - enrichment[condition == 'dark'],
    facet = 'light-specificity'
  ) |>
  ungroup() |>
  mutate(
    round = if_else(objective == 'light-specific', -round, round)
  )

# temperature-specificity
temp_spec <- eVal_evolution |>
  filter(objective %in% c('start', 'warm-specific', 'cold-specific') & condition %in% c('warm', 'cold') & muts_per_round <= 1) |>
  group_by(id, objective, round) |>
  filter(n() == 2) |>
  summarise(
    specificity = enrichment[condition == 'warm'] - enrichment[condition == 'cold'],
    facet = 'temperature-specificity'
  ) |>
  ungroup() |>
  mutate(
    round = if_else(objective == 'warm-specific', -round, round)
  )

# combine specificities, normalize data, and export plot
bind_rows(species_spec, light_spec, temp_spec) |>
  group_by(facet) |>
  mutate(
    specificity = specificity - median(specificity[round == 0])
  ) |>
  ungroup() |>
  mutate(
    color = str_replace(objective, '-specific', ''),
    fade = abs(round) == 6,
    facet = ordered(facet, levels = paste0(c('species', 'light', 'temperature'), '-specificity')),
    round = ordered(round),
  ) |>
  pggfplot(
    filename = 'evolution_specificity',
    x = round,
    y = specificity,
    facet_col = facet
  ) |>
  pggf_violin(
    extra_cols = c(color, fade),
    signif = list(test = 'Tukey', save_test_results = cur_fig, order_by = 'sample')
  )


# Supplementary Figure 23 -------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### in silico evolution with multiple mutations per round ###
eVal_evolution |>
  filter(objective %in% c('start', 'constitutive')) |>
  mutate(
    mutations = round * muts_per_round,
    mutations = if_else(mutations == 0, 0, (1.1 * mutations / 2) - 1),
    mutations = round(mutations / 3, 5),
    muts_per_round = if_else(objective == 'start', 2, muts_per_round),
    muts_per_round = paste0('mpr', muts_per_round)
  ) |>
  pggfplot(
    filename = 'evolution_muts_per_round',
    x = mutations,
    y = enrichment,
    facet_col = condition,
    preview_plot = FALSE
  ) |>
  pggf_violin(
    group = muts_per_round,
    signif = list(test = 'Wilcox', p_adjust = 'bonferroni', comparisons = 'sample groups', save_test_results = cur_supp_fig)
  )


# Supplementary Figure 24 -------------------------------------------------
cur_supp_fig <- nextfigure(cur_supp_fig)

### in silico evolution with different objectives ###
objective_order <- c('no evolution', 'constitutive', 'tobacco-specific', 'maize-specific', 'light-specific', 'dark-specific', 'warm-specific', 'cold-specific')

eVal_evolution |>
  mutate(
    objective = if_else(objective == 'start', 'no evolution', objective),
    objective = ordered(objective, levels = objective_order)
  ) |>
  drop_na(objective) |>
  filter(objective == 'no evolution' | (muts_per_round == 1 & round == 12)) |>
  mutate(
    target = case_when(
      condition == 'light' & objective %in% c('constitutive', 'tobacco-specific', 'light-specific') ~ 'high',
      condition == 'light' & objective %in% c('constitutive', 'maize-specific', 'dark-specific') ~ 'low',
      condition == 'dark' & objective %in% c('constitutive', 'tobacco-specific', 'dark-specific') ~ 'high',
      condition == 'dark' & objective %in% c('constitutive', 'maize-specific', 'light-specific') ~ 'low',
      condition == 'warm' & objective %in% c('constitutive', 'tobacco-specific', 'warm-specific') ~ 'high',
      condition == 'warm' & objective %in% c('constitutive', 'maize-specific', 'cold-specific') ~ 'low',
      condition == 'cold' & objective %in% c('constitutive', 'tobacco-specific', 'cold-specific') ~ 'high',
      condition == 'cold' & objective %in% c('constitutive', 'maize-specific', 'warm-specific') ~ 'low',
      condition == 'maize' & objective %in% c('constitutive', 'maize-specific') ~ 'high',
      condition == 'maize' & objective %in% c('constitutive', 'tobacco-specific') ~ 'low',
      .default = 'any'
    )
  ) |>
  pggfplot(
    filename = 'evolution_objectives',
    x = objective,
    y = enrichment,
    facet_col = condition
  ) |>
  pggf_violin(
    extra_cols = target,
    signif = list(test = 'Tukey', save_test_results = cur_supp_fig, order_by = 'sample')
  )


# Figure 5 ----------------------------------------------------------------
cur_fig <- nextfigure(cur_fig)

### SlCLV3 results from Wang et al., 2021, Nature Plants, https://doi.org/10.1038/s41477-021-00898-x ###
# load SlCLV3 data
SlCLV3_alleles <- read_tsv('data/extra_files/SlCLV3_alleles.tsv') |>
  mutate(
    id = ordered(id, levels = rev(c('gRNAs', 'WT', 'CLV3p13', 'CLV3p15', 'CLV3p24', 'CLV3p26')))
  )

SlCLV3_locules <- read_excel('data/extra_files/41477_2021_898_MOESM3_ESM.xlsx', range = 'A4:E3651') |>
  mutate(
    id = ordered(Genotype, levels = rev(c('gRNAs', 'WT', 'SlCLV3pro-13', 'SlCLV3pro-15', 'SlCLV3pro-24', 'SlCLV3pro-26')), labels = rev(c('gRNAs', 'WT', 'CLV3p13', 'CLV3p15', 'CLV3p24', 'CLV3p26')))
  ) |>
  drop_na(id)

# export plot data
SlCLV3_alleles |>
  mutate(
    start = if_else(start > end, start + 0.5, start - 0.5),
    end = if_else(end > start, end + 0.5, end - 0.5),
    length = end - start,
    across(c(start, end, length), ~ .x/1000),
    panel = 'genotype'
  ) |>
  bind_rows(
    tibble(
      panel = 'phenotype',
      id = ordered('WT', levels = levels(SlCLV3_alleles$id)),
      feature = 'ignore',
      start = c(0, 1.05 * max(SlCLV3_locules$LoculeNumber)),
      length = 0
    )
  ) |>
  mutate(
    panel = ordered(panel, c('genotype', 'phenotype'))
  ) |>
  pggfplot(
    filename = 'SlCLV3_alleles',
    x = start,
    y = id,
    facet_col = panel,
    scales = 'free_x',
    axis_annotation = list(xaxis = c('genotype', 'phenotype'))
  ) |>
  pggf_quiver(
    u = length,
    split = feature
  )

SlCLV3_locules |>
  pggfplot(
    filename = 'SlCLV3_locules',
    x = LoculeNumber,
    y = id
  ) |>
  pggf_boxplot(
    orientation = 'y'
  )


### enhancer strength of fragments upstream of SlCLV3 ###
cur_fig <- nextsubfig(cur_fig)

SlCLV3_ATG <- 58038742

SlCLV3_STARR <- tamsACR_data |>
  filter(str_detect(id, 'Solyc11g071380') & condition == 'dark') |>
  mutate(
    across(c(start, end), ~ SlCLV3_ATG - .x - 0.5),
    tmp = if_else(orientation == 'rev', end, start),
    end = if_else(orientation == 'rev', start, end),
    start = tmp,
    length = end - start,
    across(c(start, end, length), ~ .x/1000)
  )

SlCLV3_ids <- SlCLV3_STARR |>
  distinct(id) |>
  pull()

SlCLV3_seqs <- setNames(
  c(tamsACR_seqs[SlCLV3_ids], Biostrings::reverseComplement(tamsACR_seqs[SlCLV3_ids])),
  paste(SlCLV3_ids, rep(c('fwd', 'rev'), each = length(SlCLV3_ids)), sep = '_')
)

Biostrings::writeXStringSet(SlCLV3_seqs, 'data/modelling_data/SlCLV3_sequences.fa')

## predict enhancer strength of SlCLV3 sequences using plantGREP
# `python code/plantGREPcli/plantGREPcli.py predict -i data/modelling_data/SlCLV3_sequences.fa -o data/modelling_data/SlCLV3_predictions.tsv`

# load predictions
SlCLV3_predictions <- read_tsv('data/modelling_data/SlCLV3_predictions.tsv') |>
  select(name, prediction_dark) |>
  pivot_longer(
    starts_with('prediction_'),
    names_to = c('.value', 'condition'),
    names_pattern = '(.*)_(.*)'
  ) |>
  mutate(
    condition = ordered(condition, levels = condition_order)
  ) |>
  separate_wider_delim(
    name,
    delim = '_',
    names = c('id', 'orientation')
  )

# get coordinates of strongest SlCLV3 promoter fragment
SlCLV3_ROI <- SlCLV3_STARR |>
  slice_max(enrichment, n = 1) |>
  select('start' = end, 'end' = start)

# combine experimental and predicted data and export
SlCLV3_STARR |>
  inner_join(
    SlCLV3_predictions,
    by = join_by(condition, orientation, id)
  ) |>
  select(condition, id, orientation, start, length, 'experiment' = enrichment, prediction) |>
  pivot_longer(
    c(experiment, prediction),
    names_to = 'data_source',
    values_to = 'enrichment'
  ) |>
  pggfplot(
    filename = 'SlCLV3_PlantSTARRseq',
    x = start,
    y = enrichment,
    axis_annotation = list(ROIstart = SlCLV3_ROI$start, ROIend = SlCLV3_ROI$end)
  ) |>
  pggf_quiver(
    u = length,
    split = data_source
  )


### ZmCLE7 results from Liu et al., 2021, Nature Plants, https://doi.org/10.1038/s41477-021-00858-5 ###
cur_fig <- nextsubfig(cur_fig)

# load ZmCLE7 data
ZmCLE7_alleles <- read_tsv('data/extra_files/ZmCLE7_alleles.tsv') |>
  mutate(
    id = ordered(id, levels = rev(c('gRNAs', 'WT', 'CLE7p1', 'CLE7p4', 'CLE7p6')))
  )

ZmCLE7_expression <- read_tsv('data/extra_files/ZmCLE7_expression.tsv') |>
  mutate(
    id = ordered(id, levels = rev(c('gRNAs', 'WT', 'CLE7p1', 'CLE7p4', 'CLE7p6')))
  )

ZmCLE7_expression_points <- ZmCLE7_expression |>
  filter(str_detect(feature, 'point'))
ZmCLE7_expression_bar <- ZmCLE7_expression |>
  filter(! str_detect(feature, 'point'))

# export plot data
ZmCLE7_alleles |>
  mutate(
    start = if_else(start > end, start + 0.5, start - 0.5),
    end = if_else(end > start, end + 0.5, end - 0.5),
    length = end - start,
    across(c(start, end, length), ~ .x/1000),
    panel = 'genotype'
  ) |>
  bind_rows(
    tibble(
      panel = 'phenotype',
      id = ordered('WT', levels = levels(ZmCLE7_alleles$id)),
      feature = 'ignore',
      start = c(0, 1.05 * max(ZmCLE7_expression_points$expression)),
      length = 0
    )
  ) |>
  mutate(
    panel = ordered(panel, c('genotype', 'phenotype'))
  ) |>
  pggfplot(
    filename = 'ZmCLE7_alleles',
    x = start,
    y = id,
    facet_col = panel,
    scales = 'free_x',
    axis_annotation = list(xaxis = c('genotype', 'phenotype'))
  ) |>
  pggf_quiver(
    u = length,
    split = feature
  )

ZmCLE7_expression_points |>
  pggfplot(
    filename = 'ZmCLE7_expression',
    x = expression,
    y = id
  ) |>
  pggf_scatter()

ZmCLE7_expression_bar |>
  pivot_wider(
    names_from = feature,
    values_from = expression
  ) |>
  pggfplot(
    filename = 'ZmCLE7_expression',
    x = mean,
    y = id
  ) |>
  pggf_bar(
    orientation = 'y',
    error = se
  )


### predicted enhancer strength of fragments upstream of ZmCLE7 ###
cur_fig <- nextsubfig(cur_fig)

## predict enhancer strength of ZmCLE7 fragments using plantGREP
# `python code/plantGREPcli/plantGREPcli.py predict -i data/extra_files/ZmCLE7_sequence.fa -o data/modelling_data/ZmCLE7_predictions.tsv`

# load data
ZmCLE7_predictions <- read_tsv('data/modelling_data/ZmCLE7_predictions.tsv') |>
  select(name, start, end, prediction_maize) |>
  mutate(
    start = start - 2484.5,
    end = end - 2483.5,
    length = end - start,
    pos = (start + end) / 2,
    across(c(start, end, length, pos), ~ .x/1000)
  )

# calculate rolling mean
ZmCLE7_predictions <- ZmCLE7_predictions |>
  mutate(
    across(
      prediction_maize,
      list(
        mean = ~ roll_mean(.x, n = 20, align = 'center', fill = NA),
        sd = ~ roll_sd(.x, n = 20, align = 'center', fill = NA)
      )
    )
  )

# find region of highest average strength
ZmCLE7_ROI <- ZmCLE7_predictions |>
  select(pos, start, end, length, prediction_maize ) |>
  arrange(pos) |>
  mutate(
    maize_mean = roll_mean(prediction_maize , n = 170, align = 'center', fill = NA)
  ) |>
  slice_max(maize_mean, n = 1) |>
  select(start, end)

# export plot data
ZmCLE7_predictions |>
  drop_na(prediction_maize_mean) |>
  pggfplot(
    filename = 'ZmCLE7_prediction_lines',
    x = pos,
    y = prediction_maize_mean
  ) |>
  pggf_line()

ZmCLE7_predictions |>
  pggfplot(
    filename = 'ZmCLE7_prediction_points',
    x = pos,
    y = prediction_maize,
    axis_annotation = list(ROIstart = ZmCLE7_ROI$start, ROIend = ZmCLE7_ROI$end)
  ) |>
  pggf_scatter()


### DeepLIFT attributions for most active SlCLV3 promoter fragment ###
cur_fig <- nextsubfig(cur_fig)

## run DeepLIFT using plantGREP for the strongest SlCLV3 promoter fragment (Solyc11g071380-9_fwd)
# `grep -A 3 'Solyc11g071380-9_fwd' data/modelling_data/SlCLV3_sequences.fa | python code/plantGREPcli/plantGREPcli.py deeplift -c dark -i - -o data/modelling_data/SlCLV3_deeplift.tsv`

# export plot data
read_tsv('data/modelling_data/SlCLV3_deeplift.tsv') |>
  mutate(
    position = 171 - position - 1531,
    base = rev_comp(base)
  ) |>
  arrange(position) |>
  pggfplot(
    filename = 'SlCLV3_DeepLIFT',
    x = position,
    y = deeplift_dark,
    preview_plot = FALSE
  ) |>
  pggf_logo(
    letters = base
  )


### DeepLIFT predictions for most active ZmCLE7 promoter region ###
cur_fig <- nextsubfig(cur_fig)

## run DeepLIFT using plantGREP for the ZmCLE7 promoter
# `python code/plantGREPcli/plantGREPcli.py deeplift -c maize -i data/extra_files/ZmCLE7_sequence.fa -o data/modelling_data/ZmCLE7_deeplift.tsv`

# export plot data
read_tsv('data/modelling_data/ZmCLE7_deeplift.tsv') |>
  mutate(
    position = position - 2484
  ) |>
  filter(between(position, ZmCLE7_ROI$start * 1000, ZmCLE7_ROI$end * 1000)) |>
  pggfplot(
    filename = 'ZmCLE7_DeepLIFT',
    x = position,
    y = deeplift_maize,
    preview_plot = FALSE
  ) |>
  pggf_logo(
    letters = base
  )


# Supplementary Table 1 ---------------------------------------------------
cur_supp_table <- nexttable(cur_supp_table)

### transcription factor motifs ###
# get length of longest motif
TF_max_len <- sapply(TF_motifs, slot, 'consensus') |> nchar() |> max()

# export all motifs
dir.create(file.path(pggf_defaults$out_dir, 'TF_logos'), showWarnings = FALSE)

TF_logo_data <- tibble(id = seq_along(TF_motifs)) |>
  group_by(id) |>
  mutate(
    logo = convert_type(TF_motifs[[id]], 'ICM')['motif'] |>
      t() |>
      as_tibble() |>
      list()
  ) |>
  unnest(logo) |>
  mutate(
    pos = seq_len(n())
  ) |>
  ungroup() |>
  pivot_longer(
    c(A, C, G, T),
    names_to = 'base',
    values_to = 'IC'
  )

TF_logo_data |>
  nest_by(id) |>
  summarise(
    data |>
      pggfplot(
        filename = paste0('TF_logos/TF_', id),
        x = pos,
        y = IC,
        xmax = TF_max_len + 0.5,
        preview_plot = FALSE,
        print_code = FALSE
      ) |>
      pggf_logo(
        letters = base
      ) |>
      print() |>
      list()
  )

TF_logo_data |>
  group_by(id) |>
  mutate(
    pos = max(pos) - pos + 1,
    base = rev_comp(base)
  ) |>
  arrange(id, pos, base) |>
  ungroup() |>
  nest_by(id) |>
  summarise(
    data |>
      pggfplot(
        filename = paste0('TF_logos/TF_', id, '_RC'),
        x = pos,
        y = IC,
        xmax = TF_max_len + 0.5,
        preview_plot = FALSE,
        print_code = FALSE
      ) |>
      pggf_logo(
        letters = base
      ) |>
      print() |>
      list()
  )


# Supplementary Data 1 ---------------------------------------------------
cur_supp_data <- nexttable(cur_supp_data)

### Plant STARR-seq results for tamsACR library ###
# select and preprocess data
table_data <- tamsACR_data |>
  pivot_wider(
    names_from = condition,
    values_from = enrichment,
    id_cols = c(species, id, orientation, GC, chromosome, start, end)
  ) |>
  left_join(
    tamsACR_seqs |>
      as.character() |>
      enframe(
        name = 'id',
        value = 'sequence'
      ),
    by = 'id'
  ) |>
  mutate(
    sequence = replace_na(sequence, ''),
    sequence = if_else(orientation == 'rev', rev_comp(sequence), sequence),
    GC = GC * 100,
    start = start + 1
  ) |>
  select(species, id, orientation, light, dark, warm, cold, maize, chromosome, start, end, GC, sequence) |>
  arrange(species, chromosome, start, end, orientation) |>
  rename('GC content (%)' = GC)

# create excel table
wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

addWorksheet(wb, sheetName = 'ACR sequence library')

writeData(
  wb,
  sheet = 1,
  paste0('Supplementary Data ', cur_supp_data, ' | Plant STARR-seq results for ACR sequence library.'),
  startCol = 1,
  startRow = 1
)

addStyle(wb, sheet = 1, style = xlsx_bold, cols = 1, rows = 1)
mergeCells(wb, sheet = 1, rows = 1, cols = seq_len(ncol(table_data)))

writeData(
  wb,
  sheet = 1,
  paste(
    "Test sequences (170 bp) were derived from accessible chromatin regions (ACRs) in the genomes of Arabidopsis (At), tomato (Sl), maize (Zm), and sorghum (Sb).",
    "The array-synthesized test sequences were cloned in the forward and reverse orientation upstream of a 35S minimal promoter driving the expression of a barcoded GFP reporter gene.",
    "The plasmid library was subjected to Plant STARR-seq in transiently transformed tobacco leaves and maize leaf protoplasts (maize).",
    "After transformation, the tobacco plants were subjected to different light (light and dark) and temperature (warm and cold) conditions before RNA extraction.",
    "Enhancer strength was determined as the enrichment of reporter mRNA over input DNA normalized to a control construct without an enhancer (noEnh; log2 set to 0)."
  ),
  startCol = 1,
  startRow = 2
)

addStyle(wb, sheet = 1, style = xlsx_wrap, rows = 2, cols = 1)
mergeCells(wb, sheet = 1, rows = 2, cols = seq_len(ncol(table_data)))
setRowHeights(wb, sheet = 1, rows = 2, heights = 12.75 * 5)

xlsx_header <- names(table_data)
enrichment_ids <- which(xlsx_header %in% condition_order)
xlsx_header[enrichment_ids] <- 'log2(enhancer strength)'

writeData(
  wb,
  sheet = 1,
  matrix(xlsx_header, nrow = 1),
  startCol = 1,
  startRow = 4,
  colNames = FALSE
)

addStyle(wb, sheet = 1, style = xlsx_bold, cols = seq_along(xlsx_header), rows = 4)

writeData(
  wb,
  sheet = 1,
  table_data,
  startCol = 1,
  startRow = 5,
  headerStyle = xlsx_bold,
  keepNA = TRUE
)

mergeCells(wb, sheet = 1, rows = 4, cols = enrichment_ids)
addStyle(wb, sheet = 1, style = xlsx_center, rows = 4, cols = enrichment_ids, stack = TRUE)

for (i in seq_along(xlsx_header)[-enrichment_ids]) {
  mergeCells(wb, sheet = 1, rows = 4:5, cols = i)
}

addStyle(wb, sheet = 1, style = xlsx_2digit, cols = enrichment_ids, rows = 6:(nrow(table_data) + 5), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = xlsx_2digit, cols = which(xlsx_header == 'GC content (%)'), rows = 6:(nrow(table_data) + 5), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = xlsx_seq_font, cols = which(xlsx_header == 'sequence'), rows = 6:(nrow(table_data) + 5))

setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'chromosome'), widths = 12)
setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'GC content (%)'), widths = 13.6)

freezePane(wb, sheet = 1, firstActiveRow = 6)

saveWorkbook(wb, paste0(supp_data_dir, '/SupplementaryData', cur_supp_data, '.xlsx'), overwrite = TRUE)


# Supplementary Data 2 ---------------------------------------------------
cur_supp_data <- nexttable(cur_supp_data)

### GO terms for the genes linked to the strongest ACRs ###
# order GO terms
GOterm_order <- GOterms |>
  distinct(term_id, term_name) |>
  arrange(term_id) |>
  mutate(
    term_name = str_replace(term_name, '_', ' '),
    GO_term = paste0(term_name, ' (', term_id, ')')
  ) |>
  select(term_id, GO_term) |>
  deframe()

# select and preprocess data
table_data <- GOterms |>
  mutate(
    GO_term = ordered(term_id, levels = names(GOterm_order), labels = GOterm_order)
  ) |>
  select(condition, species, GO_term, p_value) |>
  pivot_wider(
    names_from = c(condition, species),
    values_from = p_value
  ) |>
  mutate(
    across(-GO_term, ~ replace_na(.x, 1))
  ) |>
  arrange(GO_term)

# create excel table
wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

addWorksheet(wb, sheetName = 'GO terms')

writeData(
  wb,
  sheet = 1,
  paste0('Supplementary Data ', cur_supp_data, ' | GO term enrichment for genes linked to the strongest test sequences.'),
  startCol = 1,
  startRow = 1
)

addStyle(wb, sheet = 1, style = xlsx_bold, cols = 1, rows = 1)
mergeCells(wb, sheet = 1, rows = 1, cols = seq_len(ncol(table_data)))

writeData(
  wb,
  sheet = 1,
  paste(
    "Test sequences from Arabidopsis (At), tomato (Sl), maize (Zm), and sorghum (Sb) were ordered by their enhancer strength (average across fwd and rev orientation).",
    "GO term enrichment for genes linked to the strongest 5% of the test sequences from each species was determined for the indicated conditions and species.",
    "Adjusted p values are shown below. Significant (adjusted p value ≤ 0.05) enrichment is indicated in green."
  ),
  startCol = 1,
  startRow = 2
)

addStyle(wb, sheet = 1, style = xlsx_wrap, rows = 2, cols = 1)
mergeCells(wb, sheet = 1, rows = 2, cols = seq_len(ncol(table_data)))
setRowHeights(wb, sheet = 1, rows = 2, heights = 12.75 * 2)

xlsx_header <- names(table_data)

writeData(
  wb,
  sheet = 1,
  matrix(c(xlsx_header[1], str_replace(xlsx_header[-1], '_.*', '')), nrow = 1),
  startCol = 1,
  startRow = 4,
  colNames = FALSE
)

writeData(
  wb,
  sheet = 1,
  matrix(c(xlsx_header[1], str_replace(xlsx_header[-1], '.*_', '')), nrow = 1),
  startCol = 1,
  startRow = 5,
  colNames = FALSE
)

addStyle(wb, sheet = 1, style = xlsx_bold, cols = seq_along(xlsx_header), rows = 4:5, gridExpand = TRUE)
for (i in seq(2, length(xlsx_header), 4)) {
  mergeCells(wb, sheet = 1, rows = 4, cols = i:(i+3))
}
addStyle(wb, sheet = 1, style = xlsx_center, rows = 4, cols = 2:length(xlsx_header), stack = TRUE)
mergeCells(wb, sheet = 1, rows = 4:5, cols = 1)

writeData(
  wb,
  sheet = 1,
  table_data,
  startCol = 1,
  startRow = 6,
  colNames = FALSE,
  keepNA = TRUE
)

addStyle(wb, sheet = 1, style = xlsx_right, cols = 1, rows = 6:(nrow(table_data) + 5), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = xlsx_mixed, cols = 2:length(xlsx_header), rows = 6:(nrow(table_data) + 5), gridExpand = TRUE)
conditionalFormatting(wb, sheet = 1, style = xlsx_green, cols = 2:length(xlsx_header), rows = 6:(nrow(table_data) + 5), gridExpand = TRUE, rule = '<0.05')

setColWidths(wb, sheet = 1, cols = 1, widths = 58)

freezePane(wb, sheet = 1, firstActiveRow = 6)

saveWorkbook(wb, paste0(supp_data_dir, '/SupplementaryData', cur_supp_data, '.xlsx'), overwrite = TRUE)


# Supplementary Data 3 ---------------------------------------------------
cur_supp_data <- nexttable(cur_supp_data)

### Plant STARR-seq results for validation library ###
# helper functions to generate experiment names
TFBS_insertion_name <- function(ids) {
  exp <- str_match(ids, 'rnd-[dm][0-9]+_.*_(.*)_(.*)_(.*)')[,-1] |>
    as_tibble(.name_repair = ~ paste0('TF', seq_len(3))) |>
    mutate(
      across(everything(), ~ if_else(.x == 'NA', 'none', paste0('TF_', .x)))
    ) |>
    unite(
      'TFs',
      everything(),
      sep = ', ',
      na.rm = TRUE
    ) |>
    mutate(
      TFs = paste0('TFBS insertion: ', TFs)
    ) |>
    pull(TFs)
  
  return(exp)
}

evolution_name <- function(ids) {
  exp <- tibble(id = ids) |>
    separate_wider_regex(
      id,
      patterns = c(
        '.*_evo',
        mpr = '.*',
        '_',
        obj = '.*',
        '_',
        rnd = '.*'
      ),
      too_few = 'align_start'
    ) |>
    mutate(
      exp = paste0('evolution: round = ', rnd, ', objective = ', obj, ', mutations/round = ', mpr)
    ) |>
    pull(exp)
  
  return(exp)
}

# select and preprocess data
table_data <- eVal_data |>
  mutate(
    experiment = case_when(
      id == 'noEnh' ~ 'no-enhancer control',
      str_detect(id, '_(fwd|rev)$') ~ 'validation of ACR sequence library',
      str_detect(id, '_shuffle$') ~ paste0('TFBS shuffling: ', str_extract(eVal_data$id, '.*_(TF_[0-9]+)_.*', group = 1)),
      str_detect(id, 'rnd-[dm][0-9]+_.*_.*_.*_.*') ~ TFBS_insertion_name(id),
      str_detect(id, '_WT$') ~ 'unmodified control',
      str_detect(id, 'evo0_start_0') ~ 'unmodified control',
      str_detect(id, '_evo') ~ evolution_name(id)
    ),
    base_seq = str_replace(id, '_evo.*', ''),
    base_seq = str_replace(base_seq, '_TF_.*', '_fwd'),
    base_seq = str_replace(base_seq, '_.*_.*_.*_.*', '')
  ) |>
  pivot_wider(
    names_from = condition,
    values_from = enrichment,
    id_cols = c(id, base_seq, experiment, GC)
  ) |>
  left_join(
    validation_seqs |>
      as.character() |>
      enframe(
        name = 'id',
        value = 'sequence'
      ),
    by = 'id'
  ) |>
  mutate(
    sequence = replace_na(sequence, ''),
    GC = GC * 100
  ) |>
  select(base_seq, experiment, light, dark, warm, cold, maize, GC, sequence) |>
  arrange(base_seq, experiment) |>
  rename('base sequence' = base_seq, 'GC content (%)' = GC)

# create excel table
wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

addWorksheet(wb, sheetName = 'validation library')

writeData(
  wb,
  sheet = 1,
  paste0('Supplementary Data ', cur_supp_data, ' | Plant STARR-seq results for validation library.'),
  startCol = 1,
  startRow = 1
)

addStyle(wb, sheet = 1, style = xlsx_bold, cols = 1, rows = 1)
mergeCells(wb, sheet = 1, rows = 1, cols = seq_len(ncol(table_data)))

writeData(
  wb,
  sheet = 1,
  paste(
    "Test sequences (170 bp) with either shuffled or added transcription factor binding sites as well as sequences evolved in silico",
    "using plantGREP were array-synthesized and cloned upstream of a 35S minimal promoter driving the expression of a barcoded GFP reporter gene.",
    "The plasmid library was subjected to Plant STARR-seq in transiently transformed tobacco leaves and maize leaf protoplasts (maize).",
    "After transformation, the tobacco plants were subjected to different light (light and dark) and temperature (warm and cold) conditions before RNA extraction.",
    "Enhancer strength was determined as the enrichment of reporter mRNA over input DNA normalized to a control construct without an enhancer (noEnh; log2 set to 0)."
  ),
  startCol = 1,
  startRow = 2
)

addStyle(wb, sheet = 1, style = xlsx_wrap, rows = 2, cols = 1)
mergeCells(wb, sheet = 1, rows = 2, cols = seq_len(ncol(table_data)))
setRowHeights(wb, sheet = 1, rows = 2, heights = 12.75 * 5)

xlsx_header <- names(table_data)
enrichment_ids <- which(xlsx_header %in% condition_order)
xlsx_header[enrichment_ids] <- 'log2(enhancer strength)'

writeData(
  wb,
  sheet = 1,
  matrix(xlsx_header, nrow = 1),
  startCol = 1,
  startRow = 4,
  colNames = FALSE
)

addStyle(wb, sheet = 1, style = xlsx_bold, cols = seq_along(xlsx_header), rows = 4)

writeData(
  wb,
  sheet = 1,
  table_data,
  startCol = 1,
  startRow = 5,
  headerStyle = xlsx_bold,
  keepNA = TRUE
)

mergeCells(wb, sheet = 1, rows = 4, cols = enrichment_ids)
addStyle(wb, sheet = 1, style = xlsx_center, rows = 4, cols = enrichment_ids, stack = TRUE)

for (i in seq_along(xlsx_header)[-enrichment_ids]) {
  mergeCells(wb, sheet = 1, rows = 4:5, cols = i)
}

addStyle(wb, sheet = 1, style = xlsx_2digit, cols = enrichment_ids, rows = 6:(nrow(table_data) + 5), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = xlsx_2digit, cols = which(xlsx_header == 'GC content (%)'), rows = 6:(nrow(table_data) + 5), gridExpand = TRUE)
addStyle(wb, sheet = 1, style = xlsx_seq_font, cols = which(xlsx_header == 'sequence'), rows = 6:(nrow(table_data) + 5))

setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'base sequence'), widths = 22.5)
setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'experiment'), widths = 33)
setColWidths(wb, sheet = 1, cols = which(xlsx_header == 'GC content (%)'), widths = 13.6)

freezePane(wb, sheet = 1, firstActiveRow = 6)

saveWorkbook(wb, paste0(supp_data_dir, '/SupplementaryData', cur_supp_data, '.xlsx'), overwrite = TRUE)


# Supplementary Data 4 ---------------------------------------------------
cur_supp_data <- nexttable(cur_supp_data)

### summary of statistical tests ###
# helper to rename TFs
TF_names <- TF_motifs |>
  sapply(function(x) x['altname']) |>
  str_replace('_', '-') |>
  as_tibble() |>
  mutate(
    id = paste0('TF_', seq_len(n())),
    family = paste0(value, ' (', seq_len(n()), ')')
  ) |>
  select(id, family) |>
  arrange(-nchar(id)) |>
  deframe()

# get test results
test_results <- pggf_get_test_results()

# select and preprocess data
table_data <- test_results |>
  unite(
    'subplot',
    col,
    row,
    sep = ', ',
    na.rm = TRUE
  ) |>
  mutate(
    subplot = if_else(subplot == '' & ! is.na(facet_id), condition_order[facet_id], subplot),
    fig_type = ordered(str_sub(name, 1, 1), levels = c('F', 'E', 'S')),
    fig_num = as.integer(str_extract(name, '.*Fig. ([0-9]+).*', group = 1)),
    subfig = str_extract(name, '.*Fig. [0-9]+(.*)', group = 1),
    p_adjust_method = replace_na(p_adjust_method, ''),
    p_adjust_method = str_replace(p_adjust_method, 'none', ''),
    comparison = str_replace(comparison, '.([0-9]+),([0-9]+). vs. .([0-9]+),([0-9]+).', '\\1-\\2 vs. \\3-\\4'),
    comparison = str_replace(comparison, "'", "' "),
    comparison = str_replace_all(
      comparison,
      c('u5\\+' = '> 5 kb upstream', 'u5' = '5-1 kb upstream', 'u1' = '1-0 kb upstream', 'd5\\+' = '> 5 kb downstream', 'd5' = '1-5 kb downstream', 'd1' = '0-1 kb downstream')
    ),
    comparison = str_replace_all(comparison, TF_names),
    comparison = str_replace(comparison, 'noTFBS', 'none'),
    comparison = str_replace_all(comparison, c('acrossPos' = 'across positions', 'acrossSeq' = 'across sequences')),
    comparison = str_replace_all(comparison, c('0.76667' = '6', '1.86667' = '12', 'mpr1' = '1 mutation/round', 'mpr2' = '2 mutations/round', 'mpr3' = '3 mutations/round')),
    comparison = case_match(
      subplot,
      'species-specificity' ~ str_replace_all(comparison, c('-12' = '12:tobacco-specific', '-6' = '6:tobacco-specific', '12(?!:)' = '12:maize-specific', '6(?!:)' = '6:maize-specific')),
      'light-specificity' ~ str_replace_all(comparison, c('-12' = '12:light-specific', '-6' = '6:light-specific', '12(?!:)' = '12:dark-specific', '6(?!:)' = '6:dark-specific')),
      'temperature-specificity' ~ str_replace_all(comparison, c('-12' = '12:warm-specific', '-6' = '6:warm-specific', '12(?!:)' = '12:cold-specific', '6(?!:)' = '6:cold-specific')),
      .default = comparison
    ),
  ) |>
  arrange(fig_type, fig_num, subfig) |>
  select('figure' = name, subplot, test, comparison, 'correction' = p_adjust_method, 'p value' = p_value)

# create excel table
wb <- createWorkbook()

modifyBaseFont(wb, fontSize = 10, fontName = 'Arial')

addWorksheet(wb, sheetName = 'p values')

writeData(
  wb,
  sheet = 1,
  paste0('Supplementary Data ', cur_supp_data, ' | Exact p values for significance levels and groups shown in the figures.'),
  startCol = 1,
  startRow = 1
)

addStyle(wb, sheet = 1, style = xlsx_bold, cols = 1, rows = 1)
mergeCells(wb, sheet = 1, rows = 1, cols = seq_len(ncol(table_data)))

writeData(
  wb,
  sheet = 1,
  paste(
    "Significant (adjusted p value ≤ 0.05) results are indicated in green.",
    "Because of technical limitations, p values below 1E‑16 (Tukey tests) or 1E‑323 (other tests) appear as 0.00E+0."
  ),
  startCol = 1,
  startRow = 2
)

addStyle(wb, sheet = 1, style = xlsx_wrap, rows = 2, cols = 1)
mergeCells(wb, sheet = 1, rows = 2, cols = seq_len(ncol(table_data)))

writeData(
  wb,
  sheet = 1,
  table_data,
  startCol = 1,
  startRow = 4,
  headerStyle = xlsx_bold,
  keepNA = TRUE
)

p_col <- which(names(table_data) == 'p value')
test_col <- which(names(table_data) == 'test')

addStyle(wb, sheet = 1, style = xlsx_mixed, cols = p_col, rows = 5:(nrow(table_data) + 4), gridExpand = TRUE)
conditionalFormatting(wb, sheet = 1, style = xlsx_green, cols = p_col, rows = 5:(nrow(table_data) + 4), gridExpand = TRUE, rule = '<0.05')

setColWidths(wb, sheet = 1, cols = which(names(table_data) == 'figure'), widths = 20.5)
setColWidths(wb, sheet = 1, cols = which(names(table_data) == 'test'), widths = 62.5)
setColWidths(wb, sheet = 1, cols = which(names(table_data) == 'comparison'), widths = 67.5)

freezePane(wb, sheet = 1, firstActiveRow = 5)

saveWorkbook(wb, paste0(supp_data_dir, '/SupplementaryData', cur_supp_data, '.xlsx'), overwrite = TRUE)


# wrap-up -----------------------------------------------------------------

### reset code printing and plot previewing after sourcing the script ###
if (sys.nframe() != 0) {
  pggf_config(print_code = TRUE, preview_plot = TRUE)
  options(readr.show_col_types = NULL)
}
