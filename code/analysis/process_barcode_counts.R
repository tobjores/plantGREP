library(readr)
library(dplyr)
library(tidyr)
library(stringr)

### load and preprocess experiment data (barcode count files) ###
# input:        character or connection; data file with input read counts
# output:       named vector with conditon and output read count file; <condition 1> = <data file 1>, <condition 2> = <data file 2>, ...
# subassembly:  R object; subassembly data frame
# rcc:          numeric; read count cutoff to be applied to the input and output counts
load_experiment_bc <- function(input, output, subassembly, rcc = 1) {
  # load and merge files (only keep barcodes detected in both input and output)
  data_inp <- read_table(input, col_names = c('count', 'barcode'), col_types = 'ic') |>
    filter(count >= rcc)
  
  data_out <- tibble(condition = names(output), file = unlist(output)) |>
    nest_by(condition) |>
    reframe(read_table(data$file, col_names = c('count', 'barcode'), col_types = 'ic')) |>
    ungroup() |>
    filter(count >= rcc) |>
    inner_join(data_inp, by = 'barcode', suffix = c('_out', '_inp'))
  
  # merge with subassembly
  data_out <- inner_join(data_out, subassembly, by = 'barcode')
  
  # return result
  return(data_out)
}


### load experiment metadata ###
experiments <- read_tsv('data/experiments.tsv', comment = '#') |>
  select(library, experiment, conditions) |>
  separate_rows(
    conditions,
    sep = '/'
  )


### process barcode counts ###
for (lib in unique(experiments$library)) {
  # load subassembly
  subassembly <- read_tsv(paste0('data/subassembly/subassembly_', lib, '_filtered.tsv.gz'), col_types = if_else(lib == 'eVal', 'ccdd', 'cccdd'))
  
  # load barcode data
  barcode_counts <- experiments |>
    filter(library == lib) |>
    group_by(experiment) |>
    summarise(
      conditions = list(setNames(paste0('data/barcode_counts/', lib, '/exp', experiment,'/barcode_', lib, '_', experiment, '_', conditions, '.count.gz'), conditions))
    ) |>
    ungroup() |>
    nest_by(experiment) |>
    reframe(
      exp_helper = if_else(experiment == 'O', 'N', experiment),
      load_experiment_bc(
        input = paste0('data/barcode_counts/', lib, '/exp', exp_helper, '/barcode_', lib, '_', exp_helper, '_input.count.gz'),
        output = unlist(data$conditions),
        subassembly = subassembly,
        rcc = 5
      )
    ) |>
    select(-exp_helper) |>
    ungroup()
  
  # rename conditions in tobacco
  barcode_counts <- barcode_counts |>
    mutate(
      condition = str_replace(condition, '35C', 'warm'),
      condition = str_replace(condition, '10C', 'cold')
    )
  
  # select controls, calculate barcode enrichment and save the results
  barcode_counts |>
    filter(substr(id, 1, 3) %in% c('35S', 'AB8', 'Cab', 'rbc', 'noE', 'Bsa')) |>
    group_by(experiment, condition) |>
    mutate(
      enrichment = log2((count_out / sum(count_out)) / (count_inp / sum(count_inp))),
      enrichment = enrichment - median(enrichment[id == 'noEnh'])
    ) |>
    ungroup() |>
    select(experiment, condition, barcode, id, enrichment) |>
    saveRDS(paste0('data/RData/', lib, '_controls.rds'))
  
  # aggregate barcodes (sum of individual counts)
  barcodes_ag <- barcode_counts |>
    group_by(across(-c(count_inp, count_out, barcode))) |>
    summarise(
      n_bc = n(),
      count_inp = sum(count_inp),
      count_out = sum(count_out)
    ) |>
    ungroup()
  
  # calculate enrichment and normalize first to library median then to mean strength of noEnh-control
  barcodes_norm <- barcodes_ag |>
    group_by(experiment, condition) |>
    mutate(
      enrichment = log2((count_out / sum(count_out)) / (count_inp / sum(count_inp))),
      enrichment = enrichment - median(enrichment)
    ) |>
    group_by(condition) |>
    mutate(
      enrichment = enrichment - mean(enrichment[id == 'noEnh'])
    ) |>
    ungroup()
  
  # calculate mean across replicates
  barcodes_mean <- barcodes_norm |>
    group_by(across(-c(experiment, n_bc, count_inp, count_out, enrichment))) |>
    summarise(
      n_experiments = n(),
      min_bc = min(n_bc),
      min_ci = min(count_inp),
      min_co = min(count_out),
      enrichment = mean(enrichment)
    ) |>
    ungroup()
  
  # save data
  saveRDS(barcodes_norm, file = paste0('data/RData/', lib, '_reps.rds'))
  saveRDS(barcodes_mean, file = paste0('data/RData/', lib, '_main.rds'))
}