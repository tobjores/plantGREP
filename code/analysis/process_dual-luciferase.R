library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

### load and preprocess dual-luciferase data ###
# file:         character or connection; plate reader file with dual-luciferase  data
# luminescence: character; region of the excel file with luminescence data
# constructs:   character; region of the excel file with construct data
load_experiment_dl <- function(file, luminescence = 'B48:O72', lines = 'B75:O99') {
  # process excel cell ranges
  lum_loc <- str_match(luminescence, '([A-Z]+)([0-9]+):([A-Z]+)([0-9]+)')[-1]
  lines_loc <- str_match(lines, '([A-Z]+)([0-9]+):([A-Z]+)([0-9]+)')[-1]
  
  # load and restructure construct data
  data_constructs <- read.xlsx(
    xlsxFile = file,
    cols = seq(col2int(lines_loc[1]), col2int(lines_loc[3])),
    rows = seq(lines_loc[2], lines_loc[4])
  ) |>
    as_tibble(
      .name_repair = 'unique'
    ) |>
    rename('row' = 1, 'type' = last_col()) |>
    mutate(
      row = coalesce(row, lag(row)),
      row = coalesce(row, lag(row)),
      across(everything(), as.character)
    ) |>
    drop_na(row) |>
    pivot_longer(
      cols = as.character(1:12),
      names_to = 'column',
      values_to = 'value'
    ) |>
    unite(
      col = 'well',
      row,
      column,
      sep = ''
    ) |>
    drop_na(value) |>
    pivot_wider(
      names_from = type,
      values_from = value
    )
  
  # load and restructure luminescence data
  data_luminescence <- read.xlsx(
    xlsxFile = file,
    cols = seq(col2int(lum_loc[1]), col2int(lum_loc[3])),
    rows = seq(lum_loc[2], lum_loc[4])
  ) |>
    as_tibble(
      .name_repair = 'unique'
    ) |>
    rename('row' = 1, 'type' = last_col()) |>
    filter(type != 'NLuc/Luc ratio') |>
    mutate(
      row = coalesce(row, lag(row)),
      type = str_replace(type, ':Lum', ''),
      across(as.character(1:12), as.numeric)
    ) |>
    pivot_longer(
      cols = as.character(1:12),
      names_to = 'column',
      values_to = 'value'
    ) |>
    unite(
      col = 'well',
      row,
      column,
      sep = ''
    ) |>
    pivot_wider(
      names_from = type,
      values_from = value
    )
  
  # merge construct and luminescence data
  data_dl <- inner_join(
    data_luminescence,
    data_constructs,
    by = 'well'
  )
  
  # return result
  return(data_dl)
}

### load and process data ###
## load data
data_dl_all <- expand_grid(experiment = seq_len(3)) |>
  nest_by(experiment) |>
  reframe(
    load_experiment_dl(paste0('data/plate_reader_data/dual-luciferase_exp', experiment, '.xlsx'))
  ) |>
  ungroup()

## background-subtract luminescence and calculate NLuc/Luc ratio
data_dl_ratio <- data_dl_all |>
  group_by(experiment) |>
  mutate(
    across(
      contains('Luc'),
      ~ .x - mean(.x[line == 'blank'])
    )
  ) |>
  ungroup() |>
  filter(! line %in% c('blank', 'WT') & Luc > 0 & NLuc > 0) |>
  mutate(
    l2ratio = log2(NLuc/Luc)
  )

## normalize NLuc/Luc ratio to control construct without enhancer
data_dl_norm <- data_dl_ratio |>
  group_by(experiment) |>
  mutate(
    l2ratio = l2ratio - median(l2ratio[enhancer == 'noEnh'])
  ) |>
  ungroup()

### save data ###
saveRDS(data_dl_norm, file = 'data/RData/DL_data.rds')