### process arguments ###
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop('A data file is required.')
} else if (length(args) > 1) {
  warning(paste0('More than one arguments passed. Ignoring arguments 2 - ', length(args), '.'))
  infile <- args[1]
} else {
  infile <- args[1]
}


### test if data file exists ###
if (! file.exists(infile)) {
  stop(paste0('The file "', infile, '" does not exist.'))
}


### load libraries ###
library(readr)
library(readxl)
library(dplyr)


### read, process and write data ###
peaks <- read_xlsx(infile, sheet = 1, range = 'B8:D9399') %>%
  arrange(Chr, Start, End) %>%
  distinct() %>%
  format_tsv(col_names = FALSE) %>%
  write(stdout())