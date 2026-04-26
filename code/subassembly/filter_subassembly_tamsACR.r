### process arguments ###
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop('At least three arguments are required:
    1. subassembly of ACR sequences in forward orientation,
    2. subassembly of ACR sequences in reverse orientation,
    3. barcodes of the control constructs,
    4. length of the ACR sequences (optional; default = 170)')
} else if (length(args) == 3) {
  args[4] <- 170
} else if (length(args) > 4) {
  warning(paste0('More than four arguments passed. Ignoring arguments 5 - ', length(args), '.'))
}

subass_fwd <- args[1]
subass_rev <- args[2]
control_bcs <- args[3]
length <- args[4]


### test if files exist ###
if (! file.exists(subass_fwd)) {
  stop(paste0('The file "', subass_fwd, '" does not exist.'))
}
if (! file.exists(subass_rev)) {
  stop(paste0('The file "', subass_rev, '" does not exist.'))
}
if (! file.exists(control_bcs)) {
  stop(paste0('The file "', control_bcs, '" does not exist.'))
}


### load libraries ###
library(readr)
library(dplyr)
library(tidyr)


### load subassembly ###
subassembly <- bind_rows('fwd' = read_tsv(subass_fwd), 'rev' = read_tsv(subass_rev), .id = 'orientation') %>%
  mutate(
    truncation = case_when(
      start > 1 & stop < length ~ paste0('1-', start - 1, '&', stop + 1,'-', length),
      start > 1 ~ paste0('1-', start - 1),
      stop < length ~ paste0(stop + 1,'-', length),
      TRUE ~ NA_character_
    ),
    length = stop - start + 1
  ) %>%
  select(-start, -stop)


### add controls to subassembly ###
controls <- read_tsv(control_bcs) %>%
  rename('id' = variant) %>%
  mutate(
    variant = 'control'
  )

subassembly <- subassembly %>%
  bind_rows(controls)


### find barcodes linked to a single ACR ###
subassembly_okay <- subassembly %>%
  group_by(barcode) %>%
  filter(n() == 1) %>%
  ungroup() %>%
  select(-assembly_count)


### find barcodes with more than one associated fragment ###
subassembly_to_fix <- subassembly %>%
  filter(! barcode %in% subassembly_okay$barcode)


### combine barcodes with multiple variants of the same ACR ###
# use most common variant if it is at least 10 times more common than the second most common one
subassembly_to_fix <- subassembly_to_fix %>%
  group_by(across(c(-variant, -truncation, -assembly_count))) %>%
  arrange(desc(assembly_count)) %>%
  summarise(
    across(
      c(variant, truncation),
      function(x) case_when(
        n_distinct(x) == 1 ~ first(x),
        first(assembly_count) >= 10 * nth(assembly_count, 2) ~ first(x),
        TRUE ~ paste('multiple (', paste(x, collapse = '/'), ')', sep = '')
      )
    ),
    assembly_count = first(assembly_count)
  ) %>%
  ungroup()


### combine barcodes with multiple ACRs ###
# use most common ACR if it is at least 20 times more common than the second most common one
subassembly_to_fix <- subassembly_to_fix %>%
  group_by(barcode) %>%
  arrange(desc(assembly_count)) %>%
  summarise(
    id = case_when(
      n() == 1 ~ first(id),
      first(assembly_count) >= 20 * nth(assembly_count, 2) ~ first(id),
      TRUE ~ paste('multiple (', paste(id, collapse = '/'), ')', sep = '')
    ),
    across(c(-id, -assembly_count), function(x) ifelse(grepl('multiple', id), NA, first(x)))
  ) %>%
  ungroup()


### combine subassembly_okay and subassembly_to_fix ###
subassembly <- bind_rows(
  subassembly_okay,
  subassembly_to_fix
)


### filter out barcodes linked to multiple ACRs ###
subassembly <- subassembly %>%
  filter(! grepl('multiple', id, fixed = TRUE) & ! grepl('multiple', variant, fixed = TRUE))


### filter out mutated or truncated sequences ###
subassembly <- subassembly %>%
  filter(variant %in% c('WT', 'control') & is.na(truncation)) |>
  select(-variant, -truncation)


### write filtered subassembly to standard out ###
subassembly %>%
  format_tsv() %>%
  write(stdout())

