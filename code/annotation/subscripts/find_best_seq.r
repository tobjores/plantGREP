### process arguments ###
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("At least four arguments required: 1. species id; 2. fasta file for original ACRs; 3. fasta file for variant ACRs (first option); 4. fasta file for variant ACRs (second option); 5. Levenshtein distance cut-off (optional)", call. = FALSE)
} else if (length(args) == 4) {
  args[5] <- NA
} else if (length(args) > 5) {
  warning(paste0('More than five arguments passed. Ignoring arguments 6 - ', length(args), '.'))
}


### load libraries ###
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(Biostrings)


### load sequences ###
spid <- args[1]
orig <- readDNAStringSet(args[2])
opt1 <- readDNAStringSet(args[3])
opt2 <- readDNAStringSet(args[4])
cutoff <- as.integer(args[5])


### find sequences without alternatives ###
opt1.shortname <- substr(names(opt1), 1, str_locate(names(opt1), '::') - 1)
opt2.shortname <- substr(names(opt2), 1, str_locate(names(opt2), '::') - 1)

opt1.only <- opt1[! opt1.shortname %in% opt2.shortname]
opt2.only <- opt2[! opt2.shortname %in% opt1.shortname]

opt1 <- c(opt1, opt2.only)
opt2 <- c(opt2, opt1.only)


### filter original sequences to only the optional sequences and bring all in same order ###
opts.shortname <- unique(c(opt1.shortname, opt2.shortname))
opts.shortname <- opts.shortname[order(as.integer(str_sub(opts.shortname, 4)))]

orig <- orig[opts.shortname]
opt1 <- opt1[match(opts.shortname, substr(names(opt1), 1, str_locate(names(opt1), '::') - 1))]
opt2 <- opt2[match(opts.shortname, substr(names(opt2), 1, str_locate(names(opt2), '::') - 1))]


### remove sequences with degenerate bases ###
opt1.deg <- which(! grepl('^[ACGT]*$', opt1))
opt2.deg <- which(! grepl('^[ACGT]*$', opt2))

both.deg <- opt1.deg[opt1.deg %in% opt2.deg]
opt1.deg.only <- opt1.deg[! opt1.deg %in% opt2.deg]
opt2.deg.only <- opt2.deg[! opt2.deg %in% opt1.deg]

opt1[opt1.deg.only] <- opt2[opt1.deg.only]
names(opt1)[opt1.deg.only] <- names(opt2)[opt1.deg.only]
opt2[opt2.deg.only] <- opt1[opt2.deg.only]
names(opt2)[opt2.deg.only] <- names(opt1)[opt2.deg.only]

if (length(both.deg) > 0) {
  opt1 <- opt1[-both.deg]
  opt2 <- opt2[-both.deg]
  orig <- orig[-both.deg]
  opts.shortname <- opts.shortname[-both.deg]
}


### find identical sequences ###
same <- which(opt1 == opt2)

same.seq <- opt1[same]

same.orig <- orig[same]

same.name <- names(same.seq)
same.dist <- sapply(seq_along(same.orig), function(x) stringDist(c(same.orig[x], same.seq[x]))[1])

if (length(same) > 0) {
  opt1 <- opt1[-same]
  opt2 <- opt2[-same]
  orig <- orig[-same]
  opts.shortname <- opts.shortname[-same]
}


### find option with lowest Levenshtein distance to original sequence ###
opt.dist <- lapply(seq_along(orig), function(x) stringDist(c(orig[x], opt1[x], opt2[x]))[-3])

best.opt <- sapply(opt.dist, which.min)


### for options with identical Levenshtein distance use the Hamming distance as a tie-breaker (uses first option if there is still a tie) ###
opt.dist.same <- sapply(opt.dist, function(x) x[1] == x[2])

opt1.same <- opt1[opt.dist.same]
opt2.same <- opt2[opt.dist.same]
orig.same <- orig[opt.dist.same]

best.opt.same <- sapply(seq_along(orig.same), function(x) which.min(stringDist(c(orig.same[x], opt1.same[x], opt2.same[x]), method = 'hamming')[-3]))

best.opt[opt.dist.same] <- best.opt.same


### get best sequence ###
best.seq <- c(opt1[best.opt == 1], opt2[best.opt == 2])
best.seq <- best.seq[match(opts.shortname, substr(names(best.seq), 1, str_locate(names(best.seq), '::') - 1))]

best.name <- names(best.seq)
best.dist <- sapply(opt.dist, min)


### combine data ###
all.data <- tibble(
  name = as.character(c(same.name, best.name)),
  dist = as.integer(c(same.dist, best.dist)),
  sequence = as.character(c(same.seq, best.seq))
)


### remove sequences with high Levenshtein distance (optional) ###
if (! is.na(cutoff)) {
  all.data <- all.data %>%
    filter(dist <= cutoff)
}


### rename sequences ###
all.data <- all.data %>%
  separate(
    name,
    into = c('name', 'strand', NA),
    sep = '[()]'
  ) %>%
  separate(
    name,
    into = c('sp', 'id', NA, 'chr', 'start', 'end'),
    sep = '[:-]'
  ) %>%
  arrange(as.numeric(id)) %>%
  mutate(
    id = paste0(sp, '-', id, '(', spid, ')')
  ) %>%
  select(-sp)


### save data ###
all.data %>%
  format_tsv(col_names = TRUE) %>%
  write(stdout())
