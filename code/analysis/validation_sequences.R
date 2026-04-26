library(tidyverse)
library(universalmotif)

### load tamsACR data and sequences ###
tamsACR_data <- readRDS('data/RData/tamsACR_main.rds')
  
tamsACR_seqs <- Biostrings::readDNAStringSet('data/refseq/tamsACR_sequences.fa.gz')


### load transcription factor data ###
TF_motifs <- read_meme('data/extra_files/TF-clusters.meme')

target_TFs <- read_tsv('data/extra_files/TFs_for_validation.txt')


### select sequences for TF validation (wildtype sequences and sequences with shuffled motifs) ###
# function to select target sequences
get_TFval_seqs <- function(TF, n_seqs) {
  if (TF %in% c(6, 19)) {# palindromic motifs
    n_hits <- 2
  } else {# non-palindromic motifs
    n_hits <- 1
  }
  
  # find sequences with only one high-confidence binding site
  one_hit_seqs <- scan_sequences(TF_motifs[[TF]], tamsACR_seqs, RC = TRUE, threshold = 0.0001, nthreads = 0) |>
    as_tibble() |>
    count('id' = sequence) |>
    filter(n == n_hits) |>
    pull(id)
  
  # select sequences with only very weak secondary binding sites
  target_seqs <- tamsACR_seqs[names(tamsACR_seqs) %in% one_hit_seqs]
  
  target_seqs <- scan_sequences(TF_motifs[[TF]], target_seqs, RC = TRUE, threshold = 0.1, nthreads = 0) |>
    as_tibble() |>
    group_by(sequence) |>
    filter(n() > n_hits) |> # for each sequence select the second best hit
    slice_max(score, n = n_hits + 1, with_ties = FALSE) |>
    slice_min(score, n = 1, with_ties = FALSE) |>
    ungroup() |>
    slice_sample(prop = 1) |> # shuffle sequences
    slice_min(score, n = n_seqs, with_ties = FALSE) |> # select sequences with weakest secondary hits
    pull(sequence)
  
  return(target_seqs)
}

# function to identify motifs and replace them with a shuffled version
motif_shuffle <- function(seq, TF) {
  # find motif
  hit <- scan_sequences(TF_motifs[[TF]], Biostrings::DNAStringSet(seq), RC = TRUE, threshold = 0.0001, nthreads = 0) |>
    _[1, c('start', 'stop', 'strand', 'match'), drop = TRUE]
  
  if (hit$strand == '-') {
    hit[c('start', 'stop')] <- hit[c('stop', 'start')]
    hit$match <- hit$match |>
      Biostrings::DNAString() |>
      Biostrings::reverseComplement() |>
      as.character()
  }
  
  # create sequences with shuffled motif
  shuf_match <- sapply(rep(hit[['match']], 100), shuffle_string, USE.NAMES = FALSE)
  shuf_match <- shuf_match[shuf_match != hit$match]
  
  shuf_seqs <- str_replace(as.character(seq), hit[['match']], shuf_match)
  
  # remove shuffled sequences with a BsaI site
  shuf_seqs <- shuf_seqs[! str_detect(shuf_seqs, '(GGTCT(C|$))|(GAGAC(C|$))')]
  
  # score shuffled sequences and select the one with the weakest secondary match
  seq_id <- scan_sequences(TF_motifs[[TF]], Biostrings::DNAStringSet(shuf_seqs), RC = TRUE, threshold = 0.1, nthreads = 0) |>
    as_tibble() |>
    group_by(sequence.i) |>
    summarise(
      score = max(score)
    ) |>
    ungroup() |>
    slice_min(score, n = 1, with_ties = FALSE) |>
    pull(sequence.i)
  
  return(
    tibble(
      shuffle_seq = shuf_seqs[seq_id],
      shuffle_match = shuf_match[seq_id],
      as_tibble(hit)
    )
  )
}

# get target sequences
set.seed(928)

TFval_seqs <- tibble(TF = target_TFs$TF) |>
  group_by(TF) |>
  reframe(
    id = get_TFval_seqs(TF, 55)
  ) |>
  filter(! id %in% id[duplicated(id)]) |>
  group_by(TF) |>
  slice_sample(n = 50) |>
  ungroup() |>
  mutate(
    WT_seq = tamsACR_seqs[id] |>
      as.character()
  ) |>
  rowwise() |>
  mutate(
    motif_shuffle(WT_seq, TF)
  ) |>
  ungroup() |>
  rename('WT_match' = match) |>
  pivot_longer(
    starts_with(c('WT_', 'shuffle_')),
    names_to = c('variant', '.value'),
    names_pattern = '(.*)_(.*)'
  ) |>
  rename('sequence' = seq) |>
  mutate(
    TF = paste0('TF_', TF),
    match = if_else(
      strand == '-',
      match |>
        Biostrings::DNAStringSet() |>
        Biostrings::reverseComplement() |>
        as.character(),
      match
    )
  ) |>
  unite(
    'name',
    id,
    TF,
    variant,
    remove = FALSE
  )

# save data to files
TFval_seqs |>
  select(-sequence) |>
  write_tsv('data/annotation/validation/TFval_annotation.tsv.gz')

### create random sequences for TF validation and in silico evolution ###
# function to create a random DNA sequence
random_seq <- function(len, nuc_freq = c(0.25, 0.25, 0.25, 0.25), seed = NULL) {
  if (! is.null(seed)) {
    set.seed(seed)
  }
  
  sample(c('A', 'C', 'G', 'T'), len, replace = TRUE, prob = nuc_freq) |>
    str_flatten()
}

# function to create random DNA sequences
create_random_seqs <- function(n, len, nuc_freq = c(0.25, 0.25, 0.25, 0.25), TFs = NULL, TF_threshold = 0.001, RE = NULL, factor = 1000) {
  random_seqs <- sapply(seq_len(n * factor), function(x) random_seq(len = len, nuc_freq = nuc_freq, seed = x))
  
  if (! is.null(RE)) {
    random_seqs <- random_seqs[! grepl(RE, random_seqs)]
  }
  
  if (! is.null(TFs)) {
    hits <- scan_sequences(TF_motifs[TFs], Biostrings::DNAStringSet(random_seqs), RC = TRUE, threshold = TF_threshold, nthreads = 0)
    
    random_seqs <- random_seqs[-unique(hits$sequence.i)]
  }
  
  random_seqs <- random_seqs |>
    head(n)
  
  if (length(random_seqs) < n) {
    cli::cli_abort('Did not find enough suitable random sequences. Try increasing `factor`.')
  }
  
  return(random_seqs)
}

# function to calculate background nucleotide frequencies
get_freqs <- function(sp) {
  tamsACR_seqs[str_sub(names(tamsACR_seqs), 1, 2) %in% sp] |>
    Biostrings::letterFrequency(letters = c('A', 'C', 'G', 'T'), as.prob = TRUE) |>
    apply(2, mean)
}

# get background nucleotide frequencies for dicot and monocot sequences
bkg_freqs <- lapply(
  list('dicot' = c('At', 'Sl', 'So'), 'monocot' = c('Zm', 'Sb')), 
  get_freqs
)

# create random sequences with a nucleotide compostion similar to an average dicot or monocot ACR
random_seqs <- tibble(bkg = c('dicot', 'monocot')) |>
  group_by(bkg) |>
  reframe(
    sequence = create_random_seqs(
      n = 30,
      len = 170,
      nuc_freq = bkg_freqs[[bkg]],
      TFs = target_TFs$TF,
      RE = '(GGTCT(C|$))|(GAGAC(C|$))',
      factor = 10000
    )
  ) |>
  mutate(
    id = paste0(bkg, seq_len(n()))
  ) |>
  select(id, sequence) |>
  deframe()

## insert motifs into random sequences
# function to get best motif match
best_match <- function(motif) {
  motif['motif'] |>
    apply(2, function(x) names(x)[which(x == max(x))])|>
    str_flatten()
}

# get best motif matches
motif_strings <- TF_motifs[target_TFs$TF] |>
  sapply(best_match)
names(motif_strings) <- paste0('TF_', target_TFs$TF)

# function to insert motifs in sequnces
insert_motifs <- function(seqs, motif_combi, positions) {
  seq_names <- names(seqs)
  
  motifs <- str_split(motif_combi, '_') |>
    unlist()
  
  for (i in seq_along(motifs)) {
    if (motifs[i] == 'NA') {
      next
    }
    TF_motif <- TF_motifs[[as.integer(motifs[i])]]
    motif_len <- ncol(TF_motif['motif'])
    
    start_pos <- positions[i] - as.integer(motif_len / 2)
    
    str_sub(seqs, start_pos, start_pos + motif_len - 1) <- motif_strings[paste0('TF_', motifs[i])]
  }
  
  if (! is.null(seq_names)) {
    seqs <- tibble(name = seq_names, sequence = seqs)
  }
  
  return(seqs)
}

# insert motif combinations
motif_sets <- list(
  'maize' = target_TFs |> filter(type == 'maize_strong') |> pull(TF),
  'tobacco' = target_TFs |> filter(type == 'tobacco_strong') |> pull(TF),
  'light' = target_TFs |> filter(type %in% c('tobacco_strong', 'light_high')) |> pull(TF),
  'dark' = target_TFs |> filter(type %in% c('tobacco_strong', 'light_low')) |> pull(TF),
  'warm' = target_TFs |> filter(type == 'temp_high') |> pull(TF),
  'cold' = target_TFs |> filter(type == 'temp_low') |> pull(TF)
)

TFins_seqs <- tibble(set = names(motif_sets)) |>
  group_by(set) |>
  reframe( # create all possible combinations of TFs in a given set
    rep(list(c(motif_sets[[set]], NA)), 3) |>
      expand.grid() |>
      as_tibble() |>
      unite(
        'combi',
        everything()
      )
  ) |>
  mutate(
    set = ordered(set, levels = names(motif_sets))
  ) |>
  arrange(set) |>
  filter(combi != 'NA_NA_NA' & ! duplicated(combi)) |>
  group_by(pick(everything())) |>
  reframe(
    insert_motifs(
      seqs = random_seqs,
      motif_combi = combi,
      positions = c(29, 85, 141)
    )
  ) |>
  ungroup()

# identify sequences with a BsaI site after motif insertion
bad_seqs <- TFins_seqs |>
  filter(str_detect(sequence, '(GGTCT(C|$))|(GAGAC(C|$))')) |>
  distinct(name) |>
  pull()

# select final set of random sequences
random_seqs <- random_seqs |>
  enframe(name = 'tmpid', value = 'sequence') |>
  filter(! tmpid %in% bad_seqs) |>
  separate_wider_regex(
    tmpid,
    c(bkg = '.*cot', id = '[0-9]+'),
    cols_remove = FALSE
  ) |>
  group_by(bkg) |>
  slice_head(n = 25) |>
  mutate(
    id = seq_len(n()),
    name = paste0('rnd-', str_sub(bkg, 1, 1), id)
  ) |>
  ungroup()

tmpid_conversion <- random_seqs |>
  select(tmpid, name) |>
  deframe()

random_seqs <- random_seqs |>
  select(name, sequence)

TFins_seqs <- TFins_seqs |>
  mutate(
    name = tmpid_conversion[name]
  ) |>
  drop_na(name) |>
  unite(
    'name',
    name,
    set,
    combi
  )


### combine all TF validation sequences ###
validation_seqs <- TFval_seqs |>
  select(name, sequence) |>
  bind_rows(
    random_seqs,
    TFins_seqs
  )

validation_seqs |>
  mutate(
    name = paste0('>', name)
  ) |>
  write_delim(
    'data/refseq/TFval_sequences.fa',
    delim = '\n',
    col_names = FALSE
  )


### load in silico evolution sequences ###
if (! file.exists('data/CNN_data/evolution_data.tsv.gz')) {
  cli::cli_abort('Cannot find data from in silico evolution. Run `code/CNN/in-silico_evolution.py` to generate it.')
}

evolution_data <- read_tsv('data/CNN_data/evolution_data.tsv.gz') |>
  mutate(
    max_mutations = round * muts_per_round,
    across(c(light, dark, warm, cold, maize), ~ round(.x, 5))
  ) |>
  group_by(id) |>
  mutate(
    mutations = adist(sequence, sequence[round == 0])[,1]
  ) |>
  ungroup()

evo_seqs <- evolution_data |>
  filter(max_mutations %in% c(0, 6, 12)) |>
  group_by(sequence) |>
  summarise(
    across(c(id, light, dark, warm, cold, maize, BsaI, max_mutations, mutations), unique),
    across(c(round, objective, muts_per_round), ~ paste0(.x, collapse = '|'))
  ) |>
  ungroup() |>
  mutate(
    name = paste(id, paste0('evo', muts_per_round), objective, round, sep = '_')
  )

# # to decode the evolution sequence name use:
# evo_seqs |>
#   select(name) |>
#   separate_wider_delim(
#     name,
#     delim = '_',
#     names = c('id', 'orientation', 'muts_per_round', 'objective', 'round'),
#     too_few = 'align_end'
#   ) |>
#   unite(
#     'id',
#     id,
#     orientation,
#     na.rm = TRUE
#   ) |>
#   mutate(
#     muts_per_round = str_replace(muts_per_round, 'evo', '')
#   ) |>
#   separate_longer_delim(
#     c(round, objective, muts_per_round),
#     delim = '|'
#   ) |>
#   mutate(
#     across(c(round, muts_per_round), as.numeric)
#   )


### combine TF validation and in silico evolution sequences ###
validation_seqs <- Biostrings::readDNAStringSet('data/refseq/TFval_sequences.fa') |>
  as.character() |>
  enframe(value = 'sequence')

validation_seqs <- evo_seqs |>
  select(name, sequence) |>
  bind_rows(validation_seqs)


### remove unmodified random sequences (currently present twice; TF insertions and in silico evolution) ###
validation_seqs <- validation_seqs |>
  filter(! str_detect(name, '^rnd-[md][0-9]+$'))


### add enhancer controls ###
ctrl_enh <- tamsACR_seqs[c('35S', '35S(GG)', 'AB80', 'Cab-1', 'rbcS-E9')] |>
  as.character() |>
  enframe(value = 'sequence') |>
  mutate(
    name = paste0(name, '_fwd')
  )

validation_seqs <- validation_seqs |>
  bind_rows(ctrl_enh)


### fill array up with sequences from original dataset ###
array_size <- 35700

n_extra <- array_size - nrow(validation_seqs)

val_ids <- validation_seqs |>
  select(name) |>
  separate_wider_delim(
    name,
    delim = '_',
    names = 'id',
    too_many = 'drop'
  ) |>
  distinct(id) |>
  filter(str_sub(id, 1, 2) != 'rn') |>
  pull(id)

set.seed(928)

extra_seqs <- tamsACR_data |>
  distinct(id, orientation) |>
  filter(! id %in% val_ids) |>
  slice_sample(n = n_extra) |>
  mutate(
    orientation = if_else(# some sequences form a BsaI when used in the rev orientation; use the fwd orientation for them instead
      str_detect(as.character(tamsACR_seqs[id]), '^(AGACC)|(GTCTC)'),
      'fwd',
      orientation
    ),
    sequence = if_else(
      orientation == 'rev',
      as.character(Biostrings::reverseComplement(tamsACR_seqs[id])),
      as.character(tamsACR_seqs[id])
    )
  ) |>
  unite(
    'name',
    id,
    orientation
  )

validation_seqs <- validation_seqs |>
  bind_rows(extra_seqs)


### save sequences to file ###
validation_seqs |>
  mutate(
    name = paste0('>', name)
  ) |>
  write_delim(
    'data/refseq/validation_sequences.fa',
    delim = '\n',
    col_names = FALSE
  )

# sanity checks
validation_seqs |>
  summarise(
    all_seqs = n(),
    correct_length = sum(nchar(sequence) == 170),
    unique_names = n_distinct(name),
    unique_seqs = n_distinct(sequence),
    no_BsaI = sum(! str_detect(sequence, '(GGTCT(C|$))|(GAGAC(C|$))')),
    no_ambiguous_bases = sum(str_detect(sequence, '^[ACGT]+$'))
  ) |>
  pivot_longer(
    everything(),
    names_to = 'test',
    values_to = 'n'
  ) 


### add cloning adapters ###
validation_seqs |>
  mutate(
    sequence = paste0('GAGCGGTCTCCACTC', sequence, 'CTGTAGAGACCGGGC')
  ) |>
  write_delim(
    'data/extra_files/validation_array.fa',
    delim = '\n',
    col_names = FALSE
  )
