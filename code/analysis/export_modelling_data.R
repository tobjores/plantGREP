library(tidyverse)
library(universalmotif)

### split data into training, testing and validation sets ###
# ids:  character; identifiers of sequences to split; sequences that should be kept together (i.e. both orientations of a sequence) should have the same identifier
# p:    numeric of length 3; probabilities for sequences to end up in the train (p[1]), test (p[2]), and validation (p[3]) set
# seed: seed used for random sampling
train_test_val_split <- function(ids, p, seed = NULL) {
  unique_ids <- unique(ids)
  if (! is.null(seed)) {
    set.seed(seed)
  }
  ttv <- sample(c('train', 'test', 'val'), length(unique_ids), replace = TRUE, prob = p)
  names(ttv) <- unique_ids
  return(ttv[ids])
}


### load tamsACR data ###
tamsACR_reps <- readRDS('data/RData/tamsACR_reps.rds')
tamsACR_data <- readRDS('data/RData/tamsACR_main.rds')


### load sequences ###
tamsACR_seqs <- Biostrings::readDNAStringSet('data/refseq/tamsACR_sequences.fa.gz')
tamsACR_seqs_df <- tamsACR_seqs |>
  as.character() |>
  enframe(
    name = 'id',
    value = 'sequence'
  )


### load transcription factor motifs ###
TF_motifs <- read_meme('data/extra_files/TF-clusters.meme')


### select maize replicate data and reshape ###
maize_reps <- tamsACR_reps |>
  filter(condition == 'maize') |>
  select(experiment, id, orientation, enrichment) |>
  mutate(
    experiment = as.integer(as.ordered(experiment))
  ) |>
  group_by(id, orientation) |>
  mutate(
    enrichment_maize = mean(enrichment, na.rm = TRUE)
  ) |>
  pivot_wider(
    names_from = experiment,
    values_from = enrichment,
    names_prefix = 'enrichment_maize_rep'
  )


### prepare data for modelling ###
tamsACR_data |>
  filter(condition != 'maize') |>
  select(condition, id, orientation, enrichment) |>
  mutate(
    condition = str_replace(condition, '35C', 'warm'),
    condition = str_replace(condition, '10C', 'cold')
  ) |>
  pivot_wider(
    names_from = condition,
    values_from = enrichment,
    names_prefix = 'enrichment_'
  ) |>
  drop_na(starts_with('enrichment')) |>
  inner_join(
    maize_reps,
    by = c('id', 'orientation')
  ) |>
  mutate(
    across(
      starts_with('enrichment_maize_rep'),
      ~ coalesce(.x, enrichment_maize)
    )
  ) |>
  inner_join(
    tamsACR_seqs_df,
    by = 'id'
  ) |>
  mutate(
    sequence = if_else(
      orientation == 'rev',
      as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(sequence))),
      sequence
    ),
    set = train_test_val_split(id, p = c(0.8, 0.1, 0.1), seed = 928)
  ) |>
  unite(
    'id',
    id,
    orientation
  ) |>
  write_tsv('data/modelling_data/modelling_data_tamsACR.tsv.gz')


### scan sequences and count transcription factor binding sites ###
TF_hits <- scan_sequences(TF_motifs, tamsACR_seqs, RC = TRUE, threshold = 0.0001, nthreads = 0, no.overlaps = TRUE, no.overlaps.by.strand = FALSE, no.overlaps.strat = 'score') |>
  as_tibble() |>
  count(motif, 'id' = sequence) |>
  mutate(
    motif = str_replace(motif, '-cluster', '')
  ) |>
  pivot_wider(
    names_from = motif,
    values_from = n
  )

TF_hits <- TF_hits |>
  bind_rows(
    tibble(
      id = names(tamsACR_seqs)[! names(tamsACR_seqs) %in% TF_hits$id]
    )
  ) |>
  mutate(
    across(starts_with('TF_'), ~ replace_na(.x, 0))
  ) |>
  arrange(id)

write_tsv(TF_hits, 'data/extra_files/TF_hits.tsv.gz')


### scan sequences for transcription factor binding sites and calculate score ###
# we are searching for motif matches with a logodds score of at least 0
# for each sequence we summarize all matches by using the sum of the logodds scores
TF_motif_scores <- scan_sequences(TF_motifs, tamsACR_seqs, RC = TRUE, threshold = 0, threshold.type = 'logodds.abs', nthreads = 0) |>
  as_tibble() |>
  mutate(
    motif = str_replace(motif, '-cluster', ''),
    rel_score = (score - thresh.score) / (max.score - thresh.score)
  ) |>
  group_by(motif, 'id' = sequence) |>
  summarise(
    score_sum = sum(rel_score)
  ) |>
  ungroup() |>
  pivot_wider(
    names_from = motif,
    values_from = score_sum
  ) |>
  mutate(
    across(-id, ~ replace_na(.x, 0))
  )

write_tsv(TF_motif_scores, 'data/extra_files/TF_motif_scores.tsv.gz')


### generate and save kmer table ###
kmer_table <- Biostrings::oligonucleotideFrequency(tamsACR_seqs, 6) |>
  as_tibble() |>
  mutate(
    id = names(tamsACR_seqs)
  ) |>
  select(id, everything())


bind_rows(
  'fwd' = kmer_table,
  'rev' = kmer_table |>
    rename_with(rev_comp, -id),
  .id = 'orientation'
) |>
  unite(
    'id',
    id,
    orientation
  ) |>
  write_tsv('data/extra_files/kmer_table.tsv.gz')
