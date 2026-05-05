#!/usr/bin/env python
# plantGREP command-line tool
# Author: Tobias Jores <tobias.jores@hhu.de>

# Import modules:
import os
import sys
import argparse
import warnings
import torch
import pytorch_lightning
import numpy as np
import pandas as pd
import seqpro as sp
from tqdm.auto import tqdm
from captum.attr import DeepLift

# Get path of this script file
script_dir = os.path.dirname(__file__)

# Warning class for invalid sequences:
class SequenceWarning(Warning):
  pass

# Change formatting of SeqeunceWarning messages
original_formatwarning = warnings.formatwarning
def custom_formatwarning(message, category, filename, lineno, line = None):
  if issubclass(category, SequenceWarning):
    return f'\033[36m{category.__name__}\033[0m: {message}\n'
  else:
    return original_formatwarning(message, category, filename, lineno, line)
warnings.formatwarning = custom_formatwarning

# Define model output names
model_outputs = ['light', 'dark', 'warm', 'cold', 'maize']

# Function to parse command-line arguments:
def parse_arguments():
  model_ref = 'Jores et al., 2026, bioRxiv (https://doi.org/10.64898/2026.04.26.720828)'

  parser = argparse.ArgumentParser(
    description = f'plantGREPcli is a command-line tool for using the deep learning model developed in {model_ref}.'
  )
  parser.add_argument(
    '-d',
    '--debug',
    action = 'store_true',
    default = False,
    help = 'Show traceback for errors and warnings'
  )

  io_parser = argparse.ArgumentParser(add_help = False)
  io_parser.add_argument(
    '-i',
    '--input',
    type = str,
    required = True,
    help = 'A fasta file with the input sequences. Use "-" to read from stdin'
  )
  io_parser.add_argument(
    '-o',
    '--output',
    type = str,
    required = False,
    default = sys.stdout,
    help = 'Name of the output file (optional). Will print to stdout if no file name is provided'
  )
  io_parser.add_argument(
    '--no-duplicates',
    action = 'store_true',
    default = False,
    help = argparse.SUPPRESS
  )

  subparsers = parser.add_subparsers(
    required = True,
    help = 'Action to perform',
    dest = 'cmd'
  )

  predict_parser = subparsers.add_parser(
    'predict',
    parents = [io_parser],
    description = f'Enhancer strength prediction using the deep learning model developed in {model_ref}',
    help = 'Predict enhancer strength'
  )

  deeplift_parser = subparsers.add_parser(
    'deeplift',
    parents = [io_parser],
    description = f'DeepLIFT attribution analysis using the deep learning model developed in {model_ref}',
    help = 'Perform DeepLIFT attribution analysis'
  )
  deeplift_parser.add_argument(
    '-c',
    '--condition',
    type = str,
    nargs = '+',
    required = True,
    choices = model_outputs,
    help = 'Experimental condition/assay system for which to run DeepLIFT'
  )

  evolution_parser = subparsers.add_parser(
    'evolve',
    parents = [io_parser],
    description = f'In silico evolution of enhancers using the deep learning model develop in {model_ref}',
    help = 'Perform in silico evolution'
  )
  evolution_parser.add_argument(
    '-r',
    '--rounds',
    type = int,
    required = True,
    help = 'Number of rounds for in silico evolution'
  )
  evolution_parser.add_argument(
    '-s',
    '--strong',
    type = str,
    nargs = '+',
    required = False,
    choices = model_outputs + ['tobacco'],
    default = [],
    help = 'Experimental condition/assay system for which to increase enhancer strength'
  )
  evolution_parser.add_argument(
    '-w',
    '--weak',
    type = str,
    nargs = '+',
    required = False,
    choices = model_outputs + ['tobacco'],
    default = [],
    help = 'Experimental condition/assay system for which to decrease enhancer strength'
  )


  return(parser.parse_args())

# Function to read contents of a fasta file:
def parse_fasta(lines):
  sequences = pd.DataFrame(columns = ['name', 'sequence'])
  seq_id = seq = ''
  for line in lines:
    line = line.strip()
    if len(line) > 0: # ignore empty lines
      if line[0] == '>': # line with sequence name
        if seq != '': # save the last sequence if it is not empty
          sequences.loc[len(sequences)] = [seq_id, seq]
        seq_id = line[1:] # save sequence name
        seq = '' # start new empty sequence
      else: # line with sequence
        seq = seq + line.upper() # append current line to sequence
  if seq != '': # save the last sequence if it is not empty
    sequences.loc[len(sequences)] = [seq_id, seq]
  return(sequences)

# Function to extract sliding windows of a sequence:
def sliding_window_sequences(seq, window_size = 170):
  return([seq[i:i + window_size] for i in range(len(seq) - window_size + 1)])

# Function to load the input sequences:
def read_sequences(*, file = None, string = None, no_dups = False):
  # convert input to sequence dataframe
  if (file is None and string is None) or (file is not None and string is not None):
    raise ValueError('Exactly one of "file" or "string" must be provided')
  if string is not None:
    sequences = parse_fasta(string.split('\n'))
  else:
    # read stdin or open file
    if file == '-':
      sequences = parse_fasta(sys.stdin.readlines())
    else:
      with open(file) as fafile:
        sequences = parse_fasta(fafile.readlines())
  # remove duplicates
  duplicates = sequences.duplicated()
  if any(duplicates):
    warnings.warn(
      f'The following sequences are duplicated and were removed: {", ".join([n if n != "" else "unnamed_seq" for n in sequences["name"][duplicates].unique()])}',
      stacklevel = 2,
      category = SequenceWarning
    )
    sequences = sequences.drop_duplicates()
  # replace emtpy names
  sequences['name'] = sequences['name'].where(sequences['name'] != '', [f'unnamed_seq{i}' for i in range(len(sequences))])
  # check for duplicated names and sequences
  dup_names = sequences['name'].duplicated(keep = False)
  if any(dup_names):
    warnings.warn(
      f'The following sequence names are duplicated and were altered to make them unique: {", ".join(sequences["name"][dup_names].unique())}',
      stacklevel = 2,
      category = SequenceWarning
    )
    name_counts = sequences.groupby('name').cumcount() + 1
    sequences['name'][dup_names] = sequences['name'][dup_names] + ' (' + name_counts[dup_names].astype(str) + ')'
  # abort if there are duplicated sequences (only if `strict` is True)
  if no_dups:
    dup_seqs = sequences['sequence'].duplicated(keep = False)
    if any(dup_seqs):
      raise ValueError(f'Sequences must be unique. Duplicated sequence(s): {", ".join(sequences["name"][dup_seqs])}')
  # filter out short sequences
  short_seqs = [len(s) < 170 for s in sequences['sequence']]
  if any(short_seqs):
    warnings.warn(
      f'The following sequences are shorter than 170 bp and were removed: {", ".join(sequences["name"][short_seqs])}',
      stacklevel = 2,
      category = SequenceWarning
    )
    sequences = sequences[[not b for b in short_seqs]]
  # filter out invalid sequences
  ambig_seqs = sequences['sequence'].str.contains(r'[^ACGTacgt]')
  if any(ambig_seqs):
    warnings.warn(
      f'The following sequences contain non-ACGT characters and were removed: {", ".join(sequences["name"][ambig_seqs])}',
      stacklevel = 2,
      category = SequenceWarning
    )
    sequences = sequences[[not b for b in ambig_seqs]]
  # stop if there are no valid sequnces
  if len(sequences) == 0:
    raise ValueError('No valid sequences supplied')
  # break each sequence into 170-bp sliding windows
  sequences['sequence'] = sequences['sequence'].apply(sliding_window_sequences)
  sequences = sequences.explode('sequence', ignore_index = True)
  sequences['start'] = sequences.groupby('name').cumcount() + 1
  sequences['end'] = sequences['start'] + 169
  # return sequences
  return(sequences)

# Function to load the plantGREP model:
def load_model():
  imp = torch.package.PackageImporter(os.path.join(script_dir, 'plantGREP_package.pt'))
  model = imp.load_pickle('plantGREP', 'model.pkl')
  return(model)

# Function to load the baseline sequences for DeepLIFT:
def load_baselines(device):
  baselines = np.load(os.path.join(script_dir, 'baselines.npy'))
  baselines_tensor = torch.tensor(baselines, dtype = torch.float32).to(device)
  return(baselines_tensor)

# Function to generate single-nucleotide variants:
def generate_variants(sequence):
  seq_variants = [sequence]
  for pos in range(len(sequence)):
    for base in ['A', 'C', 'G', 'T']:
      if sequence[pos] != base:
        seq_var = sequence[:pos] + base + sequence[pos+1:]
        seq_variants.append(seq_var)
  return(seq_variants)

# Function to one-hot encode sequences:
def one_hot_encode(seqs):
  ohe_seqs = sp.ohe(seqs, alphabet = sp.DNA).transpose(0, 2, 1)
  return(ohe_seqs)

# Define a Dataset for DNA sequences:
class DNAds(torch.utils.data.Dataset):
  def __init__(self, sequences):
    self.ohe_seqs = one_hot_encode(sequences)
  def __len__(self):
    return self.ohe_seqs.shape[0]
  def __getitem__(self, idx):
    ohe_seq = self.ohe_seqs[idx]
    seq_tensor = torch.tensor(ohe_seq, dtype = torch.float32)
    return seq_tensor
    
# Function to set up a Dataloader:
def create_dataloader(sequences, batchsize):
  seqDS = DNAds(sequences)
  seqDL = torch.utils.data.DataLoader(
    dataset = seqDS,
    batch_size = batchsize,
    shuffle = False,
    pin_memory = torch.cuda.is_available()
  )
  return(seqDL)

# Function to predict enhancer strength:
def predict(model, sequence_df, device, quiet = False):
  batchsize = 1024
  dataloader = create_dataloader(sequence_df['sequence'].to_list(), batchsize)
  predictions = np.empty((0, 5))
  with torch.no_grad():
    for batch in tqdm(dataloader, desc = 'Predicting enhancer strength', disable = quiet):
      batch = batch.to(device)
      pred = model(batch).detach().cpu().numpy()
      predictions = np.append(predictions, pred, axis = 0)
  predictions = pd.DataFrame(predictions, columns = [f'prediction_{c}' for c in model_outputs])
  predictions = pd.concat([sequence_df, predictions], axis = 1)
  return(predictions)

# Function to run DeepLIFT:
def deeplift(model, sequence_df, device, conditions):
  batchsize = 128
  # Repeat each sequence at least 10 times (to average the DeepLIFT results in the end):
  if len(sequence_df) >= 12:
    rep_factor = 10
  else:
    rep_factor = batchsize // len(sequence_df)
  rep_seqs = [s for s in sequence_df['sequence'] for _ in range(rep_factor)]
  # Fill up last batch with dummy sequences:
  n_seqs = len(rep_seqs)
  missing_seqs = batchsize - n_seqs % batchsize
  if missing_seqs < 128:
    rep_seqs = rep_seqs + ['A' * 170] * missing_seqs
  # Transform sequence dataframe (1 row per nucleotide):
  sequence_df['base'] = sequence_df['sequence'].apply(list)
  sequence_df = sequence_df.explode('base', ignore_index = True)
  sequence_df['position'] = sequence_df['start'] + [i for i in range(1, 171)] * (len(sequence_df) // 170) - 1
  sequence_df = sequence_df.drop(['sequence', 'start', 'end'], axis = 1)
  # Set up dataloader and load baselines:
  dataloader = create_dataloader(rep_seqs, batchsize)
  baselines = load_baselines(device)
  # Run DeepLIFT:
  for condition in conditions:
    output_idx = model_outputs.index(condition)
    attributions = []
    for batch in tqdm(dataloader, desc = f'Running DeepLIFT ({condition})'):
      batch = batch.to(device)
      batch_attributions = DeepLift(model).attribute(batch, baselines = baselines, target = output_idx)
      batch_attributions = batch_attributions * batch
      batch_attributions = batch_attributions.sum(dim = 1)
      attributions.append(batch_attributions.detach().cpu().numpy())
    attributions = np.concatenate(attributions, axis = 0)
    attributions = attributions[:n_seqs].reshape(n_seqs // rep_factor, rep_factor, 170).mean(axis = 1)
    sequence_df[f'deeplift_{condition}'] = attributions.reshape(-1)
  # Summarize results (per-position mean of DeepLIFT score):
  sequence_df = sequence_df.groupby(['name', 'position', 'base'], as_index = False).mean()
  return(sequence_df)

# Function to perform in silico evolution:
def evolve(model, sequence_df, device, rounds, strong, weak):
    if len(sequence_df) != 1:
      raise ValueError('Only a single, 170-bp sequence can be supplied for in silico evolution')
    if len(strong + weak) == 0:
      raise ValueError('At least one evolution objective must be supplied with "--strong <condition>" or "--weak <condition>"')
    elif len(strong + weak) > len(set(strong + weak)):
      raise ValueError('Conditions/assay systems cannot be repeated in "--strong" and "--weak"')
    elif 'tobacco' in strong + weak and len(set(model_outputs[:4]).intersection(strong + weak)) > 0:
      raise ValueError(f'Conditions `{"`, `".join(model_outputs[:4])}` cannot be used together with `tobacco`')
    strong = [f'prediction_{c}' for c in strong]
    weak = [f'prediction_{c}' for c in weak]
    evolution_results = pd.DataFrame()
    best_seq = sequence_df['sequence'][0]
    for round in tqdm(range(1, rounds + 1), desc = 'Performing in silico evolution'):
      sequence_variants = pd.DataFrame({'sequence' : generate_variants(best_seq)})
      predictions = predict(model, sequence_variants, device, quiet = True)
      if round == 1:
        evolution_results = pd.concat([pd.DataFrame({'round' : [0]}), predictions.iloc[[0]]], axis = 1)
      scores = predictions.copy()
      scores['prediction_tobacco'] = scores[[f'prediction_{c}' for c in model_outputs[:4]]].mean(axis = 1)
      scores[weak] = -scores[weak]
      scores = scores[strong + weak].apply(sum, axis = 1)
      best_id = scores.idxmax()
      best_seq = predictions['sequence'][best_id]
      evolution_results = pd.concat([evolution_results, pd.concat([pd.DataFrame({'round' : [int(round)]}), predictions.iloc[[best_id]].reset_index(drop = True)], axis = 1,)])
    return(evolution_results)

# Main:
if __name__ == "__main__":
  # Parse arguments:
  args = parse_arguments()
  # Hide tracebacks from errors and warnings:
  if not args.debug:
    sys.tracebacklimit = 0
    warnings.simplefilter('ignore')
    warnings.simplefilter('always', category = SequenceWarning)
  # Load input sequences and the plantGREP model:
  sequence_df = read_sequences(file = args.input, no_dups = args.no_duplicates)
  model = load_model()
  # set up GPU if available
  if torch.cuda.is_available():
    device = torch.device('cuda')
    model.to(device)
  else:
    device = torch.device('cpu')
  # execute main function
  if args.cmd == 'predict':
    results = predict(model, sequence_df, device)
  elif args.cmd == 'deeplift':
    results = deeplift(model, sequence_df, device, args.condition)
  elif args.cmd == 'evolve':
    results = evolve(model, sequence_df, device, args.rounds, args.strong, args.weak)
  else:
    raise ValueError('Unknown command. Available commands: "predict" or "deeplift"')
  # save results or return as stdout
  results.to_csv(args.output, sep = '\t', index = False)