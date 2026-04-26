### Train a CNN model ensemble to predict enhancer activity in Plant STARR-seq ###

# Load the required modules:
import os
import tempfile
import glob
import re
import copy
import torch
from torch import nn
import seqdata as sd
import pandas as pd
from eugene import settings
from eugene import preprocess as pp
from eugene import train
from eugene import evaluate
from eugene.models import load_config

# Import the BiDen model architecture and register it with the EUGENe model zoo:
from BiDen import BiDen
import eugene.models.zoo
eugene.models.zoo.BiDen = BiDen

# Set torch precision to 'medium' or 'high' to properly utilize the GPU Tensor Cores:
torch.set_float32_matmul_precision('high')

# Get number of CPUs for dataloader workers:
dl_workers = os.cpu_count() - 1

# Set up directories:
model_version = 'plantGREP'

settings.config_dir = '.'
data_dir = os.path.join('..', '..', 'data', 'modelling_data')
output_dir = os.path.join(data_dir, model_version + '_ensemble')
temp_dir = tempfile.TemporaryDirectory()

## Load and preprocess data
# Load data into SeqData object:
enhSD = sd.read_table(
  name = 'seq',
  out = os.path.join(temp_dir.name, 'enhSD.zarr'),
  tables = os.path.join(data_dir, 'modelling_data_tamsACR.tsv.gz'),
  seq_col = 'sequence',
  batch_size = 4096,
  fixed_length = True,
  overwrite = True
)

# One-hot encode sequences:
pp.ohe_seqs_sdata(enhSD)

# Transpose ohe_seq dimensions to make them compatible with pytorch:
enhSD['ohe_seq'] = enhSD.ohe_seq.transpose('_sequence', '_ohe', 'length')

# Split data into train and test sets:
enhSD_train = enhSD.sel(_sequence=(enhSD['set'] != 'test').compute())
enhSD_test = enhSD.sel(_sequence=(enhSD['set'] == 'test').compute())

# Split train and validation sets:
enhSD_train = enhSD_train.assign(train_val = enhSD_train['set'] == 'train')

# Select samples for training:
samples = ['enrichment_' + s for s in ['light', 'dark', 'warm', 'cold', 'maize', 'maize_rep1', 'maize_rep2', 'maize_rep3']]

## Define a function for model training and evaluation
def train_model(model_version, output_name, samples, seed):
  # Build BiDen model from configuration file:
  model = load_config(f'{model_version}.yaml')

  # # Train the model:
  # train.fit_sequence_module(
  #   model = model,
  #   sdata = enhSD_train,
  #   seq_var = 'ohe_seq',
  #   target_vars = samples,
  #   in_memory = True,
  #   epochs = 50,
  #   batch_size = 128,
  #   num_workers = dl_workers,
  #   log_dir = output_dir,
  #   name = '',
  #   version = output_name,
  #   seed = seed
  # )

  # Load model checkpoint with the lowest validation loss:
  model_file = glob.glob(os.path.join(output_dir, output_name, 'checkpoints', '*'))[0]
  best_model = eugene.models.SequenceModule.load_from_checkpoint(model_file, arch = model.arch)

  # Create a version of the model with only 5 outputs (i.e. without the individual maize reps):
  trimmed_model = copy.deepcopy(best_model)

  c_in = best_model.arch.fcnet.weight.shape[1]
  c_out = 5

  trimmed_model.output_dim = c_out
  trimmed_model.arch.fcnet = nn.Linear(in_features = c_in, out_features = c_out)
  trimmed_model.arch.fcnet.weight = nn.Parameter(best_model.arch.fcnet.weight[:c_out,])
  trimmed_model.arch.fcnet.bias = nn.Parameter(best_model.arch.fcnet.bias[:c_out])

  # Use this model to predict the activity of the held-out test sequences:
  evaluate.predictions_sequence_module(
    trimmed_model,
    sdata = enhSD_test,
    seq_var = 'ohe_seq',
    target_vars = samples[:c_out],
    batch_size = 1024,
    in_memory = True,
    out_dir = output_dir,
    name = '',
    file_label = output_name,
    num_workers = dl_workers
  )

## Train model ensemble
# Set number of iterations:
iterations = 50

# Continue from last iteration if the exceution was interrupted:
pattern = re.compile(r'.*iteration_(\d+)_predictions\.tsv')
prediction_files = glob.glob(os.path.join(output_dir, 'iteration_*_predictions.tsv'))

if len(prediction_files) > 0:
  last_iteration = max([int(pattern.search(fname).group(1)) for fname in prediction_files])
  
  # delete unfinished checkpoints
  for file in glob.glob(os.path.join(output_dir, f'iteration_{last_iteration + 1}', 'checkpoints', '*')):
    os.remove(file)
else:
  last_iteration = 0

# Train models:
for i in range(last_iteration + 1, iterations + 1):
  train_model(
    model_version = model_version,
    output_name = f"iteration_{i}",
    seed = int(i),
    samples = samples
  )

## Combine predictions
# Set up dataframe to receive predictions:
all_predictions = enhSD_test[['id'] + samples[:5]].to_dataframe().reset_index(names = 'seq_id')

# Reshape dataframe to long format:
all_predictions = pd.wide_to_long(
  df = all_predictions,
  stubnames = 'enrichment',
  i = 'id',
  j = 'condition',
  sep = '_',
  suffix = '\w+'
).reset_index().set_index(['seq_id', 'condition'])

# Dictionary to map conditions to the corresponding predictions 
condition_map = {i: samples[i].split('_')[1] for i in range(5)}

# Load predictions and add them to the dataframe:
for i in range(1, iterations + 1):
  # load predictions
  predictions = pd.read_csv(os.path.join(output_dir, f'iteration_{i}_predictions.tsv'), sep = '\t').reset_index(names = 'seq_id')
  
  # reshape predictions
  predictions = pd.wide_to_long(
    df = predictions,
    stubnames = ['predictions', 'target'],
    i = 'seq_id',
    j = 'condition',
    sep = '_',
    suffix = '\d+'
  ).reset_index()
  
  # rename conditions
  predictions['condition'] = predictions['condition'].map(condition_map)

  # select and rename final columns
  predictions = predictions.drop(columns = 'target').rename(columns = {'predictions' : f'prediction_{i}'})

  # merge predictions
  all_predictions = all_predictions.join(predictions.set_index(['seq_id', 'condition']))

# Save combined predictions to file:
all_predictions.reset_index(level = 'condition').to_csv(os.path.join(output_dir, 'ensemble_predictions.tsv.gz'), sep = '\t', na_rep = 'NA', index = False)