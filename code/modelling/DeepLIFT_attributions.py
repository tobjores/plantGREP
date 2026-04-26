### Calculate DeepLIFT attributions with the BiDen enhancer model ###
# Load the required modules:
import os
import glob
import random
import torch
import numpy as np
import pandas as pd
import seqpro as sp
from captum.attr import DeepLift
from eugene import settings
from eugene.models import load_config
from seqexplainer.preprocess._preprocess import dinuc_shuffle_seq

# Set up directories and random seed:
settings.config_dir = '.'
data_dir = os.path.join('..', '..', 'data', 'modelling_data')

random.seed('928')

## Load the model and test sequences
# Import the BiDen model architecture and register it with the EUGENe model zoo:
from BiDen import BiDen
import eugene.models.zoo
eugene.models.zoo.BiDen = BiDen

# Build model and load weights from checkpoint:
model = load_config('BiDen_enh.yaml')
model_file = glob.glob(os.path.join(data_dir, 'BiDen_enh', 'checkpoints', '*'))[0]
model = eugene.models.SequenceModule.load_from_checkpoint(model_file, arch = model.arch)

model_outputs = ['light', 'dark', 'warm', 'cold', 'maize']

# Set the model to evaluation mode and send it to the GPU if available:
model.eval()

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
device_name = 'GPU' if torch.cuda.is_available() else 'CPU'

model.to(device)

print('Model sent to:', device_name)

# Load sequences of test set:
modelling_data = pd.read_csv(os.path.join(data_dir, 'modelling_data_tamsACR.tsv.gz'), sep = '\t')
test_data = modelling_data[modelling_data['set'] == 'test'].reset_index(drop = True)
test_seqs = test_data.sequence

# Create 128 dinucleotide-shuffled sequences as baseline:
random_seqs = random.choices(list(test_seqs), k = 128)

shuffled_seqs = []
for i in random_seqs:
    shuffled_seqs.append(dinuc_shuffle_seq(i, rng = np.random.RandomState(928)))

# One-hot encode sequences:
ohe_test_seqs = sp.ohe(test_seqs, alphabet = sp.DNA).transpose(0, 2, 1)
ohe_shuffled_seqs = sp.ohe(shuffled_seqs, alphabet = sp.DNA).transpose(0, 2, 1)

# Save baseline sequences:
np.save(os.path.join('..', 'BiDenCLI', 'baselines'), ohe_shuffled_seqs)

## Calculate DeepLIFT attributions
# Preparations:
dl = DeepLift(model)

baselines = torch.tensor(ohe_shuffled_seqs).to('cuda')

DeepLIFT_attributions = {'ohe_seqs': ohe_test_seqs.transpose(0, 2, 1)}

batch_size = 128

print('Running DeepLIFT:')

# Calculate attributions for each condition:
for i, condition in enumerate(model_outputs):
  print(f'Condition: {condition}')

  condition_attributions = []

  for batch in range(0, len(ohe_test_seqs), batch_size):
    for j in range(10):
      batch_seqs = torch.tensor(np.roll(ohe_test_seqs[batch:batch + batch_size], j, axis = 0), dtype = torch.float32).to('cuda')

      if j == 0:
        attributions = dl.attribute(batch_seqs, baselines = baselines[:len(batch_seqs)], target = i)
      else:
        attributions += torch.roll(dl.attribute(batch_seqs, baselines = baselines[:len(batch_seqs)], target = i), -j, dims = 0)

    condition_attributions.append(attributions / 10)

  condition_attributions = torch.cat(condition_attributions, dim = 0).detach().cpu().numpy()
  DeepLIFT_attributions[f'DeepLIFT_{condition}'] = condition_attributions.transpose(0, 2, 1)

# Save test sequences and DeepLIFT attributions to files:
deeplift_dir = os.path.join(data_dir, 'DeepLIFT')

if not os.path.isdir(deeplift_dir):
  os.mkdir(deeplift_dir)
        
for file in DeepLIFT_attributions:
  np.savez_compressed(os.path.join(deeplift_dir, file), DeepLIFT_attributions[file])

## Export DeepLIFT attributions for an example sequence (At-12806_fwd)
# Get index:
idx = test_data.index[test_data['id'] == 'At-12806_fwd'][0]

# Extract DeepLIFT attributions for all model outputs:
DeepLIFT_scores = pd.DataFrame()

for condition in model_outputs:
  contrib_scores = DeepLIFT_attributions['ohe_seqs'][idx] * DeepLIFT_attributions['DeepLIFT_' + condition][idx]
  contrib_scores = pd.DataFrame(
    data = {
      'condition' : condition,
      'pos' : [i + 1 for i in range(len(contrib_scores))],
      **{b : contrib_scores[:, i] for i, b in enumerate(['A', 'C', 'G', 'T'])}
    }
  )

  DeepLIFT_scores = pd.concat([DeepLIFT_scores, contrib_scores])

# Save to file:
DeepLIFT_scores.to_csv(os.path.join(deeplift_dir, 'DeepLIFT_example.tsv'), sep = '\t', na_rep = 'NA', index = False)