### Train a CNN to predict enhancer activity in Plant STARR-seq ###

# Load the required modules:
import os
import tempfile
import glob
import copy
import torch
from torch import nn
import seqdata as sd
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
settings.config_dir = '.'
data_dir = os.path.join('..', '..', 'data', 'modelling_data')
model_dir = os.path.join('..', '..', 'data', 'CNN_model')
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

## Build the model
# Build BiDen model from configuration file:
model_version = 'plantGREP'
model = load_config(f'{model_version}.yaml')

# Train the model:
train.fit_sequence_module(
    model = model,
    sdata = enhSD_train,
    seq_var = 'ohe_seq',
    target_vars = samples,
    in_memory = True,
    epochs = 50,
    batch_size = 128,
    num_workers = dl_workers,
    log_dir = data_dir,
    name = '',
    version = model_version
)

# Load model checkpoint with the lowest validation loss:
model_file = glob.glob(os.path.join(data_dir, model_version, 'checkpoints', '*'))[0]
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
    out_dir = data_dir,
    name = '',
    file_label = model_version,
    num_workers = dl_workers
)

## Trace and save the model 
# Prepare an exemplary input:
test_seq = torch.tensor(enhSD_test.ohe_seq.data[:1], dtype = torch.float32).to('cpu')

# Send model to the CPU and set to evaluation mode:
trimmed_model.to('cpu')
trimmed_model.eval()

# Trace and save model:
if not os.path.exists(model_dir):
    os.makedirs(model_dir)

model_file = os.path.join(model_dir, f'{model_version}.pt')

_ = trimmed_model.to_torchscript(method = 'trace', example_inputs = test_seq, file_path = model_file)

## Save a model compatible with Google Colab
plantGREPcli_dir = os.path.join('..', 'plantGREPcli')

if not os.path.exists(plantGREPcli_dir):
    os.makedirs(plantGREPcli_dir)

with torch.package.PackageExporter(os.path.join(plantGREPcli_dir, 'plantGREP_package.pt')) as exp:
  exp.intern('eugene.**')
  exp.intern('BiDen')
  exp.extern('**')
  exp.save_pickle('plantGREP', 'model.pkl', trimmed_model)