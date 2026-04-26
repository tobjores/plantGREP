### In silico evolution of enhancers

## Import modules
# Import the required modules:
import os
import re
import numpy as np
import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from tqdm import tqdm, trange

# Set torch precision to 'medium' or 'high' to properly utilize the GPU Tensor Cores:
torch.set_float32_matmul_precision('high')

# Define a simple Dataset class for DNA sequences:
class DNAds(Dataset):
    def __init__(self, ohe_seqs):
        self.ohe_seqs = ohe_seqs
    
    def __len__(self):
        return self.ohe_seqs.shape[0]

    def __getitem__(self, idx):
        ohe_seq = self.ohe_seqs[idx]
        seq_tensor = torch.tensor(ohe_seq, dtype = torch.float32)
        return seq_tensor

## Define functions to transform, mutagenize, and score sequences
# Define a function to ohe-hot encode DNA sequences:
def one_hot_encode(seqs):
    if isinstance(seqs, str):
        seqs = [seqs]

    ohe_seqs = []

    for seq in seqs:
        mapping = dict(zip("ACGT", range(4)))    
        seq2 = [mapping[i] for i in seq]
        ohe_seqs.append(np.eye(4, dtype = np.uint8)[seq2])
    
    ohe_seqs = np.stack(ohe_seqs).swapaxes(-1, -2)
    
    return ohe_seqs

# Define a function to decode ohe-hot encoded sequences:
def decode_seqs(ohe_seqs):
    if ohe_seqs.ndim == 2:
        ohe_seqs = np.expand_dims(ohe_seqs, axis = 0)

    seqs = []

    for seq in ohe_seqs:
        seqs.append(''.join(['ACGT'[b] for b in seq.argmax(axis = 0)]))

    if len(seqs) == 1:
        return seqs[0]
    else:
        return seqs

# Define a function to insert mutations in one-hot encoded sequence:
def mutate_seq(seq, mut, win):
    if not 1 <= mut <= 3:
        raise ValueError('`mut` must be between 1 and 3')
    if win < mut:
        raise ValueError('`win` must be greater than or equal to `mut`')

    # singles
    n_choices, seq_len = seq.shape
    idxs = seq.argmax(axis = 0)
    n = seq_len * (n_choices - 1)
    singles = np.tile(seq, (n, 1))
    singles = singles.reshape(n, n_choices, seq_len)
    for k in range(1, n_choices):
        i = np.arange(seq_len) * (n_choices - 1) + (k - 1)
        singles[i, idxs, np.arange(seq_len)] = 0
        singles[i, (idxs + k) % n_choices, np.arange(seq_len)] = 1

    WT = seq.reshape(1, n_choices, seq_len)

    if mut == 1:
        return np.concatenate((WT, singles), axis = 0)
    
    # doubles
    n_singles = singles.shape[0]
    padded_singles = np.pad(singles, ((0, 0), (0, 0), (0, win - 1)))
    doubles = np.repeat(padded_singles, (win - 1) * (n_choices - 1), axis = 0)
    for m in range(1, win):
        idxs_m = np.repeat(np.concatenate((idxs[m:], idxs[:m])), (n_choices - 1))
        for k in range(1, n_choices):
            i = np.arange(n_singles) * (win - 1) * (n_choices - 1) + (n_choices - 1) * (m - 1) + (k - 1)
            doubles[i, idxs_m,  np.repeat(np.arange(seq_len), (n_choices - 1)) + m] = 0
            doubles[i, (idxs_m + k) % n_choices, np.repeat(np.arange(seq_len), (n_choices - 1)) + m] = 1

    if mut == 2:
        valid_seqs = np.all(doubles[:,:,170:] == 0, axis = (1, 2)) # find sequences without out-of-bounds mutations
        doubles = doubles[valid_seqs, :, :seq_len]
        return np.concatenate((WT, singles, doubles), axis = 0)
    
    # triples
    n_triples = (win - 1) - np.repeat([k for k in range(1, win)], (n_choices - 1))
    n_triples = np.tile(n_triples * (n_choices - 1), 170 * (n_choices - 1))
    triples = np.repeat(doubles, n_triples, axis = 0)
    loop_id = 0
    for o in range(2, win):
        offset = loop_id * (n_choices - 1)**2
        for m in range(o, win):
            idxs_m = np.repeat(np.concatenate((idxs[m:], idxs[:m])), (n_choices - 1))
            loop_id += 1
            for j in range(n_choices - 1):
                for k in range(1, n_choices):
                    i = np.arange(n_singles) * ((win - 2) * (win - 1)) // 2 * (n_choices - 1)**2 + offset + (m - o) * (n_choices - 1)  + (k - 1) + j * (win - o) * (n_choices - 1)
                    triples[i, idxs_m,  np.repeat(np.arange(seq_len), (n_choices - 1)) + m] = 0
                    triples[i, (idxs_m + k) % n_choices, np.repeat(np.arange(seq_len), (n_choices - 1)) + m] = 1

    multiples = np.concatenate((doubles, triples), axis = 0)
    valid_seqs = np.all(multiples[:,:,170:] == 0, axis = (1, 2)) # find sequences without out-of-bounds mutations
    multiples = multiples[valid_seqs, :, :seq_len]
    
    return np.concatenate((WT, singles, multiples), axis = 0)

# Define a function to predict enhancer strength:
def predict_enh_strength(sequences, model, batchsize = 1024, pred_names = ['prediction']):
    if next(model.parameters()).is_cuda:
        device = torch.device('cuda')
        dl_kwargs = {'pin_memory': True}
    else:
        device = torch.device('cpu')
        dl_kwargs = {}
    seqDS = DNAds(sequences)
    seqDL= DataLoader(
        dataset = seqDS,
        batch_size = batchsize,
        shuffle = False,
        **dl_kwargs
    )
    predictions = np.empty((0, 5))

    with torch.no_grad():
        for batch in seqDL:
            batch = batch.to(device)
            pred = model(batch).detach().cpu().numpy()
            predictions = np.append(predictions, pred, axis = 0)

    predictions = pd.DataFrame(predictions, columns = pred_names)
    
    return(predictions)

# Define a function to score sequences and pick the best one:
def score_seqs(sequences, model, mode = 'max', columns = None, fn = sum, batchsize = 1024, pred_names = ['predictions']):
    if mode not in ['max', 'min']:
        raise ValueError('`mode` must be "max" (default) or "min"')

    predictions = predict_enh_strength(sequences, model, batchsize = batchsize, pred_names = pred_names)

    columns = predictions.columns if columns is None else columns

    best_id = predictions[columns].apply(fn, axis = 1)
    
    if mode == 'max':
        best_id = best_id.idxmax()
    else:
        best_id = best_id.idxmin()
    
    best_seq = sequences[best_id].copy()
    prediction = predictions.iloc[best_id]

    return(best_seq, prediction)

# Helper function to mutate sequences between 'start' and 'stop':
def mutate_region(sequences, start, stop):
    mut_seqs = []
    
    for sequence in sequences:
        for pos in range(start, stop):
            for base in ['A', 'C', 'G', 'T']:
                if sequence[pos] != base:
                    mut_seq = sequence[:pos] + base + sequence[pos+1:]
                    mut_seqs.append(mut_seq)
                
    return (mut_seqs)

# Define a function to remove restriction enzyme sites from a sequence by mutating the site and picking the best performing sequence:
def fix_REsite(sequence, model, REsite, mode = 'max', columns = None, fn = sum, batchsize = 1024, pred_names = ['predictions']):
    if not isinstance(REsite, re.Pattern):
        REsite = re.compile(REsite)

    matches = REsite.finditer(sequence)

    mut_seqs = [sequence]
    
    # introduce mutations to destroy RE sites
    for match in matches:
        mut_seqs = mutate_region(mut_seqs, match.start(), match.end())
        
    # filter out sequences where a new RE site was introduced
    mut_seqs = [seq for seq in mut_seqs if not REsite.search(seq)]
    
    # one-hot encode mutated sequences
    one_hot_seqs = one_hot_encode(mut_seqs)
    
    # score sequences and pick best
    fixed_seq, prediction = score_seqs(one_hot_seqs, model, mode = mode, columns = columns, fn = fn, batchsize = batchsize, pred_names = pred_names)
    
    # reverse one-hot encoding of fixed sequence
    fixed_seq = decode_seqs(fixed_seq)
        
    return(fixed_seq, prediction)

## Load models and prepare data
# Load the model:
enh_model = torch.jit.load(os.path.join('..', '..', 'data', 'CNN_model', 'BiDen_enh.pt'))

model_outputs = ['light', 'dark', 'warm', 'cold', 'maize']

# Set the model to evaluation mode and send it to the GPU if available:
enh_model.eval()

device = torch.device('cuda') if torch.cuda.is_available() else torch.device('cpu')
device_name = 'GPU' if torch.cuda.is_available() else 'CPU'

enh_model.to(device)

print('Model sent to:', device_name)

# Load validation sequences:
with open(os.path.join('..', '..', 'data', 'refseq', 'TFval_sequences.fa')) as fasta_file:
    fasta_lines = fasta_file.read().splitlines()

fasta_names = [name[1:] for name in fasta_lines[0::2]]
fasta_seqs = fasta_lines[1::2]

# Extract random sequences:
random_seqs = [(name, seq) for name, seq in zip(fasta_names, fasta_seqs) if 'rnd-' in name and not '_' in name]
random_seqs = pd.DataFrame(random_seqs, columns = ['id', 'sequence'])

# Get names of ACR sequences in the validation sequences:
valACRs = [name.split('_')[0] for name in fasta_names if not 'rnd-' in name and 'WT' in name]

# Load CNN data, select test set, and remove sequences already present in the validation sequences:
test_seqs = pd.read_csv(os.path.join('..', '..', 'data', 'modelling_data', 'modelling_data_tamsACR.tsv.gz'), sep = '\t')

test_seqs = test_seqs.query('set == "test"')

test_seqs[['id', 'orientation']] = test_seqs['id'].str.split('_', expand = True)

test_seqs = test_seqs[~test_seqs['id'].isin(valACRs)]

# Select ACR sequences for in silico evolution:
ACR_seqs = test_seqs.sample(n = 100, random_state = 928)

ACR_seqs['id'] = ACR_seqs['id'] + '_' + ACR_seqs['orientation']

ACR_seqs = ACR_seqs[['id', 'sequence']]

# Combine ACR sequences and random sequences:
evo_seqs = pd.concat([ACR_seqs, random_seqs], ignore_index = True)

# One-hot encode start sequences:
start_seqs = one_hot_encode(evo_seqs['sequence'])

# Set up dataframe for results:
evolution_data = pd.DataFrame(data = {
    'round' : 0,
    'objective' : 'start',
    'muts_per_round' : 0,
    'id' : evo_seqs['id'],
    'sequence' : [seq for seq in start_seqs],
})

# Predict enhancer strength of start sequences and combine with results dataframe:
predictions = predict_enh_strength(start_seqs, enh_model, pred_names = model_outputs)

evolution_data = pd.concat([evolution_data, predictions], axis = 1)

# Define evolution objectives:
objectives = {
    'light' : {'mode' : 'max', 'columns' : ['light'], 'fn' : sum},
    'dark' : {'mode' : 'max', 'columns' : ['dark'], 'fn' : sum},
    'warm' : {'mode' : 'max', 'columns' : ['warm'], 'fn' : sum},
    'cold' : {'mode' : 'max', 'columns' : ['cold'], 'fn' : sum},
    'maize' : {'mode' : 'max', 'columns' : ['maize'], 'fn' : sum},
    'constitutive' : {'mode' : 'max', 'columns' : None, 'fn' : sum},
    'maize-specific' : {'mode' : 'max', 'columns' : None, 'fn' : lambda x: x['maize'] - x.drop('maize').mean()},
    'tobacco-specific' : {'mode' : 'min', 'columns' : None, 'fn' : lambda x: x['maize'] - x.drop('maize').mean()},
    'light-specific' : {'mode' : 'max', 'columns' : None, 'fn' : lambda x: x['light'] - x['dark']},
    'dark-specific' : {'mode' : 'min', 'columns' : None, 'fn' : lambda x: x['light'] - x['dark']},
    'warm-specific' : {'mode' : 'max', 'columns' : None, 'fn' : lambda x: x['warm'] - x['cold']},
    'cold-specific' : {'mode' : 'min', 'columns' : None, 'fn' : lambda x: x['warm'] - x['cold']},
}

## Perform the evolution
# Perform iterative evolution of the promoter sequences:
n_rounds = 12

print('Running in silico evolution:')

for obj in tqdm(objectives, desc = 'objectives'):
    for n_muts in trange(1, 4, desc = 'mutations per round', leave = False):
        for this_round in trange(1, 1 + n_rounds // n_muts, desc = 'rounds', leave = False):
            lastround = this_round - 1

            if lastround == 0:
                data = evolution_data[evolution_data['round'] == lastround]
            else:
                data = evolution_data[(evolution_data['round'] == lastround) & (evolution_data['objective'] == obj) & (evolution_data['muts_per_round'] == n_muts)]
                
            if data.empty:
                break
                
            ids = data['id']
            sequences = data['sequence']
            
            for this_id, sequence in tqdm(zip(ids, sequences), total = len(ids), desc = 'sequences', leave = False):
                mut_seqs = mutate_seq(sequence, n_muts, 8)
                opt_seq, prediction = score_seqs(mut_seqs, enh_model, pred_names = model_outputs, **objectives[obj])

                if not(np.array_equal(opt_seq, sequence)):
                    evolution_data.loc[len(evolution_data)] = [this_round, obj, n_muts, this_id, opt_seq, *prediction]

        torch.cuda.empty_cache()

# Reverse one-hot encoding:
ohe_seqs = np.stack(evolution_data['sequence']).astype(np.uint8)
DNA_seqs = decode_seqs(ohe_seqs)

evolution_data['sequence'] = DNA_seqs

# Define BsaI sites as regular expression (also considers sites that would be formed with the cloning scars):
BsaI = re.compile(r'(GGTCT(C|$))|(GAGAC(C|$))')

# Annotate if the sequence contains a restriction enzyme site:
evolution_data['BsaI'] = [bool(BsaI.search(x)) for x in evolution_data['sequence']]

# Select sequences that need to be fixed (round 6 or 12 and contains a restriction enzyme site):
fix_seqs = evolution_data[((evolution_data['round'] * evolution_data['muts_per_round'] == 6) | (evolution_data['round'] * evolution_data['muts_per_round'] == 12)) & (evolution_data['BsaI'])].copy()

# Fix sequences:
for idx in fix_seqs.index:
    new_seq, new_pred = fix_REsite(fix_seqs.loc[idx].sequence, enh_model, BsaI, pred_names = model_outputs, **objectives[fix_seqs.loc[idx].objective])
    fix_seqs.loc[idx, ['sequence', *model_outputs]] = [new_seq, *new_pred]

# Update evolution data with fixed sequences:
evolution_data = fix_seqs.combine_first(evolution_data)

# Convert 'round' column to integer:
evolution_data = evolution_data.astype({'round' : int, 'muts_per_round' : int})

# Save data to a file:
evolution_data.to_csv(os.path.join('..', '..', 'data', 'modelling_data', 'evolution_data.tsv.gz'), sep = '\t', na_rep = 'NA', index = False)