### Identify motifs associated with enhancer strength using TF-MoDISco ###
# Load the required modules:
import os
import modiscolite
import h5py
import numpy as np
import pandas as pd

# Set up:
data_dir = os.path.join('..', '..', 'data', 'modelling_data')
out_dir = os.path.join(data_dir, 'TF-MoDISco')

if not os.path.isdir(out_dir):
  os.mkdir(out_dir)

conditions = ['light', 'dark', 'warm', 'cold', 'maize']

# Load one-hot encoded sequences:
ohe_seqs = np.load(os.path.join(data_dir, 'DeepLIFT', 'ohe_seqs.npz'))['arr_0']

# Run TF-MoDISco for all conditions:
for condition in conditions:
  print(f'Running TF-MoDISo for condition "{condition}":')

  DeepLIFT_data = np.load(os.path.join(data_dir, 'DeepLIFT', f'DeepLIFT_{condition}.npz'))['arr_0']

  pos_patterns, neg_patterns = modiscolite.tfmodisco.TFMoDISco(
    hypothetical_contribs = DeepLIFT_data, 
    one_hot = ohe_seqs,
    verbose = True,
    sliding_window_size = 15, 
    flank_size = 10,
    target_seqlet_fdr = 0.01
  )

  modiscolite.io.save_hdf5(
    filename = os.path.join(out_dir, f'TF-MoDISco_patterns_{condition}.h5'),
    pos_patterns = pos_patterns,
    neg_patterns = neg_patterns,
    window_size = 15
  )

# Save pattern data as .tsv files:
all_patterns = pd.DataFrame()

for condition in conditions:
  with h5py.File(os.path.join(out_dir, f'TF-MoDISco_patterns_{condition}.h5')) as modisco_file:
    for pattern in ['pos', 'neg']:
      if pattern + '_patterns' in modisco_file:
        for name, datasets in modisco_file[pattern + '_patterns'].items():
          pattern_df = pd.DataFrame(
            data = {
              'condition' : condition,
              'pattern' : pattern + '_' + name,
              'pos' : [i + 1 for i in range(len(datasets['sequence']))],
              **{'ppm_' + b : datasets['sequence'][:, i] for i, b in enumerate(['A', 'C', 'G', 'T'])},
              **{'cwm_' + b : datasets['contrib_scores'][:, i] for i, b in enumerate(['A', 'C', 'G', 'T'])}
            }
          )
        
          all_patterns = pd.concat([all_patterns, pattern_df])

all_patterns.to_csv(os.path.join(out_dir, 'TF-MoDISco_patterns.tsv.gz'), sep = '\t', na_rep = 'NA', index = False)