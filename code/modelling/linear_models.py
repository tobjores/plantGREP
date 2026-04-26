### Build linear models to predict enhancer strength

## Import modules
import os
import pandas as pd
from sklearn.linear_model import LinearRegression

## Read and preprocess data
# experimental data
exp_data = pd.read_csv(os.path.join('..', '..', 'data', 'modelling_data', 'modelling_data_tamsACR.tsv.gz'), sep = '\t')
exp_data = exp_data.drop(columns = [col for col in exp_data.columns if '_rep' in col])
exp_data[['id', 'orientation']] = exp_data['id'].str.split('_', expand = True)
exp_data = exp_data.set_index(['id', 'orientation'])

conditions = [col.split('_')[1] for col in exp_data.columns if col.startswith('enrichment_')]

# TF motif matches
TF_hits = pd.read_csv(os.path.join('..', '..', 'data', 'extra_files', 'TF_hits.tsv.gz'), sep = '\t')
TF_hits = TF_hits.set_index('id')

# TF binding scores
TF_scores = pd.read_csv(os.path.join('..', '..', 'data', 'extra_files', 'TF_motif_scores.tsv.gz'), sep = '\t')
TF_scores = TF_scores.set_index('id')

# kmer counts
kmers = pd.read_csv(os.path.join('..', '..', 'data', 'extra_files', 'kmer_table.tsv.gz'), sep = '\t')
kmers[['id', 'orientation']] = kmers['id'].str.split('_', expand = True)
kmers = kmers.set_index(['id', 'orientation'])

## Train and evaluate linear models
models = ['TF_hits', 'TF_scores', 'kmers']

for model in models:
  x = globals()[model]

  model_data = exp_data.join(x)

  model_train = model_data[model_data['set'] == 'train']
  model_test = model_data[model_data['set'] == 'test'].copy()

  coef_name = 'TF' if 'TF' in model else model[:-1]
  model_coefs = pd.DataFrame({coef_name : x.columns})

  for condition in conditions:
    model_instance = LinearRegression()
    model_instance.fit(model_train[x.columns], model_train['enrichment_' + condition])

    predictions = model_instance.predict(model_test[x.columns])
    model_test['prediction_' + condition] = predictions

    model_coefs['coefficient_' + condition] = model_instance.coef_
  
  model_test[[col for col in model_test.columns if col.startswith(('enrichment_', 'prediction_'))]].to_csv(os.path.join('..', '..', 'data', 'modelling_data', f'lm_{model}_predictions.tsv.gz'), sep = '\t')
  model_coefs.to_csv(os.path.join('..', '..', 'data', 'modelling_data', f'lm_{model}_coefficients.tsv.gz'), sep = '\t', index = False)