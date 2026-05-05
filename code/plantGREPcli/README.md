# Welcome to plantGREPcli

plantGREPcli is a command-line interface to run plantGREP (plant Gene Regulatory Element Predictor), a deep learning model for predicting enhancer strength, developed in [Jores et al., 2026, bioRxiv](https://doi.org/10.64898/2026.04.26.720828).

plantGREPcli offers an easy-to-use interface to:

- predict enhancer strength (`predict`)
- identify the underlying functional sequence motifs (`deeplift`)
- improve enhancer activity (`evolve`)

## Installation

### Prerequisites

- **python** (version >= 3.10.12) with packages:
    - **torch** (version >= 2.4.1)
    - **pytorch-lightning** (version >= 2.4.0)
    - **captum** (version >= 0.5.0)
    - **pandas** (version >= 1.5.3)
    - **seqpro** (verion >= 0.1.11)

### Installation

1. Download the latest plantGREPcli release from GitHub ([direct link](https://github.com/tobjores/plantGREP/releases/latest/download/plantGREPcli.zip)):

    ```
    wget https://github.com/tobjores/plantGREP/releases/latest/download/plantGREPcli.zip
    ```

1. Unzip the downloaded file:

    ```
    unzip plantGREPcli.zip
    ```

1. To test if the installation was successful and all dependencies are fulfilled, run plantGREPcli and display the help message:

    ```
    python plantGREPcli.py --help
    ```

## Usage

### General

Run plantGREPcli by executing the main python script followed by the desired functionality (`predict`, `deeplift`, or `evolve`) and the corresponding arguments, *e.g.*, `python plantGREPcli.py predict ...`

### Arguments

- `--input <file>` (required): path to fasta file with input sequences; use `-` instead of `<file>` to read from standard input

- `--output <file>` (optional): results will be saved to `<file>`; if not given, results are sent to standard output

- `--condition <condition>` (only `deeplift`; required): conditions/species for which to calculate the DeepLIFT attribution scores; one or multiple of `light`, `dark`, `warm`, `cold`, and `maize`

- `--rounds <n>` (only `evolve`; required): number of *in silico* evolution rounds to perform; 10 to 15 usually gives good results

- `--strong <condition>` (only `evolve`; optional): conditions/species for which to inrease enhancer strength during *in silico* evolution

- `--weak <condition>` (only `evolve`; optional): conditions/species for which to decrease enhancer strength during *in silico* evolution

For quick usage instructions use the `--help` flag (*e.g.*, `python plantGREPcli.py predict --help`).

### Usage examples

- Predict enhancer strength of all sequences in file `my_sequences.fa`:

    ```
    python plantGREPcli.py predict --input my_sequences.fa
    ```

- Calculate DeepLIFT attribution scores for the `light` and `dark` condition for all sequences in file `my_sequences.fa`; save results in file `deeplift_scores.tsv`:

    ```
    python plantGREPcli.py deeplift --condition light dark --input my_sequences.fa --output deeplift_scores.tsv
    ```

- Perform 12 rounds of *in silico* evolution for the first sequence in the file `my_sequences.fa`; increase enhancer strength in the light while reducing it in the dark and in maize:

    ```
    head -n 2 my_sequences.fa | python plantGREPcli.py evolve --rounds 12 --strong light --weak dark maize --input -
    ```

## FAQs

1. **What are the sequence requirements?**

    Sequences must be supplied in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format), contain only A, C, G, and T nucleotides, and must be at least 170-bp long (shorter sequnces will be ignored). For sequences longer than 170 bp, you will get an enhancer strength prediction for every possible 170-bp fragment of it.
    
    For the `evolve` function, only a single, 170-bp sequence can be used.

1. **What do the results mean?**

    The `predict` and `evolve` functions return the predicted enhancer strength (log<sub>2</sub>-transformed and normalized to a no-enhancer control) in the indicated condition.

    The `deeplift` function returns the contribution score for each base indicating how much this base contributes to the overall enhancer strength prediction. Positive values are associated with increased enhancer strength, negative values with decreased strength.

1. **What does `light`, `dark`, `warm`, `cold`, and `maize` mean?**

    In our experiments, we measured the enhancer strength of candidate sequences in different assay systems and under different environmental conditions. Experiments were conducted with tobacco plants kept in normal light/dark cycles (`light`), in complete darkness (`dark`), and at elevated (`warm`) or reduced (`cold`) ambient temperature. Additionally, we also conducted experiments in maize protoplasts (`maize`). The plantGREP model was trained to predict enhancer strength in all these conditions/assay systems.
    
    Please read [Jores et al., 2026, bioRxiv](https://doi.org/10.64898/2026.04.26.720828) for more information.

1. **What is *in silico* evolution?**

    *In silico* evolution is a method to generate strong enhancers. For this, all possible single-nucleotide substitution variants of a given sequence are scored with the plantGREP model. The highest-scoring sequence variant is then used as the starting point for the next round, with each iteration introducing a single mutation that should improve enhancer strength.


## Reference

Please cite [Jores et al., 2026, bioRxiv](https://doi.org/10.64898/2026.04.26.720828) when using plantGREPcli. Thank you!