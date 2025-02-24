# ArtiCull: Variant Call Refinement for scWGS Data

**ArtiCull** is a variant call refinement algorithm for single-cell whole-genome sequencing (scWGS) data. It enables single nucleotide variant (SNV) based analyses in scWGS, including mutational signature and evolutionary analyses.

## Table of Contents

- [Setup](#setup)
- [Example Usage](#example-usage)
- [Classification](#classification)
  - [Extracting Features](#extracting-features)
  - [Running Classification](#running-classification)
- [Citing ArtiCull](#citing-articull)

## Setup

### 1. Setting up the Conda Environment

The recommended way to set up ArtiCull is using conda. Alternatively, you can ensure all packages listed in `requirements.txt` are available.

First, ensure you have [conda](https://docs.conda.io/en/latest/) installed on your machine. Then run:

```bash
conda env create -f requirements.txt -n articull-env
conda activate articull-env
```

### 2. Setting up Genomic Data Tracks

Use the provided script to download and process the reference genome mappability track:

```bash
bash scripts/setup_mappability_track.bash <genome_name> [output_directory]
```

**Note**: 
- The download requires ~1GB of space and expands to ~5GB when uncompressed
- Files are saved to the `resources` directory by default unless an alternative output directory is specified
- Currently only `hg19`/`GRCh37` is supported (additional reference genomes coming soon)
- For other genome versions, please open an issue on GitHub

## Example Usage

The following example demonstrates how to run ArtiCull. A complete example script is also available at `examples/run_example.sh`.

### Required Files

- `example.maf`: Mutation Annotation Format (MAF) file containing candidate variants
- `example.bam`: BAM file(s) with sequencing data
- `example.bam.bai`: Index files for BAM
- `resources/hg19_mappability.bedGraph`: Mappability track file
  - See [Setup](#setup) section if this file is missing
  - See [Extracting Features](#extracting-features) if you saved this file to a custom location

### 0. Activate Environment

First, ensure you're using the conda environment created during setup:

```bash
conda activate articull-env
```

### 1. Extract Features

```bash
maf=example/example.maf
features_output=example/articull.features.tsv
bam=example/example.bam

python src/main.py extract_features $maf $features_output $bam --cores 8
```

### 2. Run Classification

```bash
features=example/articull.features.tsv
output_dir=example/
model_dir=models/preprint_model/

python src/main.py classify $features $output_dir $model_dir
```

## Classification

Classification involves two steps using a pretrained model (available in the `models` directory). For training new models, see [model_training.md](doc/model_training.md).

### Extracting Features

The `extract_features` command processes input variants and generates features required for classification.

#### Usage

```bash
python src/main.py extract_features <input_file> <output> <bams> [--map_bedgraph <bedgraph>] [--cores <ncores>]
```

#### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `<input_file>` | Yes | Candidate variants file (MAF or VCF format) |
| `<output>` | Yes | Output path for extracted features (tab-separated format) |
| `<bams>` | Yes | One or more BAM files containing sequencing data |
| `--map_bedgraph` | No | Path to mappability bedgraph file (default: `resources/hg19_mappability.bedGraph`) |
| `--cores` | No | Number of CPU cores for parallel processing |

### Running Classification

The `classify` command processes extracted features using a pre-trained model.

#### Usage

```bash
python src/main.py classify <features> <output_dir> <model_dir> [--chunksize <n>] [--cores <ncores>]
```

#### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `<features>` | Yes | Input file containing extracted features |
| `<output_dir>` | Yes | Directory for classification results |
| `<model_dir>` | Yes | Directory containing pre-trained model (`model.pkl`) and scaler (`scaler.pkl`) |
| `--chunksize` | No | Rows per worker for parallel processing (default: 5000) |
| `--cores` | No | Number of CPU cores for parallel processing |

## Citing ArtiCull

If you use ArtiCull in a manuscript, please cite the following paper: [update with DOI]
