 [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# ArtiCull: Variant Call Refinement for scWGS Data

**ArtiCull** (<ins>Arti</ins>fact <ins>Cull</ins>er) is a variant call refinement algorithm for single-cell whole-genome sequencing (scWGS) data. It enables single nucleotide variant (SNV) based analyses in scWGS, including mutational signature and evolutionary analyses.

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
conda env create -f requirements.yml -n articull-env
conda activate articull-env
```

### 2. Setting up Genomic Data Tracks

Use the provided script to download and process the reference genome mappability track:

```bash
bash scripts/setup_mappability_track.bash [optional: output_directory]
```

**Note**: 
- By default, files are saved to the `resources` directory unless an alternative output directory is specified
- The download requires ~1GB of space and expands to ~5GB when uncompressed
- Currently only `hg19`/`GRCh37` is supported. Support for additional reference genomes coming soon
- For other genome versions, please open an issue on [GitHub](https://github.com/shahcompbio/ArtiCull/issues)

## Example Usage

The following example demonstrates how to run ArtiCull. 

### Required Files

- `example/example.maf`: Mutation Annotation Format (MAF) file containing candidate variants
- `example/example.bam`: BAM file with sequencing data
- `example/example.bam.bai`: Index files for BAM

### 0. Follow setup steps and activate conda environment

Followed the directions in [Setup](#setup) to create a conda environment and download required genomic tracks. 

```bash
conda env create -f requirements.yml -n articull-env
conda activate articull-env
bash scripts/setup_mappability_track.bash [output_directory]
```

If you have already done this previously, ensure that the articull-env environment is activated. 

```bash
conda activate articull-env
```

### 1. Extract features

```bash
maf=example/example.maf
features_output=example/articull.features.tsv
bam=example/example.bam
mappability_file=resources/hg19_mappability.bedGraph  # update if you saved to a different location during setup

python src/main.py extract_features $maf $features_output $bam --map_bedgraph $mappability_file --cores 8 
```

### 2. Run classification

```bash
features=example/articull.features.tsv
output_dir=example/
model_dir=models/preprint_model/

python src/main.py classify $features $output_dir $model_dir
```

### Output

The output file, `result.tsv`, is saved in the specified output directory. This is a tab-separated values (TSV) file with the following columns:

| **Column**       | **Description**                                                                      |
|-----------------|------------------------------------------------------------------------------------|
| `chrm`          | Chromosome where the variant is located.                                           |
| `pos`           | Position of the variant on the chromosome.                                         |
| `ref_allele`    | Reference allele at the given position.                                            |
| `alt_allele`    | Alternate allele identified at the given position.                                 |
| `result`        | Classification of the variant (`ARTIFACT`, `PASS`, or `SKIP`).                    |
| `prob_artifact` | Probability that the variant is an artifact (only provided for classified variants). |

Example output:

| **chrm** | **pos**       | **ref_allele** | **alt_allele** | **result** | **prob_artifact**          |
|----------|---------------|----------------|----------------|------------|----------------------------|
| 1        | 201206823     | G              | A              | ARTIFACT   | 0.9729175263160109         |
| 1        | 203994039     | C              | A              | PASS       | 0.01078797299659806        |
| 1        | 201226655     | T              | C              | SKIP       |                            |

- **Classification criteria:**
  - Variants are classified as `ARTIFACT` if `prob_artifact` > 0.5.
  - Variants are classified as `PASS` otherwise.
  - Variants are classified as `SKIP` if no supporting variant reads are found in the BAM file (e.g., due to realignment during variant calling).

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

If you use ArtiCull in your research, please cite the following paper: [update with DOI]

## Issues and Feedback

If you encounter any problems, have questions, or would like to provide feedback, please open an issue on [GitHub](https://github.com/shahcompbio/ArtiCull/issues).
