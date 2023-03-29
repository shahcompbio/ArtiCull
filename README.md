# DLP SNV Filtering
## Table of Contents
- [Setup](https://github.mskcc.org/SatasG/dlp-snv-classifier#setup)
- [Usage Overview](https://github.mskcc.org/SatasG/dlp-snv-classifier#usage-overview)
- [Recipes](https://github.mskcc.org/SatasG/dlp-snv-classifier#recipes)
  - [Applying pre-trained model](https://github.mskcc.org/SatasG/dlp-snv-classifier#applying-pre-trained-model)
  - [Training new model](https://github.mskcc.org/SatasG/dlp-snv-classifier#applying-pre-trained-model)
- [Modes](https://github.mskcc.org/SatasG/dlp-snv-classifier#modes)
  - [extract_features](https://github.mskcc.org/SatasG/dlp-snv-classifier#extract-features)
  - [classify](https://github.mskcc.org/SatasG/dlp-snv-classifier#classify)
  - [mm_preprocessing](https://github.mskcc.org/SatasG/dlp-snv-classifier#mm-preprocessing)
  - [mixture_model](https://github.mskcc.org/SatasG/dlp-snv-classifier#mixture-model)
  - [train](https://github.mskcc.org/SatasG/dlp-snv-classifier#mm-preprocessing)


## Setup
Dependencies are listed in the `requirements.yml` file. To create and activate a new conda environment with these dependencies,
run
```
conda env create -n snv_filter -f requirements.yml
conda activate snv_filter
```
If requirements change, you can update environment using
```
conda activate snv_filter
conda env update -f requirements.yml
```

## Usage Overview:
The software has 5 modes: `extract_features`, `mm_preprocessing`, `mixture_model`, `train_classifier`, `classify`. We include recipes below for how to use these modes for two common use cases -- (1) applying a pre-trained classifier, and (2) training a classifier from new data. Detailed documentation of each mode follows. To run the program in a specific mode, run
```
python src/main.py [mode] [mode specific arguments]
```
To get details on mode specific arguments, see documentation below or run
```
python src/main.py [mode] --help
```

## Recipes
### Applying Pre-trained Model

To run a pretrained model (either output from `train` mode or contained in `models/`): First, extract features for the set of potential variants (listed in the maf file). Then, run the classifier.

```
python src/main.py extract_features {maf} {map_bedgraph} {output} {bam_dirs} [--cores ncores]
python src/main.py classify {features} {output_dir} {model_dir}
```

### Training New Model
To train a new model, for each individual/sample, we need a features file (output from `extract_features`) and an assignments file (output from `mixture_model`).

Generate the feature file with
```
python src/main.py extract_features {maf} {map_bedgraph} {output} {bam_dirs} [--cores ncores]
```

Running the mixture model requires a manual inspection step following preprocessing to select a pair of clones. First run preprocessing as:

```
python src/main.py mm_preprocessing {maf} {signals_dir} {output_dir} {bam_dirs} [--cores ncores]
```

`mm_preprocessing` will produce a plot named `clone_ccfs.pdf`. Each row and column corresponds to one clone in the file. Select a pair of clones (if it exists) which have a distinct clonal and shared subclonal distribution. Then run the mixture_model as:

```
python src/main.py mixture_model {input_file} {output_dir} {clone1} {clone2} [--cores ncores]
```

Finally, create a file listing the features and assignment files for the cohort (as described in mode:classify below), and train the classifer using:

```
python src/main.py train {filelist} {output_dir}
```

## Modes
### Extract features
#### Usage

```
python src/main.py extract_features {maf} {map_bedgraph} {output} {bam_dirs} [--cores ncores]
```

#### Positional arguments
- `maf`: maf format file containing candidate variants
- `mappability_bedgraph`: bedgraph containing mappability information across the genome
- `output`: full path and name of output file
- `bam_dirs`: list of bam directories (as output by isabl app `merge_bams`)

#### Optional arguments
- `--cores`: Number of workers to use for parallelization. By default will use all available cores


#### Output
Extract features will create a tsv file to location specified by `output`, containing the following columns

- `chrm`
- `pos`
- `f_var_length`. Average length over all variant reads at position
- `f_aligned_length`. Average aligned length over all variant reads at position
- `f_nm`. Average number of mismatches over all variant reads at position
- `f_softclip`. Average number of softclip bases over all variant reads at position
- `f_mapq`. Average mapping quality over all variant reads at position
- `f_tlen`. Average template length over all variant reads at position
- `f_map`. Mappability in 300bp window around position

### Classify
#### Usage

```
python src/main.py classify {features} {output_dir} {model_dir}
```

#### Positional arguments
- `features`: Input file containing variant features (as output by extract_features). All feature columns must start with 'f_'
- `output_dir`: Output directory
- `model_dir`: Directory containing model.pkl and scaler.plk (in `models/` or as output by train_classifier)

#### Output
`classify` will output a file called `result.tsv` in `output_dir` containing the following columns

- `chrm`
- `pos`
- `result`: variant classification, either `PASS` or `ARTIFACT`
- `probs`: model probability of each classification

### MM Preprocessing
#### Usage
```
python src/main.py mm_preprocessing {maf} {signals_dir} {output_dir} {bam_dirs} [--cores ncores]
```

#### Positional arguments
- `maf`: maf format file containing candidate variants
- `signals_dir`: The output directory for signals
- `output_dir`: Output directory
- `bam_dirs`: list of bam directories (as output by isabl app `merge_bams`)

#### Optional arguments
- `--cores`: Number of workers to use for parallelization. By default will use all available cores

#### Output
`mm_preprocessing` will output two files.

- `{output_dir}/var_counts.tsv`: Read counts for variant and reference alleles per clone
- `{output_dir}/clone_ccfs.pdf`: Pairwise plots of clone ccfs.

### Mixture Model
#### Usage
```
python src/main.py mixture_model {input_file} {output_dir} {clone1} {clone2} [--cores ncores]
```

#### Positional arguments
- `input_file`: `var_counts.tsv` file from `mm_preprocessing`
- `output_dir`: output_directory
- `clone1`: selected clone for mixture model
- `clone2`: selected clone for mixture model

#### Optional arguments
- `--cores`: Number of workers to use for parallelization. By default will use all available cores

#### Output
`mixture_model` will output two files

- `{output_dir}/assignments.tsv`: tsv file with variant labels from mixture model
- `{output_dir}/ccfs_labeled.pdf`: labeled pairwise plot for selected ccfs with the inferred classification

### Train classifier
#### Usage
```
python src/main.py train {filelist} {output_dir}
```

#### Positional arguments
- `filelist`: A list of files containing training data. Each line corresponds to one sample or individual, and should contain the features file (output by `extract_features`) and assignment file (output by `mixture_model`) separate by a tab. This file has no header
- `output_dir`: Output directory

#### Output
`train` will output a pickled model and data scaler `model.pkl` and `scaler.pkl` to `output_dir`
