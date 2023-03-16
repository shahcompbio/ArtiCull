# DLP SNV Filtering

## Dependencies
Dependencies are listed in the `requirements.yml` file. To create and activate a new conda environment with these dependencies,
run
```
conda create -n snv_filter -f requirements.yml
conda activate snv_filter
```
If requirements change, you can update environment using
```
conda activate snv_filter
conda env update -f requirements.yml
```

## Usage:
The software has 5 modes: `extract_features`, `mm_preprocessing`, `mixture_model`, `train_classifier`, `classify`. We include recipes below for how to use these modes for two common use cases -- (1) applying a pre-trained classifier, and (2) training a classifier from new data. Detailed documentation of each mode follows. To run the program in a specific mode, run
```
python src/main.py [mode] [mode specific arguments]
```
To get details on mode specific arguments, see documentation below or run
```
python src/main.py [mode] --help
```

### Recipe: Applying Pre-trained Model
Applying a pre-trained model is straightfoward. First, we extract features for the set of potential variants (listed in the maf file). Then run the classifier.

```
python src/main.py extract_features {maf} {map_bedgraph} {output} {bam_dirs} [--cores ncores]
python src/main.py classify {features} {output_dir} {model_dir}
```

### Recipe: Training New Model
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
python src/main.py mixture_model {maf} {signals_dir} {output} {bam_dirs} [--cores ncores]
```

Finally, create a file listing the features and assignment files for the cohort (as described in mode:classify below), and train the classifer using:

```
python src/main.py mixture model {maf} {signals_dir} {output} {bam_dirs} [--cores ncores]
```

### Mode: Extract features
#### Input

#### Output
Extract features will create a tsv file, containing the following columns

- `chrm`
- `pos`
- `f_var_length`. Average length over all variant reads at position
- `f_aligned_length`. Average aligned length over all variant reads at position
- `f_nm`. Average number of mismatches over all variant reads at position
- `f_softclip`. Average number of softclip bases over all variant reads at position
- `f_mapq`. Average mapping quality over all variant reads at position
- `f_tlen`. Average template length over all variant reads at position
- `f_tot`. Total number of reads at position
- `f_map`. Mappability in 300bp window around position

### Mode: Classify
#### Input

#### Output

### Mode: MM Preprocessing
#### Input

#### Output

### Mode: Mixture Model
#### Input

#### Output

### Mode: Mixture Model
#### Input

#### Output
