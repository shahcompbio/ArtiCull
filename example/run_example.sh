#!/bin/bash
set -e

# This script demonstrates how to run ArtiCull classification using example data.

### 1. Extract Features
maf=example/example.maf
features_output=example/articull.features.tsv
bam=example/example.bam
python -m articull extract_features $maf $features_output $bam --cores 1

### 2. Run Classification
features=example/articull.features.tsv
output_dir=example/
model_dir=models/preprint_model/
python -m articull classify $features $output_dir $model_dir --cores 1