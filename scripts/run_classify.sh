#!/bin/bash
set -eu


conda activate snv_filter
### Input Files
input=$1
output_dir=$2
bam=$3
log_file=$4
map_file=$5
cores=$6

# Get path of script
code_dir=$( cd -- "$( dirname -- "$0" )" &> /dev/null && pwd )/../src

mkdir -p $output_dir
echo "---- Extracting Features ----" > $log_file
python $code_dir/main.py  extract_features --cores $cores $input $map_file $output_dir/features.tsv None None $bam --fullbam &>> $log_file
echo "---- Classifying ----" >> $log_file
python $code_dir/main.py classify --cores $cores $output_dir/features.tsv $output_dir $code_dir/../models/0817/  &>> $log_file
echo "---- Complete ----" >> $log_file
