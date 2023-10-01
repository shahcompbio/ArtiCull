#!/bin/bash

code_dir=$( cd -- "$( dirname -- "$0" )" &> /dev/null && pwd )/../src

# General
map_file=/home/satasg/projects/resources/mappabilityEncodeAlign50mer.bedGraph

# Sample Specific
maf=/home/satasg/projects/2023/snv_filter/data/mutsigs/mafs/OV2295/OV2295_SA1090_A96213A_L5/mutect.maf
bam=/juno/work/shah/isabl_data_lake/analyses/67/51/36751/results/SHAH_H001755_T01_01_DLP01_all_cells_bulk.bam
output_dir=$code_dir/../example
log_file=$output_dir/classification.log

bash run_classify.sh $maf $output_dir $bam $log_file $map_file