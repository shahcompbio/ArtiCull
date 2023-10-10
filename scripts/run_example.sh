#!/bin/bash

code_dir=$( cd -- "$( dirname -- "$0" )" &> /dev/null && pwd )

# General
map_file=/home/satasg/projects/resources/mappabilityEncodeAlign50mer.bedGraph

# Sample Specific
vcf=/juno/work/shah/isabl_data_lake/analyses/77/41/37741/results/SHAH_H001755_T01_01_DLP01_mutect.vcf.gz
bam=/juno/work/shah/isabl_data_lake/analyses/67/51/36751/results/SHAH_H001755_T01_01_DLP01_all_cells_bulk.bam
output_dir=$code_dir/../example
log_file=$output_dir/classification.log
cores=8
echo $log_file

bash $code_dir/run_classify.sh $vcf $output_dir $bam $log_file $map_file $cores
