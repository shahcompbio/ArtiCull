#!/bin/bash
set -e

maf=example/example.maf
bam=example/example.bam
output_prefix=example/example
model_dir=models/preprint_model/

python -m articull classify $maf $output_prefix $model_dir $bam --cores 1
