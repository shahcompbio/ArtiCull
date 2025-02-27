import os
import argparse
import subprocess
import tempfile
import numpy
import random

import pandas as pd
pd.options.mode.chained_assignment = None #Suppress SettingWithACopy Warning

from utils_io import get_variants, update_progress
from utils_bams import match_variants_to_filenames,  generate_reads,  get_sam, get_sam_path
from scipy.stats import binomtest
import math

def main(maf, bams, mappability, output):
    
    validate_arguments(maf, bams, mappability, output)
    random.seed(42)
    mappability = os.path.join(mappability, 'mappability') 
    print("1. Reading Variants from: {}\n".format(maf))
    df = get_variants(maf)
    print("2. Extracting Read Features from: {}".format(bams))
    df = extract_read_features(df, cell_labels = None, subsample = None, data_dirs = False, filelist = bams)
    print("\n3. Extracting Mappability from: {}\n".format(mappability) )
    df = run_mappability_by_chrm(df, mappability)
    print("4. Outputting Result to: {}\n".format(output))
    df.to_csv(output, sep = '\t', index=False)

def input_cell_labels(filename, patient_id):
    df = pd.read_table(filename, sep='\t')
    df = df[df['patient_id']== patient_id].set_index('cell_id')
    return df[['is_normal_cell', 'is_high_quality_normal_cell']]

def validate_arguments(input_file, bams, resources_dir, output):
    # Checks if input files exist and if output files are in directories that exist and can be written to
    assert os.path.isfile(input_file), f"Input file {input_file} does not exist."
    assert os.access(input_file, os.R_OK), (
        f"Input file exists, but cannot be read due to permissions: {input_file}"
    )
    assert os.path.isdir(resources_dir), (
        f"Resources directory {resources_dir} does not exist. "
        "See setup instructions for how to download the resources."
    )
    assert os.path.isdir(os.path.join(resources_dir, 'mappability')), (
        f"Resources directory {resources_dir} does not contain a mappability directory. "
        "See setup instructions for how to download the resources."
    )
    output_dir = os.path.dirname(output)
    assert os.path.isdir(output_dir), f"Output directory {output_dir} does not exist."
    assert os.access(output_dir, os.W_OK), (
        f"Output directory exists, but cannot be written to due to permissions: {output_dir}"
    )
    for bam in bams:
        assert os.path.isfile(bam), f"Input bam file {bam} does not exist."
        assert os.access(bam, os.R_OK), (
            f"Input bam file exists, but cannot be read due to permissions: {bam}"
        )


import functools
@functools.lru_cache # memoizes
def binomtest_memo(k, n):
    return math.log(binomtest(k, n).pvalue)

def extract_read_features(df, cell_labels = None, subsample = False, data_dirs=False, filelist=False):

    def mean(data):
        try: mean = sum(data)/len(data)
        except ZeroDivisionError: mean = float("nan")
        return mean

    def extract_single_variant_features(x):

        def get_tag(read, tag):
            try:
                return read.alignment.get_tag(tag)
            except KeyError:
                return False

        if data_dirs:
            try:
                var_sams = [get_sam(data_dir, x['filename']) for data_dir in data_dirs]
            except:
                return [None] * len(features)
        elif filelist:
            var_sams = [get_sam_path(file) for file in filelist]
        else:
            raise Exception("Either a list of files or directories must be provided")

        reads = []
        for sam in var_sams: reads+=[r for r in generate_reads(sam, x, cell_labels)]

        if subsample:
            nreads = random.randint(1, len(reads))
            reads = random.sample(reads, nreads)

        alt_reads = [read for read, is_alt, is_norm in reads if is_alt]

        var_length = mean([len(read.alignment.query_sequence) for read in alt_reads])
        tlen = min(1000, mean([abs(read.alignment.template_length) for read in alt_reads]))
        num_mm = mean([read.alignment.get_cigar_stats()[0][10] for read in alt_reads])
        softclip = mean([read.alignment.get_cigar_stats()[0][4] > 0 for read in alt_reads])
        mapq = mean([read.alignment.mapping_quality for read in alt_reads])
        num_ins = mean([read.alignment.get_cigar_stats()[0][1] > 0 for read in alt_reads ])
        num_del = mean([read.alignment.get_cigar_stats()[0][2] > 0 for read in alt_reads ])
        dist_readend = mean([min(read.query_position, len(read.alignment.query_sequence) -  read.query_position)  for read in alt_reads])
        prop_normal = mean([is_alt for read, is_alt, is_norm in reads if is_norm])
        prop_XA = mean([get_tag(read, 'XA') != False for read in alt_reads])
        mean_XS = mean([get_tag(read, 'XS') for read in alt_reads])


        def get_aligned_block(read):
            blocks = read.alignment.get_blocks()
            return blocks[0][0], blocks[-1][1]

        start_variance = numpy.std([get_aligned_block(read)[0] for read in alt_reads], ddof = 1) # bessel's correction
        end_variance = numpy.std([get_aligned_block(read)[1] for read in alt_reads], ddof = 1)

        mate_mapped = mean([read.alignment.mate_is_mapped for read in alt_reads])
        forward = sum([read.alignment.is_forward for read in alt_reads])
        try:
            directionality = max(-20, binomtest_memo(forward, len(alt_reads)))
        except:
            # If there's only one alt read
            directionality = 0

        #return var_length, num_mm, softclip, mapq, tlen, num_ins, num_del, dist_readend, start_variance, end_variance, mate_mapped, directionality, prop_normal#, tot
        return var_length, tlen, softclip, mapq, num_mm, num_ins, num_del, dist_readend, start_variance, end_variance, mate_mapped, directionality, prop_normal, prop_XA, mean_XS

    if data_dirs: df = match_variants_to_filenames(df, data_dirs)

    features = ['f_mean_len', 'f_mean_tlen', 'f_p_softclip', 'f_mean_mapq', 'f_mean_mm', 'f_p_ins', 'f_p_del', 'f_mean_readend', 'f_std_start', 'f_std_end', 'f_p_matemapped', 'f_directionality', 'f_p_normal', 'f_p_XA', 'f_mean_XS'] #, 'f_tot']

    df[features] = df.parallel_apply(lambda x: extract_single_variant_features(x), axis=1, result_type="expand")

    try: del df['filename']
    except: pass

    return df

def run_mappability_by_chrm(df, map_dir):
    chrm_files = {f'chr{chrm}': os.path.join(map_dir, f'chr{chrm}.bedGraph') for chrm in range(1, 23)}  # Assuming chromosomes 1-22
    chrm_files.update({'chrX': os.path.join(map_dir, 'chrX.bedGraph'), 'chrY': os.path.join(map_dir, 'chrY.bedGraph')})

    def get_map_filename(chrm):
        if 'chr' not in chrm:
            chrm = 'chr' + chrm
        return os.path.join(map_dir, f'{chrm}.bedGraph')
    
    df = df.groupby('chrm').parallel_apply(lambda x: run_mappability(x, get_map_filename(x.iloc[0]['chrm']))).reset_index(drop=True)

    return df

def run_mappability(df, mappability):
    """
    Gets the average mappability of the region using bedtools to intersect with
    the mappability bedGraph
    Note: This is fairly slow in practice. In a future version, I want to split up the processing. 
    Subdivide the mappability file into chunks and run the intersection in parallel. 
    """
    # mappability is in chrN chromosome format. 

    def write_bedfile(df, bed):
        bed_df = df[['chrm', 'pos']]
        bed_df['chrm'] = bed_df['chrm'].apply(lambda x: x if str(x).startswith('chr') else 'chr' + str(x)) # N -> chrN
        bed_df['start'] = bed_df['pos'] - 150
        bed_df['end'] = bed_df['pos'] + 150
        del bed_df['pos']
        bed_df.to_csv(bed, index=False, header=False, sep = '\t')
        bed.close()

    def intersect_mappability(mappability, bed, bed_sorted, map):

        sort_bed = "bedtools sort -i {}".format(bed.name)
        intersect_map = "bedtools map -a {} -b {} -c 4 -o mean -header".format(bed_sorted.name, mappability)
        subprocess.check_call(sort_bed.split(' '), stdout=bed_sorted)
        bed_sorted.close()
        subprocess.check_call(intersect_map.split(' '), stdout=map)
        map.close()

    def join_mappability(df, map):
        map_df = pd.read_table(map, names=['chrm', 'start', 'end', 'f_map'])
        map_df['pos'] = ((map_df['start'] + map_df['end'])/2).astype(int)
        map_df['chrm'] = map_df['chrm'].astype(str)
        if not df['chrm'].iloc[0].startswith('chr'):
            map_df['chrm'] = map_df['chrm'].apply(lambda x: x.replace('chr', '')) # chrN -> N
        df = df.merge(map_df[['chrm', 'pos', 'f_map']], on = ['chrm','pos'], how = 'inner')
        if len(df) != len(map_df):
            print("Warning: {} variants were not correctly matched with mappability.\
                This shouldn't happen. Please raise an issue on github with this error".format(len(map_df)-len(df)))
        return df

    # Create sorted bed file
    bed, bed_sorted, map = [tempfile.NamedTemporaryFile(delete=False, mode='w') for i in range(3)]

    try:
        # Close all files since subprocesses will be writing to them
        write_bedfile(df, bed)
        intersect_mappability(mappability, bed, bed_sorted, map)
        df = join_mappability(df, map.name)

    finally:
        # Remove the created temporary files
        os.remove(bed.name); os.remove(bed_sorted.name); os.remove(map.name)

    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)
