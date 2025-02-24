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

def main(args):
    validate_arguments(args)
    random.seed(42)
    maf, bams, mappability, output = args.input_file, args.bams, args.map_bedgraph, args.output
    
    #cell_labels = None if cell_labels.lower() == 'none' else cell_labels

    print("1. Reading Variants from: {}\n".format(maf))
    df = get_variants(maf)
    
    #labels=None
    #if cell_labels: labels = input_cell_labels(args.cell_labels, args.patient_id)
    print("2. Extracting Read Features from: {}".format(bams))
    #if fullbam:
    df = extract_read_features_new(df, cell_labels = None, subsample = None, data_dirs = False, filelist = bams)
    #else:
    #    df = extract_read_features_new(df, labels, subsample, data_dirs = bam_dirs, filelist = False)
    print("\n3. Extracting Mappability from: {}\n".format(mappability) )
    df = run_mappability(df, mappability)
    print("4. Outputting Result to: {}\n".format(output))
    df.to_csv(output, sep = '\t', index=False)

def input_cell_labels(filename, patient_id):
    df = pd.read_table(filename, sep='\t')
    df = df[df['patient_id']== patient_id].set_index('cell_id')
    return df[['is_normal_cell', 'is_high_quality_normal_cell']]

def add_parser_arguments(parser):
    parser.add_argument(dest='input_file', type = str, help = '<Required> file containing candidate variants')
    parser.add_argument(dest='output', type = str, help = '<Required> Full path and name of output file')
    parser.add_argument(dest='bams', nargs="+", type = str, help = '<Required> list of bam files')

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir))
    default_resources_path = os.path.join(repo_root, 'resources', 'hg19_mappability.bedGraph')
    parser.add_argument('--map_bedgraph', type=str, default=default_resources_path, help='<Optional> Mappability bedgraph (default: resources/hg19_mappability.bedGraph)')

def validate_arguments(args):
    # Checks if input files exist and if output files are in directories that exist and can be written to
    for arg in vars(args):
        print(arg, getattr(args, arg))

    assert os.path.isfile(args.input_file)
    assert os.path.isfile(args.map_bedgraph)
    #assert os.path.isfile(args.cell_labels)
    for bam in args.bams:
        #print(dir)
        assert os.path.isfile(bam)
    assert os.access(os.path.dirname(args.output), os.W_OK)

#def get_sam(data_dir, filename):
#    # Could be shared with extract features?
#    return pysam.AlignmentFile(os.path.join(data_dir, filename), "rb")
import functools
@functools.lru_cache # memoizes
def binomtest_memo(k, n):
    return math.log(binomtest(k, n).pvalue)

def extract_read_features_new(df, cell_labels = None, subsample = False, data_dirs=False, filelist=False):

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

def extract_read_features(df, data_dirs, cell_labels = None, subsample = False, fullbam=False):

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

        try:
            var_sams = [get_sam(data_dir, x['filename']) for data_dir in data_dirs]
        except:
            return [None] * len(features)
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

    df = match_variants_to_filenames(df, data_dirs, fullbam)
    features = ['f_mean_len', 'f_mean_tlen', 'f_p_softclip', 'f_mean_mapq', 'f_mean_mm', 'f_p_ins', 'f_p_del', 'f_mean_readend', 'f_std_start', 'f_std_end', 'f_p_matemapped', 'f_directionality', 'f_p_normal', 'f_p_XA', 'f_mean_XS'] #, 'f_tot']

    df[features] = df.parallel_apply(lambda x: extract_single_variant_features(x), axis=1, result_type="expand")


    del df['filename']
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
