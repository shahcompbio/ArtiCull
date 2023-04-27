import os
import argparse
import subprocess
import tempfile
import numpy



import pandas as pd
pd.options.mode.chained_assignment = None #Suppress SettingWithACopy Warning

from utils_io import get_variants, update_progress
from utils_bams import match_variants_to_filenames,  generate_reads,  get_sam
from scipy.stats import binomtest
import math

def main(args):
    validate_arguments(args)
    maf, bam_dirs, mappability, output, cell_labels, patient_id = args.maf, args.bam_dirs, args.map_bedgraph, args.output, args.cell_labels, args.patient_id

    print("1. Reading Variants from: {}\n".format(maf))
    df = get_variants(maf)
    labels = input_cell_labels(args.cell_labels, args.patient_id)
    print("2. Extracting Read Features from: {}".format(bam_dirs))
    df = extract_read_features(df, bam_dirs, labels)
    print("\n3. Extracting Mappability from: {}\n".format(mappability) )
    df = run_mappability(df, mappability)
    print("4. Outputting Result to: {}\n".format(output))
    df.to_csv(output, sep = '\t', index=False)

def input_cell_labels(filename, patient_id):
    df = pd.read_table(filename, sep='\t')
    df = df[df['patient_id']== patient_id].set_index('cell_id')
    return df[['is_normal_cell', 'is_high_quality_normal_cell']]

def add_parser_arguments(parser):
    parser.add_argument(dest='maf', type = str, help = '<Required> maf file containing candidate variants')
    parser.add_argument(dest='map_bedgraph', type = str, help = '<Required> Mappability bedgraph')
    parser.add_argument(dest='output', type = str, help = '<Required> Full path and name of output file')
    parser.add_argument(dest='cell_labels', type=str, help = '<Required> table indicating which cells are normal')
    parser.add_argument(dest='patient_id', type=str, help = '<Required> patient id, corresponding to patient_id column in cell_label table')
    parser.add_argument(dest='bam_dirs', nargs="+", type = str, help = '<Required> list of bam directories')

def validate_arguments(args):
    # Checks if input files exist and if output files are in directories that exist and can be written to
    for arg in vars(args):
        print(arg, getattr(args, arg))

    assert os.path.isfile(args.maf)
    assert os.path.isfile(args.map_bedgraph)
    assert os.path.isfile(args.cell_labels)
    for dir in args.bam_dirs:
        assert os.path.isdir(dir)
    assert os.access(os.path.dirname(args.output), os.W_OK)

#def get_sam(data_dir, filename):
#    # Could be shared with extract features?
#    return pysam.AlignmentFile(os.path.join(data_dir, filename), "rb")
import functools
@functools.lru_cache # memoizes
def binomtest_memo(k, n):
    return math.log(binomtest(k, n).pvalue)

def extract_read_features(df, data_dirs, cell_labels):

    def mean(data):
        try: mean = sum(data)/len(data)
        except ZeroDivisionError: mean = float("nan")
        return mean

    def extract_single_variant_features(x):
        # Just to print progress
        try:
            var_sams = [get_sam(data_dir, x['filename']) for data_dir in data_dirs]
        except:
            return [None] * len(features)
        reads = []
        for sam in var_sams: reads+=[r for r in generate_reads(sam, x, cell_labels)]
        alt_reads = [read for read, is_alt, is_norm in reads if is_alt]

        var_length = mean([len(read.alignment.query_sequence) for read in alt_reads])
        tlen = min(1000, mean([abs(read.alignment.template_length) for read in alt_reads]))
        num_mm = mean([read.alignment.get_cigar_stats()[0][10] for read in alt_reads])
        softclip = mean([read.alignment.get_cigar_stats()[0][4] > 0 for read in alt_reads])
        mapq = mean([read.alignment.mapping_quality for read in alt_reads])
        num_ins = mean([read.alignment.get_cigar_stats()[0][1] for read in alt_reads ])
        num_del = mean([read.alignment.get_cigar_stats()[0][2] for read in alt_reads ])
        dist_readend = mean([min(read.query_position, len(read.alignment.query_sequence) -  read.query_position)  for read in alt_reads])
        prop_normal = mean([is_alt for read, is_alt, is_norm in reads if is_norm])

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

        return var_length, num_mm, softclip, mapq, tlen, num_ins, num_del, dist_readend, start_variance, end_variance, mate_mapped, directionality, prop_normal#, tot

    df = match_variants_to_filenames(df, data_dirs)
    features = ['f_mean_len', 'f_mean_tlen', 'f_p_softclip', 'f_mean_mapq', 'f_mean_tlen', 'f_p_ins', 'f_p_del', 'f_mean_readend', 'f_std_start', 'f_std_end', 'f_p_matemapped', 'f_directionality', 'f_p_normal'] #, 'f_tot']

    df[features] = df.parallel_apply(lambda x: extract_single_variant_features(x), axis=1, result_type="expand")


    del df['filename']
    return df

def run_mappability(df, mappability):
    """
    Gets the average mappability of the region using bedtools to intersect with
    the mappability bedGraph
    """

    def write_bedfile(df, bed):
        bed_df = df[['chrm', 'pos']]
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
        df = df.merge(map_df[['chrm', 'pos', 'f_map']], on = ['chrm','pos'])
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
