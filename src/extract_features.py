import os
import argparse
import subprocess
import tempfile

import pandas as pd
pd.options.mode.chained_assignment = None #Suppress SettingWithACopy Warning

from utils_io import get_variants, update_progress
from utils_bams import match_variants_to_filenames, get_sam, generate_reads

def main(args):
    validate_arguments(args)
    maf, bam_dirs, mappability, output = args.maf, args.bam_dirs, args.map_bedgraph, args.output

    print("1. Reading Variants from: {}\n".format(maf))
    df = get_variants(maf)
    print("2. Extracting Read Features from: {}".format(bam_dirs))
    df = extract_read_features(df, bam_dirs)
    print("\n3. Extracting Mappability from: {}\n".format(mappability) )
    df = run_mappability(df, mappability)
    print("4. Outputting Result to: {}\n".format(output))
    df.to_csv(output, sep = '\t', index=False)

def add_parser_arguments(parser):
    parser.add_argument(dest='maf', type = str, help = '<Required> maf file containing candidate variants')
    parser.add_argument(dest='map_bedgraph', type = str, help = '<Required> Mappability bedgraph')
    parser.add_argument(dest='output', type = str, help = '<Required> Full path and name of output file')
    parser.add_argument(dest='bam_dirs', nargs="+", type = str, help = '<Required> list of bam directories')

def validate_arguments(args):
    # Checks if input files exist and if output files are in directories that exist and can be written to
    for arg in vars(args):
        print(arg, getattr(args, arg))

    assert os.path.isfile(args.maf)
    assert os.path.isfile(args.map_bedgraph)
    for dir in args.bam_dirs:
        assert os.path.isdir(dir)
    assert os.access(os.path.dirname(args.output), os.W_OK)

def extract_read_features(df, data_dirs):

    def mean(data):
        try: mean = sum(data)/len(data)
        except ZeroDivisionError: mean = float("nan")
        return mean

    def extract_single_variant_features(x):
        # Just to print progress
        next(prog)

        var_sams = [get_sam(data_dir, x['filename']) for data_dir in data_dirs]

        reads = []
        for sam in var_sams: reads+=[r for r in generate_reads(sam, x)]
        alt_reads = [read for read, is_alt in reads if is_alt]

        var_length = mean([len(read.alignment.query_sequence) for read in alt_reads])
        aligned_length = mean([read.alignment.get_cigar_stats()[0][0] for read in alt_reads])
        nm = mean([read.alignment.get_cigar_stats()[0][10] for read in alt_reads])
        softclip = mean([read.alignment.get_cigar_stats()[0][4] for read in alt_reads])
        mapq = mean([read.alignment.mapping_quality for read in alt_reads])
        tlen = mean([abs(read.alignment.template_length) for read in alt_reads])
        tot = len(reads)

        return var_length, aligned_length, nm, softclip, mapq, tlen, tot

    df = match_variants_to_filenames(df, data_dirs)
    prog = update_progress(len(df))
    features = ['f_var_length', 'f_aligned_length', 'f_nm', 'f_softclip', 'f_mapq', 'f_tlen', 'f_tot']
    df[features] = df.apply(extract_single_variant_features, axis=1, result_type="expand")

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
