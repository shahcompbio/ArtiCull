import pysam
import os
import pandas as pd
import functools
import sys
import argparse
import subprocess
import tempfile

def parse_arguments():

    def get_arguments():
        parser = argparse.ArgumentParser()
        parser.add_argument(dest='maf', type = str, help = '<Required> maf file containing candidate variants', required=True)
        parser.add_argument(dest='output', type = str, help = '<Required> Full path and name of output file', required=True)
        parser.add_argument(dest='map_bedgraph', type = str, help = '<Required> Mappability bedgraph', required = True)
        parser.add_argument(dest='bam_dirs', nargs="+", type = str, help = '<Required> list of bam directories', required=True)
        args = parser.parse_args()
        return args

    def validate_arguments(args):
        # Checks if input files exist and if output files are in directories that exist and can be written to
        assert os.path.isfile(args['maf'])
        assert os.path.isfile(args['map_bedgraph'])
        assert os.path.isdir(args['bam_dirs'])
        assert os.access(os.access(os.path.dirname(args['output']), os.W_OK))

    args = get_arguments()
    validate_arguments(args)
    return args

def get_variants(filename):
    df = pd.read_table(var_data, skiprows=1)
    df['chrm'] = df['Chromosome'].astype('str')
    df['pos'] = df['Start_Position']
    df = df[(df['Variant_Type'] == 'SNP')&(df['FILTER'] == 'PASS')]
    df['ref_allele'] = df[' Reference_Allele']
    df['alt_allele'] = df['Tumor_Seq_Allele2']

    return df[['chrm', 'pos', 'ref_allele', 'alt_allele']]

def extract_read_features(df, data_dirs):

    def match_variants_to_filenames(df, data_dirs):
        """
        Bam files are assumed to be merged and in region format. Thus, for each variant
        we match it to the appropriate file that contains that locus
        """

        files = [(f, f.strip().split('.')[0].split('-')) for f in os.listdir(data_dirs[0]) if f.endswith('.bam')]
        files = [(filename, (region[0], int(region[1]), int(region[2]))) for filename, region in files]

        def get_file(chrm, pos):
            for filename, region in files:
                if region[0] == chrm and region[1] <= pos and region[2] > pos:
                    return filename
            else:
                raise Error('Variant does not fit in any region')

        df['filename'] = df.apply(lambda x: get_file(x['chrm'], x['pos']), axis=1)
        return df

    @functools.cache # memoizes
    def get_sam(i, filename):
        return pysam.AlignmentFile(os.path.join(data_dirs[i], x['filename']), "rb")

    def generate_reads(samfile, fasta, chrm, pos):
        chrm, pos, ref, alt = x['chrm'], x['pos'], x['ref_allele'], x['alt_allele']

        for pileupcolumn in samfile.pileup(chrm, pos-1, pos+1, min_base_quality=20, min_mapping_quality=0):
            if pileupcolumn.pos != pos-1: continue
            for i, pileupread in enumerate(pileupcolumn.pileups):
                if pileupread.is_del or pileupread.is_refskip: continue
                is_alt = pileupread.alignment.query_sequence[pileupread.query_position] == alt
                yield pileupread, is_alt

    def mean(data):
        try: mean = sum(data)/len(data)
        except ZeroDivisionError: mean = float("nan")
        return mean

    def update_progress():
        prog = 0
        while True:
            yield prog
            prog += 1
            if prog % 500 == 0: print(prog)
    prog = update_progress()

    def extract_single_variant_features(x):
        # Just to print progress
        next(prog)

        var_sams = [get_sam(i, x['filename']) for i in range(len(data_dirs))]

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
    features = ['f_var_length', 'f_aligned_length', 'f_nm', 'f_softclip', 'f_mapq', 'f_tlen', 'f_tot']
    df[features] = df.apply(extract_single_variant_features, axis=1, result_type="expand")

    return df

def run_mappability(df, mappability):
    """
    Gets the average mappability of the region using bedtools to intersect with
    the mappability bedGraph
    """

    def write_bedfile(df, filename):
        bed_df = df[['chrm', 'pos']]
        bed_df['start'] = bed['pos'] - 150
        bed_df['end'] = bed['pos'] + 150
        del bed_df['pos']
        bed_df.to_csv(filename, index=False, header=False, sep = '\t')

    def intersect_mappability(mappability, bed, sorted_bed, map):

        sort_bed = "bedtools sort -i {} > {}".format(bed, sorted_bed)
        intersect_map = "bedtools map -a {} -b {} -c 4 -o mean -header > {}"./
            format(sorted_bed, mappability, map)
        subprocess.call(sort_bed.split(' '))
        subprocess.call(intersect_map.split(' '))

    def join_mappability(df, map):
        map_df = pd.read_table(map, names=['chrm', 'start', 'end', 'f_map'])
        map_df['pos'] = ((map_df['start'] + map_df['end'])/2).astype(int)
        df = df.merge(map_df[['chrm', 'pos', 'f_map']], on = ['chrm','pos'])
        return df

    # Create sorted bed file
    bed, bed_sorted, map = [tempfile.NamedTemporaryFile(delete=False) for i in range(3)]

    try:
        # Close all files since subprocesses will be writing to them
        bed.close(); bed_sorted.close(); map.close()

        write_bedfile(df, bed.name)
        intersect_mappability(mappability, bed.name, sorted_bed.name, map.name)
        df = join_mappability(df, map.name)

    finally:
        # Remove the created temporary files
        os.remove(bed); os.remove(bed_sorted); os.remove(map)

    return df

def main():
    args = parse_arguments()
    df = get_variants(args['maf'])
    df = extract_features(df, args['bam_dirs'])
    df = run_mappability(df, mappability)
    df.to_csv(args['output'], sep = '\t', index=False)

if __name__ == '__main__':
    main()
