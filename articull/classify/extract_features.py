"""
This module provides functions for extracting features from MAF (Mutation Annotation Format) or 
VCF (Variant Calling Format) files.

Functions:
    extract_features(maf: str, bams: str, mappability: str, output: str) -> None:
        Extracts features from given MAF file, BAM files, and mappability data, and outputs the result to a specified file.
"""

import os
import subprocess
import tempfile
import numpy
import random
from scipy.stats import binomtest # type:ignore
import pandas as pd # type:ignore
pd.options.mode.chained_assignment = None #Suppress SettingWithACopy Warning

from articull._utils.io import get_variants
from articull._utils.bams import match_variants_to_filenames,  generate_reads,  get_sam, get_sam_path

import math

def extract_features(maf, bams, mappability, output_prefix):
    """
    Extracts features from given MAF (Mutation Annotation Format) file, BAM files, and mappability data, 
    and outputs the result to a specified file.

    Args:
        maf (str): Path to the MAF file containing variant information.
        bams (str): Path to the directory containing BAM files.
        mappability (str): Path to the directory containing mappability data.
        output (str): Path to the output file where the extracted features will be saved.

    Returns:
        features_file (str): Path to the output file containing the extracted features.
    """
    
    _validate_arguments(maf, bams, mappability, output_prefix)
    random.seed(42)
    mappability = os.path.join(mappability, 'mappability') 
    print("1. Reading Variants from: {}\n".format(maf))
    df = get_variants(maf)
    print("2. Extracting Read Features from: {}".format(bams))
    df = _extract_read_features(df, cell_labels = None, subsample = None, data_dirs = False, filelist = bams)
    print("\n3. Extracting Mappability from: {}\n".format(mappability) )
    df = _run_mappability_by_chrm(df, mappability)
    output = f'{output_prefix}_features.tsv'
    print("4. Outputting Features to: {}\n".format(output))

    df.to_csv(output, sep = '\t', index=False)
    return output


def _validate_arguments(input_file, bams, resources_dir, output):
    """
    Validates the input arguments for the feature extraction process.

    Args:
        input_file (str): Path to the input file.
        bams (list of str): List of paths to BAM files.
        resources_dir (str): Path to the resources directory.
        output (str): Path to the output file.

    Raises:
        AssertionError: If any of the input files do not exist or cannot be read.
        AssertionError: If the resources directory does not exist or does not contain a mappability directory.
        AssertionError: If the output directory does not exist or cannot be written to.
    """
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
def _binomtest_memo(k, n):
    """
    Perform a binomial test and return the logarithm of the p-value.
    Memoizes the results to avoid redundant calculations.

    Args:
        k (int): The number of successes.
        n (int): The number of trials.

    Returns:
        float: The natural logarithm of the p-value from the binomial test.
    """
    return math.log(binomtest(k, n).pvalue)

def _extract_read_features(df, cell_labels=None, subsample=False, data_dirs=False, filelist=False):
    """
    Extracts read features from a given DataFrame of variants.

    Args:
        df (pandas.DataFrame): DataFrame containing variant information.
        cell_labels (list, optional): List of cell labels to filter reads. Defaults to None.
        subsample (bool, optional): Whether to subsample reads. Defaults to False.
        data_dirs (list, optional): List of directories containing SAM files. Defaults to False.
        filelist (list, optional): List of file paths to SAM files. Defaults to False.

    Returns:
        pandas.DataFrame: DataFrame with extracted features.

    Raises:
        Exception: If neither data_dirs nor filelist is provided.

    Features Extracted:
        - f_mean_len: Mean length of alternative reads.
        - f_mean_tlen: Mean template length of alternative reads.
        - f_p_softclip: Proportion of alternative reads with soft clipping.
        - f_mean_mapq: Mean mapping quality of alternative reads.
        - f_mean_mm: Mean number of mismatches in alternative reads.
        - f_p_ins: Proportion of alternative reads with insertions.
        - f_p_del: Proportion of alternative reads with deletions.
        - f_mean_readend: Mean distance to read end in alternative reads.
        - f_std_start: Standard deviation of start positions of aligned blocks in alternative reads.
        - f_std_end: Standard deviation of end positions of aligned blocks in alternative reads.
        - f_p_matemapped: Proportion of alternative reads with mapped mates.
        - f_directionality: Directionality score of alternative reads.
        - f_p_normal: Proportion of normal reads.
        - f_p_XA: Proportion of alternative reads with XA tag.
        - f_mean_XS: Mean XS tag value of alternative reads.
    """
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
            directionality = max(-20, _binomtest_memo(forward, len(alt_reads)))
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

def _run_mappability_by_chrm(df, map_dir):
    """
    Calculate mappability for each chromosome in the given DataFrame.
    This function groups the input DataFrame by chromosome and applies the 
    `run_mappability` function to each group. The mappability data is read 
    from bedGraph files located in the specified directory.

    Args:
        df (pandas.DataFrame): Input DataFrame containing genomic data. 
            It must have a column named 'chrm' indicating the chromosome.
        map_dir (str): Directory path where the bedGraph files for each chromosome are stored.

    Returns:
        pandas.DataFrame: DataFrame with mappability information added for each chromosome.
    """
    chrm_files = {f'chr{chrm}': os.path.join(map_dir, f'chr{chrm}.bedGraph') for chrm in range(1, 23)}  # Assuming chromosomes 1-22
    chrm_files.update({'chrX': os.path.join(map_dir, 'chrX.bedGraph'), 'chrY': os.path.join(map_dir, 'chrY.bedGraph')})

    def get_map_filename(chrm):
        if 'chr' not in chrm:
            chrm = 'chr' + chrm
        return os.path.join(map_dir, f'{chrm}.bedGraph')
    
    df = df.groupby('chrm').parallel_apply(lambda x: _run_mappability(x, get_map_filename(x.iloc[0]['chrm']))).reset_index(drop=True)

    return df

def _get_chrm_bounds(chrm_file):
    """
    Get the highest and lowest positions in a BED file.

    Args:
        chrm_file (str): Path to the BED file.

    Returns:
        tuple: A tuple containing the lowest and highest positions (start, end) in the BED file.
    """
    with open(chrm_file, 'r') as file:
        positions = [(int(line.split('\t')[1]), int(line.split('\t')[2])) for line in file]
    min_pos = min(start for start, _ in positions)
    max_pos = max(end for _, end in positions)
    return min_pos, max_pos

def _check_variant_chrm_bounds(df, chrm_file):
    """
    Validate that the variants in the DataFrame are within the bounds of the corresponding chromosome file.

    Args:
        df (pandas.DataFrame): DataFrame containing genomic data with 'pos' column.
        chrm_file (str): Path to the BED file for the chromosome.

    Returns:
        pandas.Series: A boolean Series where True indicates the variant is within bounds, and False indicates it is out of bounds.
    """
    min_pos, max_pos = _get_chrm_bounds(chrm_file)
    in_bounds = (df['pos'] >= min_pos) & (df['pos'] <= max_pos)
    if not in_bounds.all():
        out_of_bounds_variants = df[~in_bounds]
        for _, row in out_of_bounds_variants.iterrows():
            print(f"Warning: Skipping variant out of chromosome bounds - Chromosome: {row['chrm']}, Position: {row['pos']}. "
                "There may be caused by a reference mismatch between mappability files and variant calls.")

    return in_bounds

def _run_mappability(df, mappability):
    """
    Annotates a DataFrame with mappability scores.
    This function takes a DataFrame containing genomic positions and annotates it with mappability scores
    by creating a BED file, intersecting it with a mappability file using bedtools, and then joining the
    results back to the original DataFrame.
    Variants that are out of bounds are skipped during bedtools operations but included in the final DataFrame with NaN for mappability.

    Args:
        df (pd.DataFrame): DataFrame containing genomic positions with columns ['chrm', 'pos'].
        mappability (str): Path to the mappability file.

    Returns:
        pd.DataFrame: DataFrame with an additional column 'f_map' containing the mappability scores.

    Raises:
        subprocess.CalledProcessError: If bedtools commands fail.

    Notes:
        - The function assumes that the 'chrm' column in the input DataFrame may or may not start with 'chr'.
        - Temporary files are created and deleted within the function.
    """
    def write_bedfile(df, bed):
        bed_df = df[['chrm', 'pos']].copy()
        bed_df['chrm'] = bed_df['chrm'].apply(lambda x: x if str(x).startswith('chr') else 'chr' + str(x))  # N -> chrN
        bed_df['start'] = bed_df['pos'] - 150
        bed_df['end'] = bed_df['pos'] + 150
        del bed_df['pos']
        bed_df.to_csv(bed, index=False, header=False, sep = '\t')
        bed.close()

    def intersect_mappability(mappability, bed, bed_sorted, map):
        sort_bed = f"bedtools sort -i {bed.name}"
        intersect_map = f"bedtools map -a {bed_sorted.name} -b {mappability} -c 4 -o mean -header"
        subprocess.check_call(sort_bed.split(' '), stdout=bed_sorted)
        bed_sorted.close()
        subprocess.check_call(intersect_map.split(' '), stdout=map)
        map.close()

    def join_mappability(df, map):
        map_df = pd.read_table(map, names=['chrm', 'start', 'end', 'f_map'])
        map_df['pos'] = ((map_df['start'] + map_df['end'])/2).astype(int)
        map_df['chrm'] = map_df['chrm'].astype(str)
        if not df['chrm'].iloc[0].startswith('chr'):
            map_df['chrm'] = map_df['chrm'].apply(lambda x: x.replace('chr', ''))  # chrN -> N
        df = df.merge(map_df[['chrm', 'pos', 'f_map']], on=['chrm', 'pos'], how='left')
        return df

    # Identify in-bounds and out-of-bounds variants
    in_bounds = _check_variant_chrm_bounds(df, mappability)
    valid_variants = df[in_bounds].copy()
    invalid_variants = df[~in_bounds].copy()
    invalid_variants['f_map'] = float('nan')  # Assign NaN for mappability

    # Create sorted bed file
    bed, bed_sorted, mapf = [tempfile.NamedTemporaryFile(delete=False, mode='w') for _ in range(3)]

    try:
        # Process only valid variants
        write_bedfile(valid_variants, bed)
        intersect_mappability(mappability, bed, bed_sorted, mapf)
        valid_variants = join_mappability(valid_variants, mapf.name)

    finally:
        # Remove the created temporary files
        os.remove(bed.name)
        os.remove(bed_sorted.name)
        os.remove(mapf.name)

    # Combine valid and invalid variants, preserving the original order
    df = pd.concat([valid_variants, invalid_variants]).sort_index()

    return df