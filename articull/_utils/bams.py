"""
This module provides utility functions for working with BAM files and genomic variants.
It includes functionality for matching variants to BAM filenames, determining if reads are from high-quality normal cells,
opening SAM/BAM files, and generating reads at specific positions.

Functions:
    match_variants_to_filenames(df, data_dirs):
        Matches genomic variants to BAM filenames based on their chromosomal location.

    is_normal(pileupread, labels=None):
        Determines if a given pileup read is from a high-quality normal cell.

    get_sam_path(path):
        Opens a SAM/BAM file using pysam and returns the AlignmentFile object.

    get_sam(data_dir, filename):
        Opens a SAM/BAM file using pysam.AlignmentFile.

    generate_reads(samfile, x, labels=None, min_base_quality=20, min_mapping_quality=0):
        Generates reads at a specified position and reports whether they match the variant allele or not.
"""

import functools
import os
import pysam # type: ignore

def match_variants_to_filenames(df, data_dirs):
    """
    Matches genomic variants to BAM filenames based on their chromosomal location.
    This function assumes that BAM files are merged and in region format. For each variant
    in the provided DataFrame, it matches the variant to the appropriate BAM file that 
    contains the corresponding locus.

    Args:
        df (pandas.DataFrame): A DataFrame containing genomic variants with columns 'chrm' (chromosome) 
            and 'pos' (position).
        data_dirs (list): A list of directories where BAM files are stored. Only the first directory 
            in the list is used.

    Returns:
        pandas.DataFrame: The input DataFrame with an additional column 'filename' that contains the 
            matched BAM filename for each variant. If no matching file is found, the filename will be None.

    Note:
        This function may not be used anymore since switching to merged BAMs, but it should be checked.
    """

    def get_region(f):
        v = f.strip().split('.')[0].split('_')[-1].split('-')
        return (v[0], int(v[1]), int(v[2]))

    files = [(f, get_region(f)) for f in os.listdir(data_dirs[0]) if f.endswith('.bam')]

    def get_file(chrm, pos):
        for filename, region in files:
            if region[0] == chrm and region[1] <= pos and region[2] > pos:
                return "{}-{}-{}.bam".format(region[0], str(region[1]), str(region[2]))
        return None

    try:
        temp = df.apply(lambda x: get_file(x['chrm'], x['pos']), axis=1)
        df['filename'] = temp
    except:
        print(temp)
    return df

def is_normal(pileupread, labels=None):
    """
    Determines if a given pileup read is from a high-quality normal cell.

    Args:
        pileupread (pysam.PileupRead): The pileup read to be evaluated.
        labels (pandas.DataFrame, optional): A DataFrame containing cell labels with 
            'is_high_quality_normal_cell' column. Defaults to None.

    Returns:
        bool: True if the cell is a high-quality normal cell, False otherwise.
        None: If labels are not provided.
    """
    if labels is None: return None

    RG = pileupread.alignment.get_tag('RG')
    cell_id = '_'.join(RG.split('_')[1:-2])
    try:
        is_hq_normal = labels.loc[cell_id].is_high_quality_normal_cell
    except:
        is_hq_normal = False
    return is_hq_normal

@functools.lru_cache
def get_sam_path(path):
    """
    Opens a SAM/BAM file using pysam and returns the AlignmentFile object.

    Args:
        path (str): The file path to the SAM/BAM file.

    Returns:
        pysam.AlignmentFile: An object representing the opened SAM/BAM file.

    Note:
        It should be checked if `get_sam` or `get_sam_path` is currently used and clean up the other.
    """
    pysam.set_verbosity(0)
    return pysam.AlignmentFile(path, "rb")

@functools.lru_cache
def get_sam(data_dir, filename):
    """
    Opens a SAM/BAM file using pysam.AlignmentFile.
    This function attempts to open a SAM/BAM file with the given filename in the specified directory.
    If the file is not found, it searches for a file with a suffix appended to the filename and tries to open it.

    Args:
        data_dir (str): The directory where the SAM/BAM file is located.
        filename (str): The name of the SAM/BAM file to open.

    Returns:
        pysam.AlignmentFile: An object representing the opened SAM/BAM file.

    Raises:
        FileNotFoundError: If the file is not found in the specified directory.
        ValueError: If no file with the expected suffix is found in the directory.
    """
    pysam.set_verbosity(0)
    try:
        return pysam.AlignmentFile(os.path.join(data_dir, filename), "rb")
    except:
        real_filename = [f for f in os.listdir(data_dir) if f.endswith("_"+filename)][0]
        return pysam.AlignmentFile(os.path.join(data_dir, real_filename), "rb")

def generate_reads(samfile, x, labels=None, min_base_quality=20, min_mapping_quality=0):
    """
    Generates reads at a specified position and reports whether they match the variant allele or not.

    Args:
        samfile (pysam.AlignmentFile): The SAM/BAM file to read from.
        x (dict): A dictionary containing the following keys:
            - 'chrm' (str): Chromosome name.
            - 'pos' (int): 1-based position of the variant.
            - 'ref_allele' (str): Reference allele.
            - 'alt_allele' (str): Alternate allele.
            - 'var_type' (str): Variant type ('SNP', 'DNP', 'DEL', 'INS').
        labels (dict, optional): A dictionary of labels for normal reads. Defaults to None.
        min_base_quality (int, optional): Minimum base quality for a read to be considered. Defaults to 20.
        min_mapping_quality (int, optional): Minimum mapping quality for a read to be considered. Defaults to 0.

    Yields:
        tuple: A tuple containing:
            - pileupread (pysam.PileupRead): The read from the pileup.
            - bool: True if the read matches the alternate allele, False otherwise.
            - bool: True if the read is normal according to the labels, False otherwise.

    Note:
        - Only reports reads with given alt or ref alleles and skips other alternate alleles.
        - There is processing for DEL, INS, and DNPs, but the code is not currently used. Realignment in mutect2
          affects these a lot more, so typically there are very few actual supporting reads. This would be good
          to revisit in the future.
        - It'd be nice to speed this step up, perhaps switching a portion of extract_features to directly use CLI 
          samtools. But I haven't been able to get that to work in a straightforward way. Future extension?
    """
    chrm, pos, ref, alt, var_type = x['chrm'], x['pos'], x['ref_allele'], x['alt_allele'], x['var_type']

    for pileupcolumn in samfile.pileup(chrm, pos-1, pos+1,
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
            ignore_orphans=False):

        # Pysam zero-indexes, hence pos-1. It also reports deletions at the reference base before the start
        # of the deletion, while maf reports them at the first deleted base. Hence pos-2. Need to check what VCF does
        if var_type == 'DEL' and pileupcolumn.pos != pos-2: continue
        elif var_type != 'DEL' and pileupcolumn.pos != pos-1: continue

        for i, pileupread in enumerate(pileupcolumn.pileups):

            if var_type == 'SNP':
                if pileupread.is_del or pileupread.is_refskip: continue
                allele = pileupread.alignment.query_sequence[pileupread.query_position]
                if allele != alt and allele != ref: continue
                yield pileupread, allele == alt, is_normal(pileupread, labels)

            elif var_type == 'DNP':
                if pileupread.is_del or pileupread.is_refskip: continue
                allele = pileupread.alignment.query_sequence[pileupread.query_position:pileupread.query_position+2]
                if allele != alt and allele != ref: continue
                yield pileupread, allele == alt, is_normal(pileupread, labels)

            elif var_type == 'DEL':
                if pileupread.is_refskip:
                    continue
                del_length = pileupread.indel
                if del_length == -len(ref):
                    yield pileupread, True, is_normal(pileupread, labels)
                elif del_length == 0:
                    yield pileupread, False, is_normal(pileupread, labels)
                else: continue

            elif var_type == 'INS':
                if pileupread.is_del or pileupread.is_refskip:
                    continue
                ins_length= pileupread.indel
                if ins_length == 0:
                    yield pileupread, False, is_normal(pileupread, labels)
                else:
                    seq_pos = pileupread.query_position
                    allele = pileupread.alignment.query_sequence[seq_pos+1:seq_pos+pileupread.indel+1]
                    if allele == alt:
                        yield pileupread, True, is_normal(pileupread, labels)
                    else: continue