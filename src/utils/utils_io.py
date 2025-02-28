"""
utils_io.py

This module provides utility functions for reading and processing variant data from MAF and VCF files.
It includes functionality for extracting variant information and reading VCF files.

Functions:
    get_variants(filename, file_type="maf"):
        Extracts variant information from a given file in MAF or VCF format.

    read_vcf(path):
        Reads a VCF (Variant Call Format) file, optionally compressed with gzip, and returns a pandas DataFrame.
"""

import pandas as pd # type: ignore
import io
import gzip

def get_variants(filename, file_type = "maf"):
    """
    Extracts variant information from a given file in MAF or VCF format.
    Parameters:
    filename (str): The path to the input file. The file can be in MAF (.maf or .maf.gz) or VCF (.vcf or .vcf.gz) format.
    file_type (str, optional): The type of the input file. Default is "maf".
    Returns:
    pandas.DataFrame: A DataFrame containing the following columns:
        - 'chrm': Chromosome
        - 'pos': Position
        - 'ref_allele': Reference allele
        - 'alt_allele': Alternate allele
        - 'var_type': Variant type (only 'SNP' variants are included)
    Raises:
    Exception: If the input file type is not MAF or VCF.
    """

    if filename.endswith('.maf') or filename.endswith('.maf.gz'):
        # Could be shared with extract features?
        df = pd.read_table(filename, skiprows=1)
        df['chrm'] = df['Chromosome'].astype('str')
        df['pos'] = df['Start_Position']
        df['ref_allele'] = df['Reference_Allele']
        df['alt_allele'] = df['Tumor_Seq_Allele2']
        df['var_type'] = df['Variant_Type'].apply(lambda x: 'SNP' if 'SNP' in x else x)
        df = df[df['var_type'] == 'SNP']
        return df[['chrm', 'pos', 'ref_allele', 'alt_allele', 'var_type']]
    elif filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
        df = read_vcf(filename)
        df['chrm'] = df['CHROM']
        df['pos'] = df['POS']
        df['ref_allele'] = df['REF']
        df['alt_allele'] = df['ALT']
        df = df[df['FILTER'] == 'PASS']
        # Filter for SNPs
        df=df[df.ref_allele.apply(lambda x: len(str(x))==1 and '-' not in x)]
        df=df[df.alt_allele.apply(lambda x: len(str(x))==1 and '-' not in x)]
        df['var_type'] = 'SNP'
        return df[['chrm', 'pos', 'ref_allele', 'alt_allele', 'var_type']]
    else:
        raise Exception("Invalid input filetype specified")

def read_vcf(path):
    """
    Reads a VCF (Variant Call Format) file, optionally compressed with gzip, and returns a pandas DataFrame.
    Parameters:
    path (str): The file path to the VCF file. The file can be plain text or gzip compressed.
    Returns:
    pandas.DataFrame: A DataFrame containing the VCF data with columns ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'].
    """

    if path.endswith('.gz'):
        with gzip.open(path,'rt') as f:
            lines = [l for l in f if not l.startswith('##')]
    else: 
        with open(path, 'r') as f:
            lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})