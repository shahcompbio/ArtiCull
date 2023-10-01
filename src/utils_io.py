import pandas as pd
import io
import os

import gzip

def get_variants(filename, file_type = "maf"):
    if filename.endswith('.maf'):
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
        
        

def update_progress(tot):
    # Could be shared with extract features?
    prog = 0
    while True:
        yield prog
        prog += 1
        if prog % 20 == 0: print("\t Progress: {}/{}".format(prog, tot))
