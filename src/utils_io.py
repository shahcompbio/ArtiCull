import pandas as pd

def get_variants(filename):
    # Could be shared with extract features?
    df = pd.read_table(filename, skiprows=1)
    df['chrm'] = df['Chromosome'].astype('str')
    df['pos'] = df['Start_Position']
    df['ref_allele'] = df['Reference_Allele']
    df['alt_allele'] = df['Tumor_Seq_Allele2']
    df['var_type'] = df['Variant_Type'].apply(lambda x: 'SNP' if 'SNP' in x else x)
    return df[['chrm', 'pos', 'ref_allele', 'alt_allele', 'var_type']]

def update_progress(tot):
    # Could be shared with extract features?
    prog = 0
    while True:
        yield prog
        prog += 1
        if prog % 20 == 0: print("\t Progress: {}/{}".format(prog, tot))
