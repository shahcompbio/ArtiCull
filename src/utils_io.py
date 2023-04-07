import pandas as pd

def get_variants(filename):
    # Could be shared with extract features?
    df = pd.read_table(filename, skiprows=1)
    df['chrm'] = df['Chromosome'].astype('str')
    df['pos'] = df['Start_Position']
    print(df.shape)
    df = df[(df['Variant_Type'].apply(lambda x: 'SNP' in x))] #&(df['FILTER'] == 'PASS')]
    df['ref_allele'] = df['Reference_Allele']
    df['alt_allele'] = df['Tumor_Seq_Allele2']

    print(df.shape)
    return df[['chrm', 'pos', 'ref_allele', 'alt_allele']]

def update_progress(tot):
    # Could be shared with extract features?
    prog = 0
    while True:
        yield prog
        prog += 1
        if prog % 20 == 0: print("\t Progress: {}/{}".format(prog, tot))
