import os
import argparse
from collections import Counter
import pandas as pd
pd.options.mode.chained_assignment = None #Suppress SettingWithACopy Warning

from maf_input import get_variants
from utils import update_progress
from process_bams import match_variants_to_filenames, get_sam, generate_reads

def parse_arguments():
    def get_arguments():
        parser = argparse.ArgumentParser()
        parser.add_argument(dest='maf', type = str, help = '<Required> maf file containing candidate variants')
        parser.add_argument(dest='clonemap', type = str, help = '<Required> Map of cells to clones')
        parser.add_argument(dest='clonecns', type = str, help = '<Required> Map of cells to clones')
        parser.add_argument(dest='output', type = str, help = '<Required> Full path and name of output file')
        parser.add_argument(dest='bam_dirs', nargs="+", type = str, help = '<Required> list of bam directories')
        args = parser.parse_args()
        return args

    def validate_arguments(args):
        # Checks if input files exist and if output files are in directories that exist and can be written to
        assert os.path.isfile(args.maf)
        assert os.path.isfile(args.clonemap)
        assert os.path.isfile(args.clonecns)
        for dir in args.bam_dirs:
            assert os.path.isdir(dir)
        assert os.access(os.path.dirname(args.output), os.W_OK)

    args = get_arguments()
    validate_arguments(args)
    return args.maf, args.bam_dirs, args.clonemap, args.clonecns, args.output


def get_clone_var_counts(df, data_dirs, clonemap, clone_ids):

    def get_clone_list(reads):
        cells = ['_'.join(read.split('_')[1:-2]) for read in reads]
        clones= []
        for cell in cells:
            try:
                clones.append(clonemap.loc[cell]['clone_id'])
            except:
                clones.append('0')
        return clones

    def process_variant(x):
        next(prog)

        var_sams = [get_sam(data_dir, x['filename']) for data_dir in data_dirs]

        reads = []
        for sam in var_sams: reads+=[r for r in generate_reads(sam, x)]

        read_groups = [r[0].alignment.get_tag('RG') for r in reads]
        clones = get_clone_list(read_groups)

        var_count = Counter()
        all_count = Counter()

        for c,(r, is_alt) in zip(clones, reads):
            if is_alt: var_count[c]+=1
            all_count[c] += 1

        def get_var_data(c):
            return [var_count[c], all_count[c], cns[c]]

        var_data = []
        for c in clone_ids:
            var_data += [var_count[c], all_count[c]]

        return var_data

    df = match_variants_to_filenames(df, data_dirs)
    prog = update_progress(len(df))

    columns = []
    for c in clone_ids:
            columns+=['var_{}'.format(c), 'tot_{}'.format(c)]

    df[columns] = df.apply(process_variant, axis=1, result_type="expand")

    del df['filename']
    return df

def get_clone_cns(df, clonemap, clone_cn_file, clone_ids):

    clone_cns = pd.read_table(clone_cn_file)

    clone_cn_dfs = {}
    for chrm in clone_cns['chr'].unique():
        clone_cns1 = clone_cns[clone_cns['chr'] == chrm]
        clone_cns1['CN'] = clone_cns['Maj'] + clone_cns['Min']
        clone_cns1 = clone_cns1.pivot(index=['start', 'end'], columns = 'clone', values = 'CN')
        clone_cns1.index = pd.IntervalIndex.from_tuples(clone_cns1.index, closed='both')
        clone_cn_dfs[chrm] = clone_cns1

    #for i in clone_cn_dfs['5'].index: print(i)

    columns = ['cn_{}'.format(c) for c in clone_ids]

    def get_cn(x):
        try:
            return clone_cn_dfs[x['chrm']].loc[int(x['pos'])]
        except KeyError:
            return [float('nan')]*len(columns)

    df[columns] = df.apply(get_cn, axis=1, result_type="expand")
    return df

def compute_ccfs(df, clone_ids):
    for clone in clone_ids:
        ccf = 'ccf_{}'.format(clone)
        var = 'var_{}'.format(clone)
        cn = 'cn_{}'.format(clone)
        tot = 'tot_{}'.format(clone)

        df[ccf] = df[var] *  df[cn] / df[tot]

    return df

def read_clonemap(clonemap_file):
    clonemap = pd.read_table(clonemap_file)
    clone_ids = sorted(clonemap['clone_id'].unique())
    return clonemap, clone_ids


def main():
    maf, bam_dirs, clonemap_file, clone_cns, output = parse_arguments()
    print("1. Reading variants from: {}\n".format(maf))
    df = get_variants(maf)
    clonemap, clone_ids = read_clonemap(clonemap_file)

    print("2. Extracting clone variant counts from: {}".format(bam_dirs))
    df = get_clone_var_counts(df, bam_dirs, clonemap, clone_ids)

    print("3. Getting clone copy numbers and computing CCFs from : {}\n".format(clone_cns))
    df = get_clone_cns(df, clonemap, clone_cns, clone_ids)
    df = compute_ccfs(df, clone_ids)


    print("3. Outputting result to: {}\n".format(output))
    df.to_csv(output, sep = '\t', index=False)

if __name__ == "__main__":
    main()
