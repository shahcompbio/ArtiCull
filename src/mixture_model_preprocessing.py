import os
from os import path
import argparse
from collections import Counter
import pandas as pd
pd.options.mode.chained_assignment = None #Suppress SettingWithACopy Warning
import subprocess
import tempfile
from matplotlib import pyplot
import seaborn as sns

from utils_io import get_variants, update_progress
from utils_bams import match_variants_to_filenames, get_sam, generate_reads


def main(args):
    #maf, bam_dirs, signals_dir, output = parse_arguments()
    validate_arguments(args)
    maf, bam_dirs, signals_dir, output_dir, fullbam = args.maf, args.bam_dirs, args.signals_dir, args.output_dir, args.fullbam

    print("1. Reading variants from: {}".format(maf))
    df = get_variants(maf)

    print("2. Process signals results from: {}".format(signals_dir))
    #clonemap, clone_ids, clone_cns = process_signals(signals_dir)
    df, clonemap, clone_ids = process_signals(df, signals_dir)

    print(clone_ids)

    print("3. Extracting clone variant counts from: {}".format(bam_dirs))
    df = get_clone_var_counts(df, bam_dirs, clonemap, clone_ids, fullbam)

    print("4. Computing CCFs")
    df = compute_ccfs(df, clone_ids)

    result_filename = os.path.join(output_dir, "var_counts.tsv")
    print("5. Outputting result to: {}".format(result_filename))
    df.to_csv(result_filename, sep = '\t', index=False)

    plot_filename = os.path.join(output_dir, "clone_ccfs.png")
    print("6. Creating clone CCF plot at: {}".format(plot_filename))
    plot_ccfs(df, plot_filename, clone_ids)

def add_parser_arguments(parser):
    parser.add_argument(dest='maf', type = str, help = '<Required> maf file containing candidate variants')
    parser.add_argument(dest='signals_dir', type = str, help = '<Required> directory of signals output')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Full path and name of output file')
    parser.add_argument(dest='bam_dirs', nargs="+", type = str, help = '<Required> list of bam directories')
    parser.add_argument('--fullbam', action="store_true", help ='The list of bams is provided and not in region format')

def validate_arguments(args):
    # Checks if input files exist and if output files are in directories that exist and can be written to
    for arg in vars(args):
        print(arg, getattr(args, arg))

    assert path.isfile(args.maf)
    assert path.isdir(args.signals_dir)
    #for dir in args.bam_dirs:
    #    assert path.isdir(dir)
    assert os.access(args.output_dir, os.W_OK)

def process_signals(df, signals_dir):

    signals_result = path.join(signals_dir, 'signals.Rdata')
    signals_cns = path.join(signals_dir, 'hscn.csv.gz')
    temp = tempfile.NamedTemporaryFile(delete=False, mode='w')
    temp_name = temp.name
    #temp_name = "./cn.txt"
    try:
        #temp.close()
        directory = path.dirname(path.realpath(path.expanduser(__file__)))
        command = "Rscript {}/extract_cell2clone.R {} {}".format(directory, signals_result, temp_name)
        #print(command)
        subprocess.check_call(command.split(' '))

        clonemap = pd.read_table(temp_name)
        clone_ids = sorted(clonemap['clone_id'].unique())
    finally:
        os.remove(temp_name)

    cn_df = pd.read_table(signals_cns, sep=',')

    # Cell to clone map table
    ids_map = clonemap['clone_id']
    cn_df['clone']= cn_df['cell_id'].map(ids_map)

    clone_cns = cn_df.groupby(['chr', 'start', 'end', 'clone'])[['Maj', 'Min']].median().reset_index()

    clone_cn_dfs = {}
    for chrm in clone_cns['chr'].unique():
        clone_cns1 = clone_cns[clone_cns['chr'] == chrm]
        clone_cns1['CN'] = clone_cns['Maj'] + clone_cns['Min']
        clone_cns1 = clone_cns1.pivot(index=['start', 'end'], columns = 'clone', values = 'CN')
        clone_cns1 = clone_cns1[sorted(clone_cns1.columns)]
        clone_cns1.index = pd.IntervalIndex.from_tuples(clone_cns1.index, closed='both')
        clone_cn_dfs[chrm] = clone_cns1

    columns = ['cn_{}'.format(c) for c in clone_ids]

    def get_cn(x):
        try:
            result = clone_cn_dfs[x['chrm']].loc[int(x['pos'])].values
            return result

        except KeyError:
            return [float('nan')]*len(columns)

    result = df.parallel_apply(get_cn, axis=1, result_type="expand")
    df[columns] = result

    return df, clonemap, clone_ids

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
        #next(prog)

        var_sams = [get_sam(data_dir, x['filename']) for data_dir in data_dirs]

        reads = []
        for sam in var_sams: reads+=[r for r in generate_reads(sam, x)]

        read_groups = [r[0].alignment.get_tag('RG') for r in reads]
        clones = get_clone_list(read_groups)

        var_count = Counter()
        all_count = Counter()

        for c, r in zip(clones, reads):
            if r[1]: var_count[c]+=1
            all_count[c] += 1

        def get_var_data(c):
            return [var_count[c], all_count[c], cns[c]]

        var_data = []
        for c in clone_ids:
            var_data += [var_count[c], all_count[c]]

        return var_data

    if fullbam: df['filename'] = data_dirs[0]
    else:
        df = match_variants_to_filenames(df, data_dirs)

    prog = update_progress(len(df))

    columns = []
    for c in clone_ids:
            columns+=['var_{}'.format(c), 'tot_{}'.format(c)]

    df[columns] = df.parallel_apply(process_variant, axis=1, result_type="expand")

    del df['filename']
    return df

def compute_ccfs(df, clone_ids):
    for clone in clone_ids:
        ccf = 'ccf_{}'.format(clone)
        var = 'var_{}'.format(clone)
        cn = 'cn_{}'.format(clone)
        tot = 'tot_{}'.format(clone)

        df[ccf] = df[var] *  df[cn] / df[tot]

    return df

def plot_ccfs(df, output, clone_ids):
    sns.set_context('paper', font_scale=1.5)
    df_plot = pd.DataFrame()
    for c in clone_ids:
        df_plot[c] = df['ccf_{}'.format(c)]

    g = sns.pairplot(df_plot, kind='scatter', plot_kws={'alpha':0.2, 's':8}, corner=True)
    g.set(xlim=(-0.1,2.1), ylim = (-0.1,2.1))
    pyplot.savefig(output, bbox_inches = 'tight')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    main(args)
