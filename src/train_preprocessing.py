import os
from os import path
import argparse
from collections import Counter
import pandas as pd # type: ignore
pd.options.mode.chained_assignment = None #Suppress SettingWithACopy Warning
import subprocess
from matplotlib import pyplot # type: ignore
import seaborn as sns # type: ignore

from utils_io import get_variants
from utils_bams import match_variants_to_filenames, get_sam, generate_reads


def main(args):
    validate_arguments(args)

    maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn = \
        args.maf, args.bam_dirs, args.signals_dir, args.output_dir, args.fullbam, args.cell_clones, args.hscn, args.use_cached_cn

    print("1. Reading variants from: {}".format(maf))
    df = get_variants(maf)

    print("2. Process copy number information")
    if signals_dir: 
        cellclone_file, hscn_file = process_signals(df, signals_dir, output_dir, use_cached_cn)
        df, clonemap, clone_ids = parse_copynumber(df, cellclone_file, hscn_file) 
    elif cellclone_file and hscn_file:
        df, clonemap, clone_ids = parse_copynumber(df, cellclone_file, hscn_file) 
    else:
        raise RuntimeError('Missing files: Must provide either --signals_dir argument, or --cell_clones AND --hscn')

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
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Full path and name of output file')
    parser.add_argument(dest='bam_dirs', nargs="+", type = str, help = '<Required> list of bam directories')
    parser.add_argument('--signals_dir', type = str, help = 'Directory of signals results')
    parser.add_argument('--cell_clones', type = str, help = 'Cell to clone mapping')
    parser.add_argument('--hscn', type = str, help='Signals HSCN file')
    parser.add_argument('--fullbam', action="store_true", help ='The list of bams is provided and not in region format')
    parser.add_argument('--use_cached_cn', action="store_true", help = 'Use already processed cell-to-clone map if it exists')

def validate_arguments(args):
    # Checks if input files exist and if output files are in directories that exist and can be written to
    for arg in vars(args):
        print(arg, getattr(args, arg))

    assert path.isfile(args.maf)
    #assert path.isdir(args.signals_dir)
    #for dir in args.bam_dirs:
    #    assert path.isdir(dir)
    assert os.access(args.output_dir, os.W_OK)

def process_signals(df, signals_dir, output_dir, use_cached):

    signals_result = path.join(signals_dir, 'signals.Rdata')
    signals_cns = path.join(signals_dir, 'hscn.csv.gz')
    outfile = os.path.join(output_dir, "signals_clones.tsv")

    if use_cached and os.path.isfile(outfile):
        print("Using existing cell-to-clone map: {}".format(outfile))
        return outfile, signals_cns
    else:
        directory = path.dirname(path.realpath(path.expanduser(__file__)))
        command = "Rscript {}/extract_cell2clone.R {} {}".format(directory, signals_result, outfile)
        subprocess.check_call(command.split(' '))
        return outfile, signals_cns


def parse_copynumber(df, cellclone_file, hscn_file):

    clonemap = pd.read_table(cellclone_file, index_col=0)
    clone_ids = sorted(clonemap['clone_id'].unique())
    ids_map = clonemap['clone_id']

    try:
        cn_df = pd.read_table(hscn_file, usecols=['cell_id', 'chr', 'start', 'end', 'Maj', 'Min'], sep=',')
    except ValueError:
        cn_df = pd.read_table(hscn_file, usecols=['cell_id', 'chr', 'start', 'end', 'A', 'B'], sep=',')
        cn_df['Maj'] = cn_df['A']
        cn_df['Min'] = cn_df['B']

    # Cell to clone map table    
    cn_df['clone']= cn_df['cell_id'].map(ids_map)
    clone_cns = cn_df.groupby(['chr', 'start', 'end', 'clone'])[['Maj', 'Min']].median().reset_index()
    clone_cns['CN'] = clone_cns['Maj'] + clone_cns['Min']

    del cn_df # large file
    clone_cn_dfs = {}

    for chrm in clone_cns['chr'].unique():
        clone_cns1 = clone_cns[clone_cns['chr'] == chrm]
        #clone_cns1['CN'] = clone_cns['Maj'] + clone_cns['Min']
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

def get_clone_var_counts(df, data_dirs, clonemap, clone_ids, fullbam=False):

    def get_clone_list(reads):
        #cells = ['_'.join(read.split('_')[1:-2]) for read in reads]
        cells = [read for read in reads]
        clones= []
        for cell in cells:
            try:
                clones.append(clonemap.loc[cell]['clone_id'])
            except:
                #print(cell)
                #raise
                clones.append('0')
        return clones

    def process_variant(x):
        var_sams = [get_sam(data_dir, x['filename']) for data_dir in data_dirs]

        reads = []
        for sam in var_sams: reads+=[r for r in generate_reads(sam, x)]

        read_groups = [r[0].alignment.get_tag('CB:Z:') for r in reads]
        #read_groups = [r[0].alignment.get_tag('RG') for r in reads]
        clones = get_clone_list(read_groups)

        var_count = Counter()
        all_count = Counter()

        for c, r in zip(clones, reads):
            if r[1]: var_count[c]+=1
            all_count[c] += 1

        var_data = []
        for c in clone_ids:
            var_data += [var_count[c], all_count[c]]

        return var_data

    if fullbam: df['filename'] = data_dirs[0]
    else:
        df = match_variants_to_filenames(df, data_dirs)

    columns = []
    for c in clone_ids: columns+=['var_{}'.format(c), 'tot_{}'.format(c)]

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
    g.set(xlim=(-0.1,1.1), ylim = (-0.1,1.1))
    pyplot.savefig(output, bbox_inches = 'tight')



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    main(args)
