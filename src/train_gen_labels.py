import argparse
from matplotlib import pyplot
import seaborn as sns
import numpy as np
import pandas as pd
import os
import scipy.stats

def main(input_file, output_dir, clone1, clone2, alpha):

    validate_arguments(input_file, output_dir)
    df = process_df(input_file, clone1, clone2)
    df['assignment'] = df.parallel_apply(lambda x: get_assignment(x, alpha), axis=1)
    df.to_csv("{}/assignments.tsv".format(output_dir), sep='\t', index=False)
    print(df.groupby('assignment').count())
    create_plot(df, clone1, clone2, output_dir)

def shared_subclonal(mut, alpha):
    v1, t1 = mut['var_1'], mut['tot_1']
    v2, t2 = mut['var_2'], mut['tot_2']
    c1, c2 = mut['cn_1'], mut['cn_2']

    x = 0.8
    x_adj = x/c1 - (2*x*eps)/c1 + eps
    clonal1 = scipy.stats.binomtest(v1, n=t1, p=x_adj, alternative='less').pvalue
    x = 0.8
    x_adj = x/c2 - (2*x*eps)/c2 + eps
    clonal2 = scipy.stats.binomtest(v2, n=t2, p=x_adj, alternative='less').pvalue

    private1 = scipy.stats.binomtest(v1, n=t1, p=eps, alternative='greater').pvalue

    private2 = scipy.stats.binomtest(v2, n=t2, p=eps, alternative='greater').pvalue

    return all([p < alpha for p in [clonal1, clonal2, private1, private2]]) #and mut['ccf_1'] < 0.75 and mut['ccf_2'] < 0.75

def shared_clonal(mut, alpha):
    v1, t1 = mut['var_1'], mut['tot_1']
    v2, t2 = mut['var_2'], mut['tot_2']
    c1, c2 = mut['cn_1'], mut['cn_2']
    x = 0.6
    x_adj = x/c1 - (2*x*eps)/c1 + eps
    clonal1 = scipy.stats.binomtest(v1, n=t1, p=x_adj, alternative='greater').pvalue
    x = 0.6
    x_adj = x/c2 - (2*x*eps)/c2 + eps
    clonal2 = scipy.stats.binomtest(v2, n=t2, p=x_adj, alternative='greater').pvalue

    return (clonal1 < alpha and clonal2 < alpha)

def get_assignment(mut, alpha):
    c1, c2 = mut['cn_1'], mut['cn_2']
    if c1 == 0 or c2 == 0: return -1
    if shared_subclonal(mut, alpha): return 0
    elif shared_clonal(mut, alpha): return 1
    else: return -1



def process_df(input_file, clone1, clone2):
    input = pd.read_table(input_file)
    df = pd.DataFrame()
    df['chrm'] = input['chrm'].astype('str')
    df['pos'] = input['pos']

    for val in ['var', 'tot', 'cn', 'ccf']:
        df['{}_1'.format(val)] = input['{}_{}'.format(val, clone1)]
        df['{}_2'.format(val)] = input['{}_{}'.format(val, clone2)]

    df = df.dropna()
    return df

def create_plot(df, clone1, clone2, output_dir):
    pyplot.gcf().set_size_inches(12,12)
    df = df.sort_values(by = 'assignment')
    c = df['assignment']
    col = ['r', 'g', 'b']
    for d,l in {-1:'?', 0:'Artifact', 1:'Real'}.items():
        plot_df = df[df['assignment'] == d]
        pyplot.scatter(plot_df['ccf_1'],plot_df['ccf_2'], alpha = 0.2, s = 2, label=l)

    pyplot.gcf().set_size_inches(5,5)
    sns.despine()
    pyplot.legend()
    pyplot.xlabel('CCF Clone {}'.format(clone1))
    pyplot.ylabel('CCF Clone {}'.format(clone2))
    pyplot.savefig("{}/labeled_ccf.pdf".format(output_dir), bbox_inches='tight')


def validate_arguments(input_file, output_dir):
    assert os.path.isfile(input_file), f"Input file {input_file} does not exist."
    assert os.access(input_file, os.R_OK), f"Input file exists, but cannot be read due to permissions: {input_file}"
    assert os.path.isdir(output_dir), f"Output directory {output_dir} does not exist."
    assert os.access(output_dir, os.W_OK), f"Output directory exists, but cannot be written to due to permissions: {output_dir}"

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    main(args)