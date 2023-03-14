import argparse
from matplotlib import pyplot
import seaborn as sns
import numpy as np
import pandas as pd

from mixture_model_EM import *

def main(args):

    # df = df_orig.sample(n=1000, random_state = 1)
    validate_arguments(args)
    input_file, output_dir, clone1, clone2 = args.input_file, args.output_dir, args.clone1, args.clone2

    df = process_df(input_file, clone1, clone2)
    df, weights, dists = initialize_EM(df)
    for i in range(2):
        prog = 0
        print("Iteration {}".format(i))
        dists, weights, R = EM(dists, weights, df)
    exit()
    print("Compute all responsibilities")
    df['responsibilities'] = df.apply(lambda x: compute_responsibilities(x, dists, weights), axis=1).tolist()

    df['assignment'] = df['responsibilities'].apply(lambda x: np.argmax(x))
    df.to_csv("{}/assignments.tsv".format(output_dir), sep='\t', index=False)

    create_plot(df, clone1, clone2, output_dir)

def add_parser_arguments(parser):
    parser.add_argument(dest='input_file', type = str, help = '<Required> Input file containing variant counts by clone')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument(dest='clone1', type = str, help = '<Required> Name of clone')
    parser.add_argument(dest='clone2',  type = str, help = '<Required> Name of clone')

def validate_arguments(args):
    pass

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
    c = [np.argmax(r) for r in df['responsibilities']]
    pyplot.scatter(df['ccf_1'],df['ccf_2'], c = c, alpha = 0.3)
    pyplot.gcf().set_size_inches(5,5)
    sns.despine()
    pyplot.xlabel('CCF Clone {}'.format(clone1))
    pyplot.ylabel('CCF Clone {}'.format(clone2))
    pyplot.savefig("{}/labeled_ccf.pdf".format(output_dir), bbox_inches='tight')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    main(args)

# TODO: Handle error rate correctly
