import argparse
import pickle
import pandas as pd
import sys
import os
import numpy as np
from multiprocessing import Pool

def main(args):
    validate_arguments(args)
    print("1. Loading model from: {}".format(args.model_dir))
    model, scaler = load_model(args.model_dir)
    print("2. Classifying data from: {}".format(args.features))

    nlines = sum(1 for _ in open(args.features, 'r'))-1

    if args.cores: ncores = args.cores
    else:
        import psutil
        ncores = psutil.cpu_count(logical=False)

    df_reader = pd.read_table(args.features, chunksize=args.chunksize)
    for i, df in enumerate(df_reader):
        print('\r\t{}/{} variants completed'.format(i*args.chunksize, nlines), end = '')        
        first = (i == 0)
        process_chunk(df, model, scaler, args.output_dir, ncores, first)
    print('\r\t{}/{} variants completed'.format(nlines, nlines), end = '')

def process_chunk(df, model, scaler, output_dir, ncores, first):
    input_data = scale_data(df, scaler)
    labels, probs = predict(df, input_data, model, ncores)
    write_output_chunk(df, labels, probs, output_dir, first)

def scale_data(df, scaler):
    features = [c for c in df.columns if c.startswith('f_')]
    data = df[features].values

    scaled_data = scaler.transform(data)
    return scaled_data

def write_output_chunk(df, labels, probs, output_dir, first):
    df['temp'] = labels
    df['result'] = df['temp'].apply(lambda x: 'PASS' if x else "ARTIFACT")
    df['prob'] = probs

    out_df = df[['chrm', 'pos', 'result', 'prob']]
    out_file = os.path.join(output_dir, 'result.tsv')
    if first:
        out_df.to_csv(out_file, sep='\t', index=False)
    else:
        out_df.to_csv(out_file, mode = 'a', header = False, sep='\t', index=False)

def add_parser_arguments(parser):
    parser.add_argument(dest='features', type = str, help = '<Required> Input file containing variant features')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument(dest='model_dir', type = str, help = '<Required> Directory containing model.pkl and scaler.plk')
    
    DEFAULT_CHUNKSIZE = 5000
    parser.add_argument('--chunksize', type = int, default = DEFAULT_CHUNKSIZE, required = False,
                        help = F'<Optional> Number of rows per worker (default {DEFAULT_CHUNKSIZE})')

def validate_arguments(args):
    pass

def load_model(model_dir):
    with open(os.path.join(model_dir, 'model.pkl'), 'rb') as f:
        model = pickle.load(f)
        model.set_params(n_jobs = 10)

    with open(os.path.join(model_dir, 'scaler.pkl'), 'rb') as f:
        scaler = pickle.load(f)
    return model, scaler

def predict(df, input_data, model, ncores):
    scaled_features = ['{}_s'.format(c) for c in df.columns if c.startswith('f_')]
    df[scaled_features] = input_data

    nperchunk = len(df)/ncores
    df['chunk'] = (np.arange(len(df))/nperchunk).astype(int)

    tasks = []
    for _, df_ in df.groupby('chunk'):
        tasks.append((model, df_[scaled_features].values.copy()))

    with Pool(ncores) as p:
        results = p.map(predict_task, tasks)
    
    probs = np.concatenate(results)
    labels = [p < 0.5 for p in probs]
    return labels, probs

def predict_task(params):
    model, data = params
    return model.predict_proba(data)[:, 0]

def write_output(df, labels, probs, output_dir):
    df['temp'] = labels
    df['result'] = df['temp'].apply(lambda x: 'PASS' if x == 1 else "ARTIFACT")
    df['prob'] = probs
    df[['chrm', 'pos', 'result', 'prob']].to_csv(os.path.join(output_dir, 'result.tsv'), sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    main(args)
