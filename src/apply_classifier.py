import argparse
import pickle
import pandas as pd
import sys
import os
import numpy as np
from multiprocessing import Pool

def classify(model_dir, features, output_dir, chunksize, ncores):

    validate_arguments(model_dir, features, output_dir)
    print("1. Loading model from: {}".format(model_dir))
    model, scaler = load_model(model_dir)
    print("2. Classifying data from: {}".format(features))
    nlines = sum(1 for _ in open(features, 'r'))-1

    if not ncores: 
        import psutil
        ncores = psutil.cpu_count(logical=False)

    df_reader = pd.read_table(features, chunksize=chunksize)
    for i, df in enumerate(df_reader):
        print('\r\t{}/{} variants completed'.format(i*chunksize, nlines), end = '')
        first = (i == 0)
        df['f_p_normal'] = df['f_p_normal'].fillna(0)
        process_chunk(df, model, scaler, output_dir, ncores, first)
    print('\r\t{}/{} variants completed'.format(nlines, nlines), end = '')

def process_chunk(df, model, scaler, output_dir, ncores, first):
    orig_df = df

    ###
    #   TEMPORARY: For now, only classify SNPs. Other models are not implemented
    ###
    input_data = scale_data(df[df['var_type'].isin(['SNP', 'DNP'])].dropna(), scaler)
    probs = predict(df[df['var_type'].isin(['SNP', 'DNP'])].dropna(), input_data, model, ncores)
    write_output_chunk(df, probs, output_dir, first)

def scale_data(df, scaler):
    features = [c for c in df.columns if c.startswith('f_')]
    data = df[features].values

    scaled_data = scaler.transform(data)
    return scaled_data

def write_output_chunk(df, probs, output_dir, first):
    df['prob_artifact'] = probs
    df['result'] = df['prob_artifact'].apply(lambda x: 'PASS' if x < 0.5 else "ARTIFACT" if x >= 0.5 else "SKIP")
    df['result'] = df['result'].fillna('SKIP')

    out_df = df[['chrm', 'pos',  'ref_allele', 'alt_allele', 'result', 'prob_artifact']]
    out_file = os.path.join(output_dir, 'result.tsv')
    if first:
        out_df.to_csv(out_file, sep='\t', index=False)
    else:
        out_df.to_csv(out_file, mode = 'a', header = False, sep='\t', index=False)


def validate_arguments(model_dir, features, output_dir):
    assert os.path.isfile(features), f"Input file {features} does not exist."
    assert os.access(features, os.R_OK), f"Input file exists, but cannot be read due to permissions: {features}"
    assert os.path.isdir(model_dir), f"Model dir {model_dir} does not exist."
    assert os.path.isfile(os.path.join(model_dir, 'model.pkl')), f"Model file {os.path.join(model_dir, 'model.pkl')} does not exist."
    assert os.path.isfile(os.path.join(model_dir, 'scaler.pkl')), f"Scaler file {os.path.join(model_dir, 'scaler.pkl')} does not exist."
    assert os.path.isdir(output_dir), f"Output directory {output_dir} does not exist."
    assert os.access(output_dir, os.W_OK), f"Output directory exists, but cannot be written to due to permissions: {output_dir}"

def load_model(model_dir):
    with open(os.path.join(model_dir, 'model.pkl'), 'rb') as f:
        model = pickle.load(f)
        #model.set_params(n_jobs = 10)

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
    return pd.Series(data = probs, index = df.index)

def predict_task(params):
    model, data = params
    return model.predict_proba(data)[:, 0]

def write_output(df, labels, probs, output_dir):
    df['temp'] = labels
    df['result'] = df['temp'].apply(lambda x: 'PASS' if x == 1 else "ARTIFACT")
    df['prob'] = probs
    df[['chrm', 'pos', 'result', 'prob']].to_csv(os.path.join(output_dir, 'result.tsv'), sep='\t', index=False)