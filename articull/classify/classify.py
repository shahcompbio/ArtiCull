"""
This module runs a pre-trained ArtiCull model to classify variants as artifacts or true mutations.

Functions:
    classify(model_dir, features, output_dir, chunksize, ncores):
        Classifies data using a pre-trained model and saves the results to the specified output directory.
"""

import pickle
import pandas as pd # type: ignore
import os
import numpy as np
from multiprocessing import Pool

def classify_variants(model_dir, features, output_prefix, chunksize, ncores):
    """
    Classifies data using a pre-trained model.

    Args:
        model_dir (str): Directory where the model and scaler are stored.
        features (str): Path to the file containing features to classify.
        output_prefix (str): Prefix to be used for naming the output files.
        chunksize (int): Number of lines to read at a time from the features file.
        ncores (int): Number of CPU cores to use for processing. If None, the number of physical cores will be used.

    Returns:
        None

    This function performs the following steps:
        1. Validates the input arguments.  
        2. Loads the model and scaler from the specified directory.  
        3. Reads the features file in chunks and classifies the data.  
        4. Saves the classification results to the specified output directory.  
    """
    _validate_arguments(model_dir, features, output_prefix)
    print("1. Loading model from: {}".format(model_dir))
    model, scaler = _load_model(model_dir)
    print("2. Classifying data from: {}".format(features))
    nlines = sum(1 for _ in open(features, 'r')) - 1

    if not ncores: 
        import psutil
        ncores = psutil.cpu_count(logical=False)

    df_reader = pd.read_table(features, chunksize=chunksize)
    for i, df in enumerate(df_reader):
        print('\r\t{}/{} variants completed'.format(i * chunksize, nlines), end='')
        first = (i == 0)
        df['f_p_normal'] = df['f_p_normal'].fillna(0)
        _process_chunk(df, model, scaler, output_prefix, ncores, first)
    print('\r\t{}/{} variants completed'.format(nlines, nlines), end='')

def _process_chunk(df, model, scaler, output_prefix, ncores, first):
    """
    Processes a chunk of data by classifying SNPs using the provided model and scaler.

    Args:
        df (pandas.DataFrame): The input dataframe containing the data to be processed.
        model: The classification model to be used for predicting probabilities.
        scaler: The scaler object used to scale the input data.
        output_prefix (str): The prefix to be used for naming the output files.
        ncores (int): The number of cores to be used for parallel processing.
        first (bool): A flag indicating if this is the first chunk being processed.

    Returns:
        None
    """
    orig_df = df

    ###
    #   TEMPORARY: For now, only classify SNPs. Other models are not implemented
    ###
    input_data = _scale_data(df[df['var_type'].isin(['SNP', 'DNP'])].dropna(), scaler)
    probs = _predict(df[df['var_type'].isin(['SNP', 'DNP'])].dropna(), input_data, model, ncores)
    _write_output_chunk(df, probs, output_prefix, first)

def _scale_data(df, scaler):
    """
    Scales the feature columns of a DataFrame using the provided scaler.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be scaled.
        scaler (object): The scaler object (e.g., from sklearn.preprocessing) used to scale the data.

    Returns:
        numpy.ndarray: The scaled data as a NumPy array.
    """
    features = [c for c in df.columns if c.startswith('f_')]
    data = df[features].values

    scaled_data = scaler.transform(data)
    return scaled_data

def _write_output_chunk(df, probs, output_prefix, first):
    """
    Writes a chunk of the output DataFrame to a TSV file.

    Args:
        df (pandas.DataFrame): The input DataFrame containing the data to be classified.
        probs (array-like): The probabilities of each row being an artifact.
        output_prefix (str): The prefix to be used for naming the output files.
        first (bool): A flag indicating whether this is the first chunk being written. 
                    If True, the file will be created with a header. 
                    If False, the data will be appended to the existing file without a header.

    Returns:
        None
    """
    df['prob_artifact'] = probs
    df['result'] = df['prob_artifact'].apply(lambda x: 'PASS' if x < 0.5 else "ARTIFACT" if x >= 0.5 else "SKIP")
    df['result'] = df['result'].fillna('SKIP')

    out_df = df[['chrm', 'pos', 'ref_allele', 'alt_allele', 'result', 'prob_artifact']]
    out_file = output_prefix +"_result.tsv"
    if first:
        out_df.to_csv(out_file, sep='\t', index=False)
    else:
        out_df.to_csv(out_file, mode='a', header=False, sep='\t', index=False)

def _validate_arguments(model_dir, features, output_prefix):
    """
    Validates the input arguments for the classifier application.

    Args:
        model_dir (str): Path to the directory containing the model and scaler files.
        features (str): Path to the input features file.
        output_prefix (str): Prefix to be used for naming the output files.

    Raises:
        AssertionError: If any of the following conditions are not met:
            - The features file exists and is readable.
            - The model directory exists.
            - The model file ('model.pkl') exists in the model directory.
            - The scaler file ('scaler.pkl') exists in the model directory.
            - The output directory exists and is writable.
    """
    assert os.path.isfile(features), f"Input file {features} does not exist."
    assert os.access(features, os.R_OK), f"Input file exists, but cannot be read due to permissions: {features}"
    assert os.path.isdir(model_dir), f"Model dir {model_dir} does not exist."
    assert os.path.isfile(os.path.join(model_dir, 'model.pkl')), f"Model file {os.path.join(model_dir, 'model.pkl')} does not exist."
    assert os.path.isfile(os.path.join(model_dir, 'scaler.pkl')), f"Scaler file {os.path.join(model_dir, 'scaler.pkl')} does not exist."
    output_dir = os.path.dirname(output_prefix)
    assert os.path.isdir(output_dir), f"Output directory {output_dir} does not exist."
    assert os.access(output_dir, os.W_OK), (
        f"Output directory exists, but cannot be written to due to permissions: {output_dir}"
    )

def _load_model(model_dir):
    """
    Loads a machine learning model and its corresponding scaler from the specified directory.

    Args:
        model_dir (str): The directory where the model and scaler files are stored.

    Returns:
        tuple: A tuple containing the loaded model and scaler objects.

    Raises:
        FileNotFoundError: If the model or scaler file does not exist in the specified directory.
        pickle.UnpicklingError: If there is an error unpickling the model or scaler file.
    """
    with open(os.path.join(model_dir, 'model.pkl'), 'rb') as f:
        model = pickle.load(f)
        # model.set_params(n_jobs = 10)

    with open(os.path.join(model_dir, 'scaler.pkl'), 'rb') as f:
        scaler = pickle.load(f)
    return model, scaler

def _predict(df, input_data, model, ncores):
    """
    Predicts the probabilities using the given model and input data.

    Args:
        df (pandas.DataFrame): The DataFrame containing the original features.
        input_data (pandas.DataFrame): The DataFrame containing the scaled features.
        model (object): The machine learning model used for prediction.
        ncores (int): The number of cores to use for parallel processing.

    Returns:
        pandas.Series: A Series containing the predicted probabilities, indexed by the original DataFrame index.
    """
    scaled_features = ['{}_s'.format(c) for c in df.columns if c.startswith('f_')]
    df[scaled_features] = input_data

    nperchunk = len(df) / ncores
    df['chunk'] = (np.arange(len(df)) / nperchunk).astype(int)

    tasks = []
    for _, df_ in df.groupby('chunk'):
        tasks.append((model, df_[scaled_features].values.copy()))

    with Pool(ncores) as p:
        results = p.map(_predict_task, tasks)

    probs = np.concatenate(results)
    return pd.Series(data=probs, index=df.index)

def _predict_task(params):
    """
    Predicts the probability of the positive class for the given data using the provided model.

    Args:
        params (tuple): A tuple containing the model and the data to be predicted.
            - model: The trained model with a `predict_proba` method.
            - data: The input data for which the probabilities are to be predicted.

    Returns:
        numpy.ndarray: An array of probabilities for the positive class.
    """
    model, data = params
    return model.predict_proba(data)[:, 0]

def _write_output(df, labels, probs, output_dir):
    """
    Writes the classification results to a TSV file.

    Args:
        df (pandas.DataFrame): DataFrame containing the data to be classified.
        labels (list or numpy.array): List or array of classification labels.
        probs (list or numpy.array): List or array of classification probabilities.
        output_dir (str): Directory where the output file will be saved.

    Returns:
        None
    """
    df['temp'] = labels
    df['result'] = df['temp'].apply(lambda x: 'PASS' if x == 1 else "ARTIFACT")
    df['prob'] = probs
    df[['chrm', 'pos', 'result', 'prob']].to_csv(os.path.join(output_dir, 'result.tsv'), sep='\t', index=False)