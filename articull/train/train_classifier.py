"""
train_classifier.py

This module provides functions to train, test, and save machine learning classifiers using various models.
It includes functionality for reading input data, scaling data, training models, testing models, and saving the trained models and scalers.

Functions:
    train_classifier(filelist, model_type, no_label_prop, output_dir):
        Trains a classifier model using the provided data and parameters.
"""

import pandas as pd # type: ignore
from sklearn.model_selection import train_test_split # type: ignore
from sklearn.preprocessing import StandardScaler # type: ignore
import pickle
import os
from sklearn.semi_supervised import LabelSpreading # type: ignore
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier # type: ignore
from sklearn.svm import LinearSVC # type: ignore
from sklearn.linear_model import LogisticRegression # type: ignore
from sklearn.neural_network import MLPClassifier # type: ignore
from sklearn.calibration import CalibratedClassifierCV # type: ignore

def train_classifier(filelist, model_type, no_label_prop, output_dir):
    """
    Trains a classifier model using the provided data and parameters.

    Args:
        filelist (list of str): List of file paths containing the input data.
        model_type (str): Type of model to train. Options are 'randomforest', 'logistic', 
            'linearsvc', 'gradientboosting', and 'mlp'. Default is 'gradientboosting'.
        no_label_prop (bool): If True, indicates that the model will not perform label propagation.
        output_dir (str): Directory where the trained model and scaler will be saved.

    Returns:
        None

    Raises:
        ValueError: If the provided arguments are not valid.

    This function performs the following steps:
        1. Validates the input arguments.
        2. Reads the input data and labels from the provided file list.
        3. Splits the data into training and testing sets.
        4. Scales the training data and applies the same scaling to the test data.
        5. Trains the specified model using the scaled training data and labels.
        6. Tests the trained model using the scaled test data and labels.
        7. Saves the trained model and scaler to the specified output directory.
    """
    _validate_arguments(filelist, output_dir)

    data, labels = _read_input_data(filelist)
    data_train, data_test, labels_train, labels_test = train_test_split(data, labels, test_size=0.20, random_state=42)

    # Standard scale data
    scaled_data_train, scaler = _scale_data(data_train)
    scaled_data_test = _scale_data(data_test, scaler)

    model = _train_model(scaled_data_train, labels_train, model_type, no_label_prop=no_label_prop)
    _test_model(scaled_data_test, labels_test, model, output_dir)
    _write_model(model, scaler, output_dir)

def _downsample_unlabeled(data, n=1000):
    """
    Downsamples the given DataFrame to a specified number of samples.

    Args:
        data (pandas.DataFrame): The DataFrame to downsample.
        n (int, optional): The number of samples to return. Default is 1000.

    Returns:
        pandas.DataFrame: A DataFrame with `n` samples if possible, otherwise returns the original DataFrame.
    """
    try:
        return data.sample(n, random_state=0)
    except ValueError:
        return data  # if not big enough

def _write_model(model, scaler, output_dir):
    """
    Saves the trained model and scaler to the specified output directory.

    Args:
        model (object): The trained model to be saved.
        scaler (object): The scaler used for preprocessing the data.
        output_dir (str): The directory where the model and scaler will be saved.

    Returns:
        None
    """
    pickle.dump(model, open(os.path.join(output_dir, 'model.pkl'), 'wb'))
    pickle.dump(scaler, open(os.path.join(output_dir, 'scaler.pkl'), 'wb'))

def _train_model(data, labels, model_type='gradientboosting', no_label_prop=False):
    """
    Trains a machine learning model based on the specified model type and data.

    Args:
        data (array-like or sparse matrix): The input data to train the model.
        labels (array-like): The target labels corresponding to the input data.
        model_type (str, optional): The type of model to train. Options are 'randomforest', 'logistic', 
            'linearsvc', 'gradientboosting', and 'mlp'. Default is 'gradientboosting'.
        no_label_prop (bool, optional): If True, the model is trained directly on the provided labels. 
            If False, label propagation is used to infer labels before training. Default is False.

    Returns:
        model: The trained machine learning model.

    Raises:
        Exception: If the specified model_type is not recognized.
    """
    if model_type.lower() == 'randomforest':
        model = RandomForestClassifier(max_depth=2, random_state=0)
    elif model_type.lower() == 'logistic':
        model = LogisticRegression(random_state=0)
    elif model_type.lower() == 'linearsvc':
        svm = LinearSVC(random_state=0)
        model = CalibratedClassifierCV(svm)  # To get probabilities
    elif model_type.lower() == 'gradientboosting':
        model = HistGradientBoostingClassifier(random_state=0)
    elif model_type.lower() == 'mlp':
        model = MLPClassifier(random_state=0)
    else:
        raise Exception('Model not recognized: {}'.format(model))

    if no_label_prop:
        model = model.fit(data, labels)
    else:
        label_prop_model = LabelSpreading(alpha=0.05)
        label_prop_model = label_prop_model.fit(data, labels)
        full_labels = label_prop_model.transduction_
        model = model.fit(data, full_labels)

    return model

def _test_model(data, labels, model, output_dir):
    """
    Tests the given model on the provided data and labels, and writes the performance statistics to a file.

    Args:
        data (array-like): The input data to test the model on.
        labels (array-like): The true labels corresponding to the input data.
        model (object): The trained model to be tested. Must have a `predict` method.
        output_dir (str): The directory where the performance statistics file will be saved.

    Returns:
        None
    """
    inferred_labels = model.predict(data)
    accuracy, precision, recall = _compute_performance_stats(inferred_labels, labels)
    with open(os.path.join(output_dir, 'test_performance.txt'), 'w') as out:
        out.write("\t".join(['accuracy', 'precision', 'recall']) + "\n")
        out.write("\t".join(map(str, [accuracy, precision, recall])))

def _scale_data(data, scaler=None):
    """
    Scales the input data using the provided scaler or a new StandardScaler.

    Args:
        data (array-like): The data to be scaled.
        scaler (object, optional): An instance of a scaler with a `transform` method. 
            If None, a new StandardScaler will be used. Default is None.

    Returns:
        tuple: If a new StandardScaler is used, returns a tuple containing the scaled data and the scaler.
        array-like: If a scaler is provided, returns the scaled data.
    """
    if scaler:
        return scaler.transform(data)
    else:
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)
        return scaled_data, scaler

def _compute_performance_stats(inferred, groundtruth):
    """
    Computes performance statistics for a binary classifier.

    Args:
        inferred (list of int): List of inferred labels (0 or 1) from the classifier.
        groundtruth (list of int): List of ground truth labels (0 or 1).

    Returns:
        tuple: A tuple containing:
            - accuracy (float): The accuracy of the classifier.
            - precision (float): The precision of the classifier.
            - recall (float): The recall of the classifier.
    """
    tp = sum([1 for inf, gt in zip(inferred, groundtruth) if inf == 1 and gt == 1])
    fp = sum([1 for inf, gt in zip(inferred, groundtruth) if inf == 1 and gt == 0])
    tn = sum([1 for inf, gt in zip(inferred, groundtruth) if inf == 0 and gt == 0])
    fn = sum([1 for inf, gt in zip(inferred, groundtruth) if inf == 0 and gt == 1])

    accuracy = (tp + tn) / (tp + tn + fp + fn)
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    return accuracy, precision, recall

def _validate_arguments(filelist, output_dir):
    """
    Validates the input file and output directory.

    Args:
        filelist (str): Path to the input file.
        output_dir (str): Path to the output directory.

    Raises:
        AssertionError: If the input file does not exist or cannot be read.
        AssertionError: If the output directory does not exist or cannot be written to.
    """
    assert os.path.isfile(filelist), f"Input file {filelist} does not exist."
    assert os.access(filelist, os.R_OK), f"Input file exists, but cannot be read due to permissions: {filelist}"
    assert os.path.isdir(output_dir), f"Output directory {output_dir} does not exist."
    assert os.access(output_dir, os.W_OK), f"Output directory exists, but cannot be written to due to permissions: {output_dir}"

def _read_input_data(filename):
    """
    Reads input data from a file containing pairs of feature and label files for each sample.


    The input file should have the following format:
        feature_file_path_1\tlabel_file_path_1
        feature_file_path_2\tlabel_file_path_2
        ...

    Each feature file and label file should be tab-separated files. The label file should contain columns 'chrm', 'pos', and 'assignment'.
    The function merges the feature and label data on 'chrm' and 'pos', filters for specific variant types ('SNP', 'DNP'), and handles missing values in 'f_p_normal'.
    It also downsamples unlabeled data and concatenates all samples into a single DataFrame before extracting feature columns and labels.

    Args:
        filename (str): Path to the file containing pairs of feature and label file paths, separated by a tab.

    Returns:
        tuple: A tuple containing:
            - data (numpy.ndarray): Array of feature values.
            - labels (pandas.Series): Series of label assignments.
    """
    dfs = []
    with open(filename) as f:
        for i, line in enumerate(f):
            feature_file, label_file = line.strip().split('\t')
            features = pd.read_table(feature_file, sep='\t')
            labels = pd.read_table(label_file, sep='\t')
            labels = labels[['chrm', 'pos', 'assignment']]  # x, 'responsibilities']]
            df = features.merge(labels, on=['chrm', 'pos'], how='inner')
            df['sample'] = i
            df = df[df['var_type'].isin(['SNP', 'DNP'])]
            df['f_p_normal'] = df['f_p_normal'].fillna(0)

            df_labeled = df[df['assignment'] != -1]
            df_unlabeled = _downsample_unlabeled(df[df['assignment'] == -1])

            dfs.append(df_labeled)
            dfs.append(df_unlabeled)

    df = pd.concat(dfs).dropna()

    feature_cols = [c for c in df.columns if c.startswith('f_')]
    data = df[feature_cols].values
    labels = df.assignment.copy()

    return data, labels