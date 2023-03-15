import argparse
import pandas as pd
import numpy as np
from sklearn.datasets import load_digits
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import pickle
import os
import numpy as np
from sklearn import datasets
from sklearn.semi_supervised import LabelSpreading


def main(args):
    validate_arguments(args)
    data, labels = read_input_data(filename)

    data_train, data_test, labels_train, labels_test = train_test_split(data, labels, test_size=0.20, random_state=42)

    # Standard scale data
    scaled_data_train, scaler = scale_data(data_train)
    scaled_data_test = scale_data(data_test, scaler)

    model = train_model(scaled_data_train, labels_train)
    test_model(scaled_data_test, labels_test, model, output_dir)
    write_model(model, scaler, output_dir)

def write_model(model, scaler, output_dir):
    pickle.dump(label_prop_model, open(os.path.join(output_dir, 'model.pkl'), 'wb'))
    pickle.dump(scaler, open(os.path.join(output_dir, 'scaler.pkl'), 'wb'))

def train_model(data, labels):
    label_prop_model = LabelSpreading(alpha=0.5)
    model = label_prop_model.fit(data, labels)
    return model

def test_model(data, labels, model, output_dir):
    inferred_labels = model.predict(data)
    accuracy, precision, recall = compute_performance_stats(new_labels, y_test)
    with open(os.path.join(output_dir, 'test_performance.txt'), 'w') as out:
        out.write("\t".join(['accuracy', 'precision', 'recall']) + "\n")
        out.write("\t".join(map(str, [accuracy, precision, recall])))

def scale_data(data, scaler = None):
    if scaler:
        return scaler.transform(data)
    else:
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)
        return data, scaler

def compute_performance_stats(inferred, groundtruth):
    tp = sum([1 for inf, gt in zip(inferred, groundtruth) if inf == 1 and gt == 1])
    fp = sum([1 for inf, gt in zip(inferred, groundtruth) if inf == 1 and gt == 0])
    tn = sum([1 for inf, gt in zip(inferred, groundtruth) if inf == 0 and gt == 0])
    fn = sum([1 for inf, gt in zip(inferred, groundtruth) if inf == 0 and gt == 1])

    accuracy = (tp+tn)/(tp+tn+fp+fn)
    precision = (tp)/(tp+fp)
    recall = (tp)/(tp+fn)

    return accuracy, precision, recall

def add_parser_arguments(parser):
    parser.add_argument(dest='file_list', type = str, help = '<Required> file containing list of training data')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')

def validate_arguments(args):
    # Checks if input files exist and if output files are in directories that exist and can be written to
    for arg in vars(args):
        print(arg, getattr(args, arg))

    assert os.path.isfile(args.file_list)

def read_input_data(filename):
    """
    Filename contains a list of pairs of files, one for each sample or individual: (feature_file, label_file)
    These files should be output from extract_features and mixture_model respectively
    """

    dfs = []
    with open(filename) as f:
        for i, line in enumerate(f):
            feature_file, label_file = line.strip().split(' ')
            features = pd.read_table(feature_file)
            labels = pd.read_table(label_file)[['chrm', 'pos', 'assignment', 'responsibilities']]
            df = features.merge(labels, on = ['chrm', 'pos'], how = 'inner')
            df['sample'] = i
            dfs.append(df)
    df = pd.concat(dfs).dropna()

    feature_cols = [c for c in df.columns if c.startswith('f_')]
    data = df[feature_cols].values
    labels = df.assignment.copy()
    labels[(df.assignment == 1) | (df.assignment == 2)] = -1 # labels[(df.assignment == 1) | (df.assignment == 2)].apply(lambda x: -1 if random() < 0.9 else 0)
    labels[df.assignment == 0 ] = 1
    labels[df.assignment == 3] = 0

    return data, labels

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    main(args)
