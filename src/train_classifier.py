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
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.svm import LinearSVC
from sklearn.linear_model import LogisticRegression
from sklearn.neural_network import MLPClassifier
from sklearn.calibration import CalibratedClassifierCV

def main(args):
    validate_arguments(args)
    data, labels = read_input_data(args.file_list)
    model_type = args.model
    no_label_prop = args.no_label_prop

    data_train, data_test, labels_train, labels_test = train_test_split(data, labels, test_size=0.20, random_state=42)

    print(len(data))
    # Standard scale data
    scaled_data_train, scaler = scale_data(data_train)
    print(scaled_data_train.sum())
    scaled_data_test = scale_data(data_test, scaler)
    print(scaled_data_test.sum())

    model = train_model(scaled_data_train, labels_train, model_type, no_label_prop=no_label_prop)
    test_model(scaled_data_test, labels_test, model, args.output_dir)
    write_model(model, scaler, args.output_dir)

def downsample_unlabeled(data, n=1000):
    try:
        return data.sample(n, random_state = 0)
    except ValueError:
        return data # if not big enough

def write_model(model, scaler, output_dir):
    pickle.dump(model, open(os.path.join(output_dir, 'model.pkl'), 'wb'))
    pickle.dump(scaler, open(os.path.join(output_dir, 'scaler.pkl'), 'wb'))

def train_model(data, labels, model_type='gradientboosting', no_label_prop=False):
    if model_type.lower() == 'randomforest':
        model = RandomForestClassifier(max_depth=2, random_state=0)
    elif model_type.lower() == 'logistic':
        model = LogisticRegression(random_state=0)
    elif model_type.lower() == 'linearsvc':
        svm = LinearSVC(random_state = 0)
        model = CalibratedClassifierCV(svm) # To get probabilities
    elif model_type.lower() == 'gradientboosting':
        model = HistGradientBoostingClassifier(random_state = 0)
                #class_weight={0:0.9, 1:.1})
    elif model_type.lower() == 'mlp':
        model = MLPClassifier(random_state = 0)
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

def test_model(data, labels, model, output_dir):
    inferred_labels = model.predict(data)
    accuracy, precision, recall = compute_performance_stats(inferred_labels, labels)
    with open(os.path.join(output_dir, 'test_performance.txt'), 'w') as out:
        out.write("\t".join(['accuracy', 'precision', 'recall']) + "\n")
        out.write("\t".join(map(str, [accuracy, precision, recall])))

def scale_data(data, scaler = None):
    if scaler:
        print("scaling here")
        return scaler.transform(data)
    else:
        print("fit scaling")
        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)
        return scaled_data, scaler

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
    parser.add_argument('--model', type=str, default='gradientboosting', help = 'Model used for classification')
    parser.add_argument('--no_label_prop', action='store_true', help = 'Don\'t use label propagation to fill in missing labels')

def validate_arguments(args):
    # Checks if input files exist and if output files are in directories that exist and can be written to
    for arg in vars(args):
        print(arg, ':\t', getattr(args, arg))

    assert os.path.isfile(args.file_list)

def read_input_data(filename):
    """
    Filename contains a list of pairs of files, one for each sample or individual: (feature_file, label_file)
    These files should be output from extract_features and mixture_model respectively
    """

    dfs = []
    with open(filename) as f:
        for i, line in enumerate(f):
            feature_file, label_file = line.strip().split('\t')
            features = pd.read_table(feature_file, sep='\t')
            labels = pd.read_table(label_file, sep='\t')
            labels = labels[['chrm', 'pos', 'assignment']] #x, 'responsibilities']]
            df = features.merge(labels, on = ['chrm', 'pos'], how = 'inner')
            df['sample'] = i
            df = df[df['var_type'].isin(['SNP', 'DNP'])]
            df['f_p_normal'] = df['f_p_normal'].fillna(0)

            df_labeled = df[df['assignment'] != -1]
            df_unlabeled = downsample_unlabeled(df[df['assignment'] == -1])


            dfs.append(df_labeled)
            dfs.append(df_unlabeled)


    df = pd.concat(dfs).dropna()

    feature_cols = [c for c in df.columns if c.startswith('f_')]
    data = df[feature_cols].values
    labels = df.assignment.copy()

    return data, labels

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    main(args)
