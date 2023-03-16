import argparse
import pickle
import pandas as pd
import sys

def main(args):
    validate_arguments(args)
    print("1. Loading model from: {}".format(args.model_dir))
    model, scaler = load_model(args.model_dir)
    print("2. Loading data from: {}".format(args.input_file))
    df, input_data = get_input_data(args.input_file, scaler)
    print("3. Predicting")
    labels, probs = predict(input_data, model)
    print("4. Writing results to: {}".format(args.output_dir))
    write_output(df, labels, probs, args.output_dir)

def add_parser_arguments(parser):
    parser.add_argument(dest='features', type = str, help = '<Required> Input file containing variant features')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument(dest='model_dir', type = str, help = '<Required> Directory containing model.pkl and scaler.plk')

def validate_arguments(args):
    pass

def load_model(model_dir):
    with open(os.path.join(model_dir, 'model.pkl'), 'rb') as f:
        model = pickle.load(f)

    with open(os.path.join(model_dir, 'scaler.pkl'), 'rb') as f:
        scaler = pickle.load(f)
    return model, scaler

def get_input_data(input_file, scaler):
    df = pd.read_table(input_file)

    features = [c for c in df.columns if c.starts_with('f_')]
    data = df[features].values

    scaled_data = scaler.transform(data)
    return df, scaled_data

def predict(input_data, model):
    labels = model.predict(input_data)
    probs = model.predict_proba(input_data)
    return labels, probs

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
