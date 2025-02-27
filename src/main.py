"""
This file sets up the command line interface for various modules in the project.
Functions:
- main(): Entry point of the script. Parses arguments, sets up the environment, and runs the appropriate module.
- setup_module_arguments(parser): Configures the argument parser with subparsers for each mode.
- run_module(args): Imports and runs the appropriate module based on the provided mode.
- EF_add_parser_arguments(parser): Adds arguments for the 'extract_features' mode.
- AC_add_parser_arguments(parser): Adds arguments for the 'apply_classifier' mode.
- TC_add_parser_arguments(parser): Adds arguments for the 'train_classifier' mode.
- TGL_add_parser_arguments(parser): Adds arguments for the 'train_genlabels' mode.
- TPP_add_parser_arguments(parser): Adds arguments for the 'train_preprocessing' mode.
Usage:
    python main.py <mode> [options]
Modes:
    extract_features: Extracts features from input data.
    train_classifier: Trains a classifier using the provided training data.
    train_preprocessing: Preprocesses training data.
    train_genlabels: Generates labels for training data.
    apply_classifier: Applies a trained classifier to input data.
"""

import argparse
import os
from utils_setup import setup_ncores, setup_pandarallel

def main():
    parser = argparse.ArgumentParser()
    args = setup_module_arguments(parser)
    setup_ncores(args.cores)
    setup_pandarallel(args.mode, args.cores)
    run_module(args)

def setup_module_arguments(parser):
    modes_parser_setup = {
            "extract_features" : EF_add_parser_arguments,
            "train_classifier" : TC_add_parser_arguments,
            "train_preprocessing" : TPP_add_parser_arguments,
            "train_genlabels" : TGL_add_parser_arguments,
            "classify" : AC_add_parser_arguments
            }

    subparsers = parser.add_subparsers(dest = 'mode', required = True, \
            help='<Required> Select mode')

    for mode in modes_parser_setup.keys():
        subparser = subparsers.add_parser(mode)
        modes_parser_setup[mode](subparser)

    args = parser.parse_args()
    return args

def run_module(args):
    # Note that these modules are imported here and not earlier because ncores needs to be set 
    # before numpy is imported, and therefore before the modules are imported. 
    # Otherwise, numpy uses as many cores as it can because it's dumb like that

    mode = args.mode
    print('Running mode: {}'.format(mode))
    if mode == 'extract_features':
        import extract_features 
        maf, bams, mappability, output = args.input_file, args.bams, args.resources_dir, args.output
        extract_features.main(maf, bams, mappability, output)

    elif mode == 'train_classifier':
        import train_classifier 
        filelist, model_type, no_label_prop, output_dir = args.file_list, args.model, args.no_label_prop, args.output_dir
        train_classifier.main(filelist, model_type, no_label_prop, output_dir)

    elif mode == 'train_preprocessing':
        import train_preprocessing
        maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn = \
        args.maf, args.bam_dirs, args.signals_dir, args.output_dir, args.fullbam, args.cell_clones, args.hscn, args.use_cached_cn
        train_preprocessing.main(maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn)

    elif mode == 'train_genlabels':
        import train_gen_labels
        input_file, output_dir, clone1, clone2, alpha = args.input_file, args.output_dir, args.clone1, args.clone2, args.alpha
        train_gen_labels.main(input_file, output_dir, clone1, clone2, alpha)

    elif mode == 'classify':
        import apply_classifier
        model_dir, features, output_dir, chunksize, ncores= args.model_dir, args.features, args.output_dir, args.chunksize, args.cores
        apply_classifier.main(model_dir, features, output_dir, chunksize, ncores)

    else: 
        raise KeyError('Invalid mode provided: {} \nMode must be one of {}'.format(mode, str(modes_parser_setup.keys())))

def EF_add_parser_arguments(parser):
    # setup arguments for extract_features
    parser.add_argument(dest='input_file', type = str, help = '<Required> file containing candidate variants')
    parser.add_argument(dest='output', type = str, help = '<Required> Full path and name of output file')
    parser.add_argument(dest='bams', nargs="+", type = str, help = '<Required> list of bam files')
    parser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir))
    default_resources_path = os.path.join(repo_root, 'resources')
    parser.add_argument('--resources_dir', type=str, default=default_resources_path, help='<Optional> Path to directory containing folder of mappability tracks (default: {}'.format(default_resources_path))

def AC_add_parser_arguments(parser):
    # setup arguments for apply_classifier
    parser.add_argument(dest='features', type = str, help = '<Required> Input file containing variant features')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument(dest='model_dir', type = str, help = '<Required> Directory containing model.pkl and scaler.plk')
    parser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

    DEFAULT_CHUNKSIZE = 5000
    parser.add_argument('--chunksize', type = int, default = DEFAULT_CHUNKSIZE, required = False,
                        help = F'<Optional> Number of rows per worker (default {DEFAULT_CHUNKSIZE})')

def TC_add_parser_arguments(parser):
    # setup arguments for train_classifier
    parser.add_argument(dest='file_list', type = str, help = '<Required> file containing list of training data')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument('--model', type=str, default='gradientboosting', help = 'Model used for classification')
    parser.add_argument('--no_label_prop', action='store_true', help = 'Don\'t use label propagation to fill in missing labels')
    parser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

def TGL_add_parser_arguments(parser):
    # setup arguments for train_genlabels
    parser.add_argument(dest='input_file', type = str, help = '<Required> Input file containing variant counts by clone')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument(dest='clone1', type = str, help = '<Required> Name of clone')
    parser.add_argument(dest='clone2',  type = str, help = '<Required> Name of clone')
    parser.add_argument('--alpha', default=0.1, type=float, help='Per sample FDR. Total FDR = alpha^2', required=False)
    parser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

def TPP_add_parser_arguments(parser):
    # setup arguments for train_preprocessing
    parser.add_argument(dest='maf', type = str, help = '<Required> maf file containing candidate variants')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Full path and name of output file')
    parser.add_argument(dest='bam_dirs', nargs="+", type = str, help = '<Required> list of bam directories')
    parser.add_argument('--signals_dir', type = str, help = 'Directory of signals results')
    parser.add_argument('--cell_clones', type = str, help = 'Cell to clone mapping')
    parser.add_argument('--hscn', type = str, help='Signals HSCN file')
    parser.add_argument('--fullbam', action="store_true", help ='The list of bams is provided and not in region format')
    parser.add_argument('--use_cached_cn', action="store_true", help = 'Use already processed cell-to-clone map if it exists')
    parser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

if __name__ == '__main__':
    main()

    
