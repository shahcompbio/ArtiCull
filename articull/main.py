"""
main.py

This module serves as the commandline entry point for running various modes of the application,
including feature extraction, classifier training, preprocessing, label generation, and classification.
It sets up the computational environment, parses command-line arguments, 
and dispatches tasks to the appropriate modules.

Functions:
    main():
        The main function that sets up the argument parser, initializes the computational environment, and runs the specified module.
"""

import argparse
import os
from .utils.utils_setup import setup_ncores, setup_pandarallel

def main():
    """
    The main function that sets up the argument parser, initializes the computational environment, and runs the specified module.

    Args:
        None

    Returns:
        None
    """
    parser = argparse.ArgumentParser()
    args = setup_module_arguments(parser)
    setup_ncores(args.cores)
    setup_pandarallel(args.mode, args.cores)
    run_module(args)

def setup_module_arguments(parser):
    """
    Sets up the argument parser with subparsers for each mode of operation.

    Args:
        parser (argparse.ArgumentParser): The argument parser to set up.

    Returns:
        argparse.Namespace: The parsed command-line arguments.
    """
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
    """
    Runs the specified module based on the parsed command-line arguments.

    Args:
        args (argparse.Namespace): The parsed command-line arguments.

    Raises:
        KeyError: If an invalid mode is provided.

    Returns:
        None
    """
    # Note that these modules are imported here and not earlier because ncores needs to be set 
    # before numpy is imported, and therefore before the modules are imported. 
    # Otherwise, numpy uses as many cores as it can because it's dumb like that

    mode = args.mode
    print('Running mode: {}'.format(mode))
    if mode == 'extract_features':
        from classify import extract_features 
        maf, bams, mappability, output = args.input_file, args.bams, args.resources_dir, args.output
        extract_features.extract_features(maf, bams, mappability, output)

    elif mode == 'train_classifier':
        from train import train_classifier 
        filelist, model_type, no_label_prop, output_dir = args.file_list, args.model, args.no_label_prop, args.output_dir
        train_classifier.train_classifier(filelist, model_type, no_label_prop, output_dir)

    elif mode == 'train_preprocessing':
        from train import preprocessing
        maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn = \
        args.maf, args.bam_dirs, args.signals_dir, args.output_dir, args.fullbam, args.cell_clones, args.hscn, args.use_cached_cn
        preprocessing.preprocessing(maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn)

    elif mode == 'train_genlabels':
        from train import generate_labels
        input_file, output_dir, clone1, clone2, alpha = args.input_file, args.output_dir, args.clone1, args.clone2, args.alpha
        generate_labels.generate_labels(input_file, output_dir, clone1, clone2, alpha)

    elif mode == 'classify':
        from classify import apply_classifier
        model_dir, features, output_dir, chunksize, ncores= args.model_dir, args.features, args.output_dir, args.chunksize, args.cores
        apply_classifier.classify(model_dir, features, output_dir, chunksize, ncores)

    else: 
        raise KeyError('Invalid mode provided: {} \nMode must be one of {}'.format(mode, str(modes_parser_setup.keys())))

def EF_add_parser_arguments(parser):
    """
    Adds command-line arguments for the 'extract_features' mode.

    Args:
        parser (argparse.ArgumentParser): The argument parser to set up.

    Returns:
        None
    """
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
    """
    Adds command-line arguments for the 'classify' mode.

    Args:
        parser (argparse.ArgumentParser): The argument parser to set up.

    Returns:
        None
    """
    parser.add_argument(dest='features', type = str, help = '<Required> Input file containing variant features')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument(dest='model_dir', type = str, help = '<Required> Directory containing model.pkl and scaler.plk')
    parser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

    DEFAULT_CHUNKSIZE = 5000
    parser.add_argument('--chunksize', type = int, default = DEFAULT_CHUNKSIZE, required = False,
                        help = F'<Optional> Number of rows per worker (default {DEFAULT_CHUNKSIZE})')

def TC_add_parser_arguments(parser):
    """
    Adds command-line arguments for the 'train_classifier' mode.

    Args:
        parser (argparse.ArgumentParser): The argument parser to set up.

    Returns:
        None
    """
    parser.add_argument(dest='file_list', type = str, help = '<Required> file containing list of training data')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument('--model', type=str, default='gradientboosting', help = 'Model used for classification')
    parser.add_argument('--no_label_prop', action='store_true', help = 'Don\'t use label propagation to fill in missing labels')
    parser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

def TGL_add_parser_arguments(parser):
    """
    Adds command-line arguments for the 'train_genlabels' mode.

    Args:
        parser (argparse.ArgumentParser): The argument parser to set up.

    Returns:
        None
    """
    parser.add_argument(dest='input_file', type = str, help = '<Required> Input file containing variant counts by clone')
    parser.add_argument(dest='output_dir', type = str, help = '<Required> Output directory')
    parser.add_argument(dest='clone1', type = str, help = '<Required> Name of clone')
    parser.add_argument(dest='clone2',  type = str, help = '<Required> Name of clone')
    parser.add_argument('--alpha', default=0.1, type=float, help='Per sample FDR. Total FDR = alpha^2', required=False)
    parser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

def TPP_add_parser_arguments(parser):
    """
    Adds command-line arguments for the 'train_preprocessing' mode.

    Args:
        parser (argparse.ArgumentParser): The argument parser to set up.

    Returns:
        None
    """
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