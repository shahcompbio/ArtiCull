"""
This module is the commandline entry point for running various modes of the application. 
It sets up the computational environment, parses command-line arguments, and dispatches tasks to the specified module.
For instructions on using the command-line interface, run `python articull/articull.py -h`, or see
instructions in :doc:`usage`

Functions:
    main():
        The main function that sets up the argument parser, initializes the computational environment, and runs the specified module.
"""

import argparse
import os
import sys
from articull._utils.setup import setup_ncores, setup_pandarallel

def articull():
    """
    The main function that sets up the argument parser, initializes the computational environment, and runs the specified module.

    Args:
        None

    Returns:
        None
    """
    parser = argparse.ArgumentParser()
    args = _setup_module_arguments(parser)

    setup(args.cores, args.progress_bar)
    _run_module(args)

def setup(ncores, progress_bar):
    setup_ncores(ncores)
    setup_pandarallel(progress_bar, ncores)

def download_resources(output_dir, force=False):
    """
    Downloads the genomic tracks required for running ArtiCull.

    Args:
        output_path (str): The path to the directory where the resources should be downloaded.
        force (bool): If True, the resources will be downloaded even if they already exist. <Default> False

    Returns:
        None
    """

    if not os.path.exists(output_dir):
        raise FileNotFoundError(f"Directory {output_dir} does not exist.")
    if not os.path.isdir(output_dir):
        raise NotADirectoryError(f"Path {output_dir} is not a directory.")
    if not os.access(output_dir, os.W_OK):
        raise PermissionError(f"Directory {output_dir} is not writable.")
    
    required_files = [f"chr{chrom}.bedGraph" for chrom in list(range(1, 23)) + ['X', 'Y']]
    missing_files = [os.path.join(output_dir, 'mappability', f) for f in required_files if not os.path.isfile(os.path.join(output_dir, 'mappability', f))]

    if len(missing_files) == 0 and not force:
        print(f"Mappability files already exist in {output_dir}. Skipping download and setup. Use parameter force=True to force re-download if files are incomplete or corrupted.")
        return
    else: 
        print("TEST")
        return 

    script_path = os.path.join(os.path.dirname(__file__), 'setup_mappability_track.bash')
    os.system(f'bash {script_path} {output_dir}')



def _setup_module_arguments(parser):
    """
    Sets up the argument parser with subparsers for each mode of operation.

    Args:
        parser (argparse.ArgumentParser): The argument parser to set up.

    Returns:
        argparse.Namespace: The parsed command-line arguments.
    """
    modes_parser_setup = {
            "extract_features" : _EF_add_parser_arguments,
            "train_classifier" : _TC_add_parser_arguments,
            "train_preprocessing" : _TPP_add_parser_arguments,
            "train_genlabels" : _TGL_add_parser_arguments,
            "classify_variants" : _AC_add_parser_arguments,
            "classify" : _C_add_parser_arguments
            }

    subparsers = parser.add_subparsers(dest = 'mode', required = True, \
            help='<Required> Select mode')


    for mode in modes_parser_setup.keys():
        subparser = subparsers.add_parser(mode)
        modes_parser_setup[mode](subparser)

    prog_bar_default = sys.stdout.isatty()
    parser.add_argument('--progress_bar', type=int, default=prog_bar_default, required=False,
                    help='<Optional> Show dynamic progress bar. Default is true if output is to stdout and false if output has been redirected.')

    args = parser.parse_args()
    return args

def _run_module(args):
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
        from articull.classify import extract_features 
        maf, bams, mappability, output = args.input_file, args.bams, args.resources_dir, args.output
        extract_features.extract_features(maf, bams, mappability, output)

    elif mode == 'train_classifier':
        from articull.train import train_classifier 
        filelist, model_type, no_label_prop, output_dir = args.file_list, args.model, args.no_label_prop, args.output_dir
        train_classifier.train_classifier(filelist, model_type, no_label_prop, output_dir)

    elif mode == 'train_preprocessing':
        from articull.train import preprocessing
        maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn, filter_vcf = \
        args.maf, args.bam_dirs, args.signals_dir, args.output_dir, args.fullbam, args.cell_clones, args.hscn, args.use_cached_cn, not args.no_vcf_filter
        preprocessing.preprocessing(maf, bam_dirs, signals_dir, output_dir, fullbam, cellclone_file, hscn_file, use_cached_cn, filter_vcf=filter_vcf)

    elif mode == 'train_genlabels':
        from articull.train import generate_labels
        input_file, output_dir, clone1, clone2, alpha = args.input_file, args.output_dir, args.clone1, args.clone2, args.alpha
        generate_labels.generate_labels(input_file, output_dir, clone1, clone2, alpha)

    elif mode == 'classify_variants':
        from articull.classify import classify

        model_dir, features, output_dir, chunksize, ncores, filter_vcf = args.model_dir, args.features, args.output_dir, args.chunksize, args.cores, not args.no_vcf_filter
        classify.classify_variants(model_dir, features, output_dir, chunksize, ncores, filter_vcf=filter_vcf)

    elif mode == 'classify':
        from articull.classify import extract_features, classify

        if args.features_file:
            features_file = args.features_file
        if not args.features_file:
            maf, bams, mappability, output, filter_vcf = args.input_file, args.bams, args.resources_dir, args.output_prefix, not args.no_vcf_filter
            features_file = extract_features.extract_features(maf, bams, mappability, output, filter_vcf=filter_vcf)

        if not args.extract_features_only: 
            model_dir, output_prefix, chunksize, ncores = args.model_dir, args.output_prefix, args.chunksize, args.cores
            classify.classify_variants(model_dir, features_file, output_prefix, chunksize, ncores)


def _EF_add_parser_arguments(parser):
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
    parser.add_argument('--no_vcf_filter', action='store_true', default=False,
                        help='<Optional> If set, do not filter VCF for FILTER == PASS. Default is to apply the filter.')
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir))
    default_resources_path = os.path.join(repo_root, 'resources')
    parser.add_argument('--resources_dir', type=str, default=default_resources_path, help='<Optional> Path to directory containing folder of mappability tracks (default: {}'.format(default_resources_path))

def _C_add_parser_arguments(parser):

    parser.add_argument(dest='input_file', type=str, help='<Required> maf or vcf containing candidate variants')
    parser.add_argument(dest='output_prefix', type=str,
                help='<Required> Output prefix (directory and sample name, e.g. /path/to/sample1. '
                    'Output files will be saved as /path/to/sample1_features.tsv and /path/to/sample1_result.tsv)')
    parser.add_argument(dest='model_dir', type=str, help='<Required> Directory containing model.pkl and scaler.plk')

    parser.add_argument(dest='bams', nargs="+", type=str, help='<Required> list of bam files')

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.abspath(os.path.join(script_dir, os.pardir))
    default_resources_path = os.path.join(repo_root, 'resources')
    parser.add_argument('--resources_dir', type=str, default=default_resources_path,
                        help='<Optional> Path to directory containing folder of mappability tracks (default: {})'.format(default_resources_path))


    DEFAULT_CHUNKSIZE = 5000
    parser.add_argument('--chunksize', type=int, default=DEFAULT_CHUNKSIZE, required=False,
                        help='<Optional> Number of rows per worker (default {})'.format(DEFAULT_CHUNKSIZE))

    parser.add_argument('--features_file', type=str, default=None,
                        help='<Optional> File containing features (e.g., generated by a previous run of articull). '
                            'If not provided, features will be extracted from input bam files')
    parser.add_argument('--cores', '-j', default=None, type=int,
                        help='Number of workers to use for parallelization. <Default> the number of available cores')
    parser.add_argument('--no_vcf_filter', action='store_true', default=False,
                        help='<Optional> If set, do not filter VCF for FILTER == PASS. Default is to apply the filter.')
    parser.add_argument('--extract_features_only', action='store_true', default=False,
                        help='<Optional> If set, will only extract features, but not run classification. ')


def _AC_add_parser_arguments(parser):
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

def _TC_add_parser_arguments(parser):
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

def _TGL_add_parser_arguments(parser):
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

def _TPP_add_parser_arguments(parser):
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
    parser.add_argument('--no_vcf_filter', action='store_true', default=False,
                        help='<Optional> If set, do not filter VCF for FILTER == PASS. Default is to apply the filter.')

if __name__ == '__main__':
    articull()