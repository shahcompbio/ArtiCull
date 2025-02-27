import argparse
import os

def setup_module_arguments(parser):
    from . import train_preprocessing
    from . import extract_features
    from . import train_gen_labels
    from . import train_classifier
    from . import apply_classifier

    modes = {
            "extract_features" : extract_features,
            "train_preprocessing" : train_preprocessing,
            "train_genlabels" : train_gen_labels,
            "train" : train_classifier,
            "classify": apply_classifier
            }

    subparsers = parser.add_subparsers(dest = 'mode', required = True, \
            help='<Required> Select mode')

    for mode in modes:
        subparser = subparsers.add_parser(mode)
        modes[mode].add_parser_arguments(subparser)
        subparser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

    args = parser.parse_args()
    mode = args.mode

    from pandarallel import pandarallel # type: ignore
    progress_bar = mode != 'classify'
    if args.cores:
        pandarallel.initialize(progress_bar=progress_bar, nb_workers = args.cores, use_memory_fs=True)
    else:
        pandarallel.initialize(progress_bar=progress_bar, use_memory_fs=False) #, nb_workers = args.cores)

    try:
        modes[mode]
    except KeyError:
        raise KeyError('Invalid mode provided: {} \nMode must be one of {}'.format(mode, str(modes.keys())))

    modes[mode].main(args)


def setup_ncores(parser):
    parser.add_argument('--cores', '-j', default = None, type = int, \
        help = 'Number of workers to use for parallelization. <Default> the number of available cores')
    args, arglist = parser.parse_known_args()

    import psutil
    threads_per_core = psutil.cpu_count() / psutil.cpu_count(logical=False)

    if args.cores:
        os.environ['OPENBLAS_NUM_THREADS'] = str(args.cores * threads_per_core)
        os.environ['MKL_NUM_THREADS'] = str(args.cores * threads_per_core)

    return parser, args

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    setup_ncores(parser)
    setup_module_arguments(parser)
