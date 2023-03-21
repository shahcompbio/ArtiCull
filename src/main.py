import argparse
import os

def setup_module_arguments(parser):
    import mixture_model_preprocessing
    import extract_features
    import mixture_model
    import train_classifier
    import apply_classifier

    modes = {
            "extract_features" : extract_features,
            "mm_preprocessing" : mixture_model_preprocessing,
            "mixture_model" : mixture_model,
            "train" : train_classifier,
            "classify": apply_classifier
            }

    subparsers = parser.add_subparsers(dest = 'mode', required = True, \
            help='<Required> Select mode')

    for mode in modes:
        subparser = subparsers.add_parser(mode)
        modes[mode].add_parser_arguments(subparser)
        #if mode not in ['train_classifier', 'classify']:
        subparser.add_argument('--cores', '-j', default = None, type = int, \
            help = 'Number of workers to use for parallelization. <Default> the number of available cores')

    args = parser.parse_args()
    mode = args.mode

    if mode not in ['train_classifier', 'classify']:
        from pandarallel import pandarallel
        progress_bar = mode != 'classify'
        if args.cores:
            pandarallel.initialize(progress_bar=progress_bar, nb_workers = args.cores)
        else:
            pandarallel.initialize(progress_bar=progress_bar) #, nb_workers = args.cores)

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
