import argparse

import mixture_model_preprocessing
import extract_features
import mixture_model
import train_classifier
import apply_classifier

from pandarallel import pandarallel


if __name__ == '__main__':
    modes = {
            "extract_features" : extract_features,
            "mm_preprocessing" : mixture_model_preprocessing,
            "mixture_model" : mixture_model,
            "train" : train_classifier,
            "classify": apply_classifier
            }

    parser = argparse.ArgumentParser()
     # to use. By default, will use all available ')
    subparsers = parser.add_subparsers(dest = 'mode', required = True, \
            help='<Required> Select mode')


    for mode in modes:
        subparser = subparsers.add_parser(mode)
        modes[mode].add_parser_arguments(subparser)
        if mode not in ['train_classifier', 'classify']:
            subparser.add_argument('--cores', '-j', default = None, type = int, \
                help = 'Number of workers to use for parallelization. <Default> the number of available cores')


    args = parser.parse_args()
    if args.cores:
        pandarallel.initialize(progress_bar=True, nb_workers = args.cores)
    else:
        pandarallel.initialize(progress_bar=True) #, nb_workers = args.cores)

    mode = args.mode
    try:
        modes[mode]
    except KeyError:
        raise Exception('Invalid mode provided: {} \nMode must be one of {}'.format(mode, str(modes.keys())))

    modes[mode].main(args)
