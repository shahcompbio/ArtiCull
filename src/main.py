import argparse

import mixture_model_preprocessing
import extract_features
import mixture_model
import train_classifier

if __name__ == '__main__':
    modes = {
            "extract_features" : extract_features,
            "mm_preprocessing" : mixture_model_preprocessing,
            "mixture_model" : mixture_model,
            "train_classifier" : train_model
            }

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest = 'mode', required = True, \
        help='<Required> Running mode: {}'.format(str(modes.keys()))

    for mode in modes:
        subparser = subparsers.add_parser(mode)
        modes[mode].add_parser_arguments(subparser)

    args = parser.parse_args()
    mode = args.mode
    try:
        modes[mode].main(args)
    except:
        print('Invalid mode provided: {} \nMode must be one of {}'.format(mode, str(modes.keys())))
