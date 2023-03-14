import argparse

import mixture_model_preprocessing
import extract_features
import mixture_model

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest = 'mode', required = True, \
        help='<Required> Running mode in {"extract_features", "mm_preprocessing", "mixture_model"} specifies functionality')

    a_parser = subparsers.add_parser("extract_features")
    extract_features.add_parser_arguments(a_parser)

    b_parser = subparsers.add_parser("mm_preprocessing")
    mixture_model_preprocessing.add_parser_arguments(b_parser)

    c_parser = subparsers.add_parser("mixture_model")
    mixture_model.add_parser_arguments(c_parser)

    args = parser.parse_args()

    mode = args.mode
    if mode == "extract_features":
        extract_features.main(args)
    elif mode == "mm_preprocessing":
        mixture_model_preprocessing.main(args)
    elif mode == "mixture_model":
        mixture_model.main(args)
    else:
        print('Invalid mode provided: {} \nMode must be in {"extract_features", "mm_preprocessing", "mixture_model"}')
