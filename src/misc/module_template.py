import argparse

def main(args):
    pass

def add_parser_arguments(parser):
    pass

def validate_arguments(args):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_parser_arguments(parser)
    args = parser.parse_args()

    main(args)
