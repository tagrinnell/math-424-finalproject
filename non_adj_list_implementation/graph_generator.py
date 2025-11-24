#
# Script to generate a graph and write to out to a file.
#
import argparse

def init_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("output_filename")
    parser.add_argument("num_vertices")
    return parser

def main():
    parser = init_argparser()
    args = parser.parse_args()

    for x in args:
        print(x)

    pass

if __name__ == "__main__":
    main()
