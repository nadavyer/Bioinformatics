import sys
import utils

from HSP import HSP


def main():
    args = sys.argv
    # utils.validateInput(args)
    sub_matrix_file_name, seq_names = utils.split_input(args)
    sub_matrix = utils.load_matrix(sub_matrix_file_name)
    blossomise_matrix = utils.blossomise_matrix(sub_matrix)


if __name__ == '__main__':
    main()
