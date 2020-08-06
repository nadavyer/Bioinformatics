import os
import sys
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
import matplotlib.pyplot as plt


def read_files(dir_path):
    seq_files = os.listdir(dir_path)
    seq_len_dict = {}
    for file in seq_files:
        seq_id, seq = read_seq_file(dir_path + file)
        seq_len_dict[seq_id] = seq

    return seq_len_dict


def read_seq_file(seq_file):
    seq_id = ''
    seq = ''
    with open(seq_file) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip()[1:]
            else:
                seq += line.strip()

        return seq_id, seq


# gen lower diagonal matrix from the output file generated from running pairwise.py
def gen_matrix_from_pair_ids_and_value(path):
    last_left_col_name = ''
    seqs_ids = []
    row = []
    matrix = []
    not_first = False
    need_to_add_name = True

    with open(path) as f:
        for line in f:
            splited = line.strip('\n').split('\t')
            name = splited[0]
            if not name == last_left_col_name:
                last_left_col_name = name
                if need_to_add_name and not seqs_ids.__contains__(splited[1]):
                    row.append(0)
                    seqs_ids.append(name)
                if not_first:
                    row.reverse()
                    matrix.append(row)
                    row = [0]
                    need_to_add_name = False
                not_first = True
            if need_to_add_name and not seqs_ids.__contains__(splited[1]):
                seqs_ids.append(splited[1])
            if not row.__contains__(int(splited[2])):
                row.append(int(splited[2]))
        row.reverse()
        matrix.append(row)
        matrix.append([0])
        matrix.reverse()
        seqs_ids.reverse()

        return matrix, seqs_ids

    # print matrix


def update_distance_file(id1, id2, distance):
    output_file = open("distances.txt", "a+")
    output_line = f"{id1}\t{id2}\t{distance}\n"
    output_file.write(output_line)
    output_file.close()


def gen_score_file_to_distance_file(scores_path, seqs):
    with open(scores_path) as f:
        for line in f:
            line_args = line.strip('\n').split('\t')
            update_distance_file(line_args[0], line_args[1],
                                 min(int(len(seqs[line_args[0]])), int(len(seqs[line_args[1]]))) * 5 - int(
                                     line_args[2]))


"""args: dir_of_seqs_path/"""


def main():
    seqs = read_files(sys.argv[1])
    gen_score_file_to_distance_file("scores.txt", seqs)
    matrix, seq_ids = gen_matrix_from_pair_ids_and_value(path='distances.txt')

    print(matrix)
    print(len(matrix))
    print(len(seq_ids))
    dm = DistanceMatrix(names=seq_ids, matrix=matrix)

    print(dm)

    constructor = DistanceTreeConstructor()

    tree = constructor.nj(dm)

    fig = plt.figure(figsize=(12, 5), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes)


if __name__ == '__main__':
    main()
