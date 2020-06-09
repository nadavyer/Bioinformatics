
from HSP import HSP


def validateInput(args):
    MIN_INPUT = 4
    ERROR_MESSAGE = "Input cmd should be with format:  [path_to_substitution_matrix] [path_to_seq_file_1] [path_to_seq_file_2] .... [path_to_seq_file_N]" \
                    "\nFor example try run: " \
                    "python pairwise.py sub_mat.txt A.fasta B.fasta C.fasta"
    if len(args) < MIN_INPUT:
        raise ValueError(ERROR_MESSAGE)

def split_input(args):
    sub_mat = args[1]
    input_sq = args[2:]
    return sub_mat, input_sq


def load_matrix(matrix_filename):
    with open(matrix_filename) as matrix_file:
        matrix = matrix_file.read()
    lines = matrix.strip().split('\n')

    header = lines.pop(0)
    columns = header.split()
    matrix = {}

    for row in lines:
        entries = row.split()
        row_name = entries.pop(0)
        matrix[row_name] = {}

        if len(entries) != len(columns):
            raise Exception('Improper entry number in row')
        for column_name in columns:
            matrix[row_name][column_name] = entries.pop(0)

    return matrix

def blossomise_matrix(matrix):
    bloss_matrix = {}
    for line in matrix:
        for val in matrix[line]:
            bloss_matrix[line, val] = matrix[line][val]
    return bloss_matrix

def gen_sub_matrix():


def align(seq1, seq2, scoring_matrix):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be with the same length")
    score = 0
    for i in range(0, len(seq1)):
        score = score + scoring_matrix[seq1[i], seq2[i]]

    return score


def find_neighbors(kmer, scoring_matrix, alphabet, T):
    neighbors = []
    max_score = align(kmer, kmer, scoring_matrix)

    if max_score >= T:
        find_neighbors_rec(kmer, kmer, 0, max_score, alphabet, neighbors, scoring_matrix, T)

    return neighbors


def find_neighbors_rec(kmer, neighbor, pos, curr_score, alphabet, neighbors, scoring_matrix, T):
    if pos == len(kmer):
        neighbors.append(neighbor)
    else:
        curr_score = curr_score - align(kmer[pos], kmer[pos], scoring_matrix)
        for letter in alphabet:
            curr_score = curr_score + align(kmer[pos], letter, scoring_matrix)
            if curr_score >= T: find_neighbors_rec(kmer, neighbor[:pos] + letter, pos + 1, curr_score, alphabet,
                                                   neighbors, scoring_matrix, T)

            curr_score = curr_score - align(kmer[pos], letter, scoring_matrix)


def get_hsps(query, db_dict, k, scoring_matrix, T):
    hsps = []
    kmers = {}

    for query_start in range(0, len(query) - k + 1):
        kmer = query[query_start:query_start + k]
        if kmer not in kmers.keys():
            neighbors = find_neighbors(kmer, scoring_matrix, ALPHABET, T)
            kmers[kmer] = neighbors
            for neighbor in neighbors:
                if neighbor in db_dict.keys():
                    score = align(kmer, neighbor, scoring_matrix)
                    for db_start in db_dict[neighbor]:
                        new_hsp = HSP(query_start, query_start + k - 1, db_start, db_start + k - 1, score)
                        hsps.append(new_hsp)
    return hsps