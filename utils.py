from HSP import HSP
from copy import copy


def validate_input(args):
    MIN_INPUT = 4
    ERROR_MESSAGE = "Input cmd should be with format:  [path_to_substitution_matrix] [path_to_seq_file_1] [path_to_seq_file_2] .... [path_to_seq_file_N]" \
                    "\nFor example try run: " \
                    "python pairwise.py sub_mat.txt A.fasta B.fasta C.fasta"
    if len(args) < MIN_INPUT:
        raise ValueError(ERROR_MESSAGE)


def split_input(args):
    sub_mat = args[1]
    input_seqs = args[2:]
    return sub_mat, input_seqs


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


def get_seqs(seq_files_names):
    seqs = {}
    for seq_file_name in seq_files_names:
        seq_id, seq = read_seq_file(seq_file_name)
        seqs[seq_id] = seq
    return seqs


def get_alphabet(matrix):
    alphabet = ''
    for letter in matrix:
        if letter not in alphabet:
            alphabet = alphabet + letter
    return alphabet


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
            bloss_matrix[line, val] = int(matrix[line][val])
    return bloss_matrix


def gen_sub_matrix_and_alphabet(matrix_filename):
    sub_mat = load_matrix(matrix_filename)
    alphabet = get_alphabet(sub_mat)
    bloss_matrix = blossomise_matrix(sub_mat)
    return bloss_matrix, alphabet


def build_db(db, k):
    db_dict = {}
    for i in range(0, len(db) - k + 1):
        if db[i: i + k] not in db_dict:
            db_dict[db[i: i + k]] = []  # init an empty list for the k's
        db_dict[db[i: i + k]].append(i)
    return db_dict


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


def get_hsps(query, db_dict, k, scoring_matrix, alphabet, T):
    hsps = []
    kmers = {}

    for seq1_start in range(0, len(query) - k + 1):
        kmer = query[seq1_start:seq1_start + k]
        if kmer not in kmers.keys():
            neighbors = find_neighbors(kmer, scoring_matrix, alphabet, T)
            kmers[kmer] = neighbors
            for neighbor in neighbors:
                if neighbor in db_dict.keys():
                    score = align(kmer, neighbor, scoring_matrix)
                    for seq2_start in db_dict[neighbor]:
                        new_hsp = HSP(seq1_start, seq1_start + k - 1, seq2_start, seq2_start + k - 1, score)
                        hsps.append(new_hsp)
    return hsps


def extend_left(query, db, hsp, scoring_matrix, X):
    """returns the left extension with the maximal score"""
    msp = copy(hsp)
    i = 0
    max_score = msp.score
    cur_score = msp.score
    query_mark = msp.seq1_start
    db_mark = msp.seq2_start

    while query_mark > 0 and db_mark > 0 and max_score - X <= cur_score:
        query_mark -= 1
        db_mark -= 1
        cur_score += scoring_matrix[query[query_mark], db[db_mark]]
        if max_score <= cur_score:  # save new max
            max_score = cur_score
            msp.seq1_start = query_mark
            msp.seq2_start = db_mark
            msp.score = cur_score
    return msp


def extend_right(query, db, hsp, scoring_matrix, X):
    """returns the right extension with the maximal score"""

    msp = copy(hsp)
    i = 0
    max_score = msp.score
    cur_score = msp.score
    query_mark = msp.seq1_end
    db_mark = msp.seq2_end

    while query_mark < len(query) - 1 and db_mark < len(db) - 1 and max_score - X <= cur_score:
        query_mark += 1
        db_mark += 1
        cur_score += scoring_matrix[query[query_mark], db[db_mark]]
        if max_score <= cur_score:  # save new max
            max_score = cur_score
            msp.seq1_end = query_mark
            msp.seq2_end = db_mark
            msp.score = cur_score
    return msp


def extend_hsp(query, db, hsp, scoring_matrix, X):
    msp_left = extend_left(query, db, hsp, scoring_matrix, X)
    msp = extend_right(query, db, msp_left, scoring_matrix, X)

    return msp


def get_top_msps(msps, top):
    sorted_msps = sorted(list(msps), reverse=True, key=lambda msp: msp.score)
    return sorted_msps[:top]


def get_top_score(pairs_msps):
    top_scores = {}
    for pair in pairs_msps:
        top_scores[pair] = get_top_msps(pairs_msps[pair], 1)
    return top_scores


def gen_output_file(pairs_msps):
    output_file = open("scores.txt", "w")
    for pair, top_msp in pairs_msps.items():
        output_line = f"{pair[0]}\t{pair[1]}\t{top_msp[0].score}\n"
        output_file.write(output_line)
    output_file.close()
