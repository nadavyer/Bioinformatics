import sys
import utils

from HSP import HSP

# constants
k = 4
T = 23
X = 8


def find_hsps_of_seq(k, T, sub_matrix, seq1, seq2, alphabet):
    db_dict = utils.build_db(seq2, k)
    return utils.get_hsps(seq1, db_dict, k, sub_matrix, alphabet, T)


def find_hsps_all_seqs(k, T, sub_matrix, seqs_dict, alphabet):
    pairs_seqs_hsps = {}
    for i, key1 in enumerate(seqs_dict):
        for j, key2 in enumerate(seqs_dict):
            if i >= j:
                continue
            pairs_seqs_hsps[key1, key2] = find_hsps_of_seq(k, T, sub_matrix,
                                                           seqs_dict[key1], seqs_dict[key2], alphabet)
    return pairs_seqs_hsps


def extend_hsps_to_msps(pairs_hsps, sub_matrix, seqs_dict, X):
    pairs_msps = {}
    for seqs_ids, hsps in pairs_hsps.items():
        pairs_msps[seqs_ids] = set()
        for hsp in hsps:
            msp = utils.extend_hsp(seqs_dict[seqs_ids[0]],
                                   seqs_dict[seqs_ids[1]],
                                   hsp, sub_matrix, X)
            pairs_msps[seqs_ids].add(msp)
    return pairs_msps


def main():
    args = sys.argv
    utils.validate_input(args)
    sub_matrix_file_name, seqs_names = utils.split_input(args)
    sub_matrix, alphabet = utils.gen_sub_matrix_and_alphabet(sub_matrix_file_name)
    seqs_dict = utils.get_seqs(seqs_names)
    pairs_hsps = find_hsps_all_seqs(k, T, sub_matrix, seqs_dict, alphabet)
    pairs_msps = extend_hsps_to_msps(pairs_hsps, sub_matrix, seqs_dict, X)
    pairs_top_msp = utils.get_top_score(pairs_msps)
    utils.gen_output_file(pairs_top_msp)


if __name__ == '__main__':
    main()
