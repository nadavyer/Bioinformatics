import sys
import utils
from datetime import datetime

# constants
k = 10
T = 50
X = 50
R = 25500  # we ignore the outer diagonals of the average lengths of the seqs top 15% as abs(i - j) - keep only the


# 85% in the middle


def find_hsps_of_seq(k, T, sub_matrix, seq1, seq2, alphabet):
    db_dict = utils.build_db(seq2, k)
    return utils.get_hsps(seq1, db_dict, k, sub_matrix, alphabet, T)


def find_hsps_all_seqs(k, T, sub_matrix, seqs_dict, alphabet, pairs_running_time):
    pairs_seqs_hsps = {}
    for i, key1 in enumerate(seqs_dict):
        for j, key2 in enumerate(seqs_dict):
            if i < j:
                startime = datetime.now()
                pairs_seqs_hsps[key1, key2] = find_hsps_of_seq(k, T, sub_matrix,
                                                               seqs_dict[key1], seqs_dict[key2], alphabet)

                pairs_running_time[key1, key2] = datetime.now() - startime
    return pairs_seqs_hsps


def extend_hsps_to_msps(pairs_hsps, sub_matrix, seqs_dict, X, pairs_running_time):
    pairs_msps = {}
    pairs_diagonals = {}
    for seqs_ids, hsps in pairs_hsps.items():
        startime = datetime.now()
        pairs_msps[seqs_ids] = []
        for hsp in hsps:
            if hsp.diagonal() not in pairs_diagonals.keys() or hsp.seq1_start > pairs_diagonals[hsp.diagonal()] \
                    and abs(hsp.diagonal()) < R:
                msp = utils.extend_hsp(seqs_dict[seqs_ids[0]],
                                       seqs_dict[seqs_ids[1]],
                                       hsp, sub_matrix, X)
                pairs_diagonals[hsp.diagonal()] = msp.seq1_end
                if not pairs_msps[seqs_ids].__contains__(msp):
                    pairs_msps[seqs_ids].append(msp)
        pairs_running_time[seqs_ids] += datetime.now() - startime
    return pairs_msps


def main():
    start_time = datetime.now()
    pairs_running_time = {}

    args = sys.argv
    utils.validate_input(args)
    sub_matrix_path, seqs_paths = utils.split_input(args)
    sub_matrix = utils.read_scoring_matrix(sub_matrix_path)
    alphabet = utils.read_alphabet(sub_matrix_path)
    seqs_dict = utils.get_seqs(seqs_paths)
    pairs_hsps = find_hsps_all_seqs(k, T, sub_matrix, seqs_dict, alphabet, pairs_running_time)
    pairs_msps = extend_hsps_to_msps(pairs_hsps, sub_matrix, seqs_dict, X, pairs_running_time)
    utils.print_pairs_msps_count(pairs_msps)
    pairs_graphs = utils.gen_graphs(pairs_msps)
    pairs_scores = utils.runDAG(pairs_graphs)
    utils.gen_output_file(pairs_scores)

    timedelta = datetime.now() - start_time
    print(pairs_running_time)
    print("total running time: ", timedelta)


if __name__ == '__main__':
    main()
