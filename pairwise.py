import sys
import utils
from datetime import datetime

# constants
k = 10
T = 50
X = 50
  # we ignore the outer diagonals of the average lengths of the seqs top 15% as abs(i - j) - keep only the
           # 85% in the middle


def main():
    args = sys.argv
    utils.validate_input(args)
    sub_matrix_path, seqs_paths = utils.split_input(args)
    sub_matrix = utils.read_scoring_matrix(sub_matrix_path)
    alphabet = utils.read_alphabet(sub_matrix_path)
    seqs_dict = utils.get_seqs(seqs_paths)
    i = 0
    for id1, seq1 in seqs_dict.items():
        j = 0
        for id2, seq2 in seqs_dict.items():
            if i < j:
                timestart = datetime.now()
                hsps = utils.find_hsps(k, T, sub_matrix, seq1, seq2, alphabet)
                msps = utils.extend_hsps_to_msps(hsps, sub_matrix, seq1, seq2, X)
                graph = utils.gen_graph(msps)
                score = utils.runDAG(graph)
                timedelta = datetime.now() - timestart
                utils.update_output_file(id1, id2, score)
                print(f"({id1}, {id2}):\t{timedelta}\t msps: {len(msps)}")
            j += 1
        i += 1


if __name__ == '__main__':
    main()
