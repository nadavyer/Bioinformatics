class HSP:
    def __init__(self, seq1_start, seq1_end, seq2_start, seq2_end, score):
        self.seq1_start = seq1_start
        self.seq1_end = seq1_end
        self.seq2_start = seq2_start
        self.seq2_end = seq2_end
        self.score = score

    def size(self):
        return self.seq1_end - self.seq1_start

    def __eq__(self, other):
        """equals function"""
        if isinstance(other, self.__class__):
            is_seq1_eq = self.seq1_start == other.seq1_start and self.seq1_end == other.seq1_end
            is_seq2_eq = self.seq2_start == other.seq2_start and self.seq2_end == other.seq2_end
            return is_seq1_eq and is_seq2_eq and self.score == other.score
        return False

    def __ne__(self, other):
        """Define a non-equality test"""
        return not self.__eq__(other)

    def __hash__(self):
        return hash(str((self.seq1_start, self.seq1_end, self.seq2_start, self.seq2_end, self.score)))

    def __str__(self):
        return (
                f'Score: {self.score}\n'
                f'Seq1: [{self.seq1_start}, {self.seq1_end}]\n'
                f'Seq2: [{self.seq2_start}, {self.seq2_end}]\n'
                )

    def __repr__(self):
        return self.__str__()