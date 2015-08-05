__author__ = 'Suleymanov'

from Bio import SeqIO


class HeavyChain(object):
    def __init__(self, seq_name, seqs_file, cdr_file):
        """ Initialize heavy chain sequence & CDR
        :param seq_name: string
        :param seqs_file: string
        :param cdr_file: string
        :return:
        """
        self.name = seq_name
        self.seq = self._read(self.name, seqs_file)
        self.cdr = self._read_cdr(self.name, cdr_file)

    def __str__(self):
        s = 'Heavy chain sequence name: ' + self.name + '\n'
        for i, cdr in enumerate(self.cdr):
            s += 'CDR #' + str(i + 1) + ': ' + cdr + '\n'
        return s

    @staticmethod
    def _read(name, seqs_file):
        """ Read sequence from file.
        :param name: string
        :param seqs_file: string
        :return: Bio.SeqRecord.SeqRecord
        """
        handle = open(seqs_file, 'rU')
        seqs = list(SeqIO.parse(handle, 'fasta'))
        assert name in [seq.id for seq in seqs]
        return next(seq for seq in seqs if seq.id == name)

    @staticmethod
    def _read_cdr(seq_name, cdr_file):
        """ Read heavy chain complementarity determining regions (CDR) from file.
        :param seq_name: string
        :param cdr_file: string
        :return: list of string
        """
        handle = open(cdr_file, 'rU')
        seqs = list(SeqIO.parse(handle, 'fasta'))
        assert seq_name in [seq.id for seq in seqs]
        return next(str(seq.seq).split('$') for seq in seqs if seq.id == seq_name)


def test():
    pass


if __name__ == '__main__':
    test()
