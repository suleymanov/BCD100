__author__ = 'Suleymanov'

# Kabat structure:
# Sequence id <-> FR/CDR indices
# FR indices: 0-1, 4-5, 8-9, 12-13, CDR indices: 2-3, 6-7, 10-11
# FASTA structure:
# Sequence id <-> sequence string

import string
import numpy as np
from Bio import SeqIO
from utils import AMINO_ACIDS

# X - special symbol for missing or unimportant data
# (https://en.wikipedia.org/wiki/Proteinogenic_amino_acid)


class CDRSet(object):
    def __init__(self, fasta_file, kabat_file, name):
        """ Initialize CDR set from sequences and regions files.
        :param fasta_file: sequences file, string
        :param kabat_file: regions file, string
        :param name: string (collection name)
        :return: nothing
        """
        self.name = name
        self.fasta_data = self.read_fasta(fasta_file)
        self.kabat_data = self.read_kabat(kabat_file)
        self.eliminate()
        self._validate()
        self.cdr_all = self.all_cdr()
        self.lengths = np.unique([len(item) for cdr in self.cdr_all for item in cdr])

    def eliminate(self):
        """ Remove all .kabat arrays that have no matching sequences.
        :return: nothing
        """
        all_ids = [item.id for item in self.fasta_data]
        for k in self.kabat_data.keys():
            if k not in all_ids:
                self.kabat_data.pop(k)

    def _validate(self):
        for values in self.kabat_data.itervalues():
            assert len(values) == 14

        invalid_letters = list(set(list(string.uppercase[:26])).difference(list(AMINO_ACIDS)))
        invalid_found = []
        for record in self.fasta_data:
            for letter in invalid_letters:
                if str(record.seq).count(letter) > 0:
                    invalid_found.append(letter)
        if invalid_found:
            raise Exception('Invalid symbols found in data: ' + ', '.join(np.unique(invalid_found)))

    def all_cdr(self):
        """ Collect all CDR sequences.
        :return: list
        """
        cdr_data = []
        for item in self.fasta_data:
            inds = [ind - 1 for ind in self.kabat_data[item.id]]
            cdr_slices = [slice(inds[2], inds[4]), slice(inds[6], inds[8]), slice(inds[10], inds[12])]
            cdr_data.append([str(item.seq)[sl] for sl in cdr_slices])
        return cdr_data

    def all_k_lengthers(self, k, region_ind=None):
        """ Return all CDR sequences of length k from region_ind (1, 2, 3 - if specified).
        (Zero-length is allowed).
        All k-lengthers if no region_ind
        :param k: int
        :param region_ind: int
        :return: list
        """
        assert 0 <= k <= max(self.lengths)
        if region_ind:
            assert 1 <= region_ind <= 3
            return [cdr[region_ind - 1] for cdr in self.cdr_all if len(cdr[region_ind - 1]) == k]
        return [item for cdr in self.cdr_all for item in cdr if len(item) == k]

    def from_region(self, region_ind):
        """ Return all CDR sequences from region_ind (1, 2, 3).
        :param region_ind: int
        :return: list
        """
        return [cdr[region_ind - 1] for cdr in self.cdr_all]

    def diff(self, other):
        """ Find all characters that don't occur in 'other' at same positions.
        :param other: CDRSet
        :return: dict; keys are lengths,
        values are lists with lists of unique amino acids for each region;
        empty if no unique letters
        """
        assert isinstance(other, CDRSet)
        diff_data = {}
        for length in self.lengths:
            if length not in other.lengths or length == 0:
                    continue
            diff_data[length] = []
            for region_ind in xrange(1, 4):
                lengthers = self.all_k_lengthers(length, region_ind)
                other_lengthers = other.all_k_lengthers(length, region_ind)
                diff = [
                    list(set([s[i] for s in lengthers]).difference(
                        set([s[i] for s in other_lengthers])))
                    for i in xrange(length)
                ]
                if not any(diff):
                    diff = []
                diff_data[length].append(diff)
        return diff_data

    def __str__(self):
        s = 'CDR set:' + '\n'
        s += 'Total ' + str(len(self.cdr_all)) + ' sets of regions.\n\n'
        for k in self.lengths:
            vls = [len(self.all_k_lengthers(k, pos)) for pos in xrange(1, 4)]
            s += 'Number of regions of length ' + str(k) + ': ' + str(vls) + '\n'
        return s

    def print_examples(self, num, fname):
        """ Print no more than 'num' examples from all k-lengthers.
        :param num: int > 1
        :param fname: string
        :return: nothing
        """
        with open(fname, 'w') as f:
            for k in self.lengths:
                f.write(str(k) + '-lengthers:\n')
                vls = self.all_k_lengthers(k)
                for i, val in enumerate(vls):
                    if i > num:
                        break
                    f.write(val + '\n')
                f.write('\n')

    @staticmethod
    def read_kabat(file_name):
        """ Read .kabat/.marking file
        :param file_name: string
        :return: dict
        """
        with open(file_name, 'r') as f:
            s = f.read()
        contents = s.split('\n')[:-1]
        data = {}
        for entry in contents:
            strs = entry.split('\t')
            data[strs[0]] = [int(val) for val in strs[1:]]
        return data

    @staticmethod
    def read_fasta(file_name):
        """ Read .fasta file
        :param file_name: string
        :return: list
        """
        handle = open(file_name, 'rU')
        return list(SeqIO.parse(handle, 'fasta'))


def read_all_data():
    pass

if __name__ == '__main__':
    read_all_data()
