__author__ = 'Suleymanov'

from itertools import izip
from Levenshtein import distance, editops, inverse
from input_data import blosum
from cdr_set import GAP

blosum_alphabet = blosum[0]
blosum_matcher = blosum[1]


def hamming(str1, str2):
    """ Hamming distance for 2 strings of equal size
    :param str1: string
    :param str2: string
    :return: int
    """
    assert len(str1) == len(str2)
    return sum(ch1 != ch2 for ch1, ch2 in izip(str1, str2))


def align(str1, str2):
    """ Align two strings to same length so that same symbols match.
    :param str1: string
    :param str2: string
    :return: tuple of strings
    """
    def _align(s, ed_ops):
        """ Align single string.
        :param s: string
        :param ed_ops: list of edit operations
        :return: string
        """
        result = ''
        inds = [op[1] for op in ed_ops if op[0] == 'insert']
        for i in xrange(len(s)):
            for k in xrange(inds.count(i)):
                result += GAP
            result += s[i]
        return result

    if str2 >= str1:
        e = editops(str1, str2)
        s1 = _align(str1, e)
        s2 = _align(str2, inverse(e))
    else:
        e = editops(str2, str1)
        s1 = _align(str1, inverse(e))
        s2 = _align(str2, e)
    if len(s1) != len(s2):
        print s1
        print s2
        print str1
        print str2
        raise Exception('Invalid strings')
    return s1, s2, 1.0 * distance(str1, str2) / len(s1)


def sim(ch1, ch2):
    index = blosum_alphabet.index(ch2)
    return blosum_matcher[ch1][index]
