__author__ = 'Suleymanov'

from Levenshtein import distance, hamming
from Bio.SubsMat.MatrixInfo import blosum62
from Bio import pairwise2

REGION = 'region'
DIST = 'distance'
COUNT = 'count'
VLS = 'values'
AMINO_ACIDS = ['D', 'T', 'S', 'E', 'P', 'G', 'A', 'C', 'V', 'M',
               'I', 'L', 'Y', 'F', 'H', 'K', 'R', 'W', 'Q', 'N', 'X']
GAP = '-'


def hamm(str1, str2):
    """ Levenshtein distance with option to throw None
    for strings of different lengths.
    :param str1: string
    :param str2: string
    :return: int or None
    """
    return hamming(str1, str2) if len(str1) is len(str2) else None


def adj_distance(str1, str2):
    """ Adjusted distance between two strings.
    :param str1: string
    :param str2: string
    :return: float
    """
    aligns = pairwise2.align.globalxx(str1, str2)[0]
    al_len = aligns[4]
    return 1.0 * distance(str1, str2) / al_len


def sim(ch1, ch2):
    """ Similarity between two amino acids.
    :param ch1: string
    :param ch2: string
    :return: int
    """
    if (ch1, ch2) in blosum62.keys():
        return blosum62[(ch1, ch2)]
    else:
        return blosum62[(ch2, ch1)]
