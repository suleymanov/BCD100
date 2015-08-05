__author__ = 'Suleymanov'

from itertools import izip


def hamming(str1, str2):
    """ Hamming distance for 2 strings of equal size
    :param str1:
    :param str2:
    :return:
    """
    assert len(str1) == len(str2)
    return sum(ch1 != ch2 for ch1, ch2 in izip(str1, str2))
