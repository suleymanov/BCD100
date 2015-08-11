import numpy as np
from utils import sim  # similarity function/matrix
from cdr_set import GAP

d = -5


def create_f(a, b):
    f = np.zeros((len(a), len(b)))
    for j in xrange(len(b)):
        f[0, j] = d * j
    for i in xrange(len(a)):
        f[i, 0] = d * i
    for i in xrange(1, len(a)):
        for j in xrange(1, len(b)):
            match = f[i-1, j-1] + sim(a[i], b[j])
            delete = f[i-1, j] + d
            insert = f[i, j-1] + d
            f[i, j] = max(match, delete, insert)
    return f


def align(a, b):
    f = create_f(a, b)
    al_a = ''
    al_b = ''
    i = len(a) - 1
    j = len(b) - 1
    while i > 0 or j > 0:
        sim_val = sim(a[i], b[j])
        if i > 0 and j > 0 and f[i, j] == f[i-1, j-1] + sim_val:
            al_a = a[i] + al_a
            al_b = b[j] + al_b
            i -= 1
            j -= 1
        elif i > 0 and f[i, j] == f[i-1, j] + d:
            al_a = a[i] + al_a
            al_b = GAP + al_b
            i -= 1
        else:
            al_a = GAP + al_a
            al_b = b[j] + al_b
            j -= 1
    return al_a, al_b


if __name__ == '__main__':
    a_str = 'GCATGCT'
    b_str = 'GATTACA'
    res = align(a_str, b_str)

    print res[0]
    print res[1]
