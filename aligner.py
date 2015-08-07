import numpy as np
from utils import sim  # similarity function/matrix
from cdr_set import GAP

d = -5


def create_F(a, b):
    F = np.zeros((len(a), len(b)))
    for j in xrange(len(b)):
        F[0, j] = d * j
    for i in xrange(len(a)):
        F[i, 0] = d * j
    for i in xrange(1, len(a)):
        for j in xrange(1, len(b)):
            match = F[i-1, j-1] + sim(a[i], b[j])
            delete = F[i-1, j] + d
            insert = F[i, j-1] + d
            F[i, j] = max(match, delete, insert)
    return F


def align(a, b):
    f = create_F(a, b)
    al_a = ''
    al_b = ''
    i = len(a) - 1
    j = len(b) - 1
    while i >= 0 and j >= 0:
        score = f[i, j]
        score_diag = f[i-1, j-1]
        score_up = f[i, j-1]
        score_left = f[i-1, j]
        if score == score_diag + sim(a[i], b[j]):
            al_a = a[i] + al_a
            al_b = b[j] + al_b
            i -= 1
            j -= 1
        elif score == score_left + d:
            al_a = a[i] + al_a
            al_b = GAP + al_b
            i -= 1
        else:
            al_a = GAP + al_a
            al_b = b[j] + al_b
            j -= 1
    while i >= 0:
        al_a = a[i] + al_a
        al_b = GAP + al_b
        i -= 1
    while j >= 0:
        al_a = GAP + al_a
        al_b = b[j] + al_a
        j -= 1
    return al_a, al_b
