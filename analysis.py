__author__ = 'Suleymanov'

import numpy as np
from itertools import izip
from collections import namedtuple
from Bio import pairwise2
from input_data import hchains, human_set, llama_set
from utils import distance, hamm, adj_distance
from utils import DIST, VLS, COUNT, REGION, GAP

ProcessResult = namedtuple('ProcessResult', 'closest stats')


class CDRAnalyzer(object):
    """ Perform complete analysis for single heavy chain.
        Input:
        - heavy chain sequence of interest with defined CDRs
        - human heavy chains collection
        - llama heavy chains collection
        Task:
        1. Find closest CDRs in human/llama collections
        2. Provide exploratory analysis/statistics on amino acids replacements/indels
        3. Inference on whom heavy chain of interest belongs - human or llama
        Main task:
        - Prove that chain of interest came from human;
        if this is not possible, show the most 'llama' positions in chain. """
    def __init__(self, hchain):
        self.hchain = hchain
        self.result = Result(self.hchain.name)

    def process(self, cdr_set, thresh):
        """ Main processing procedure:
        - find closest CDRs in human/llama collections for 3 distances
        (hamming, levenshtein, adjusted)
        - compute region statistics for various thresholds
        :param thresh: float
        :return: nothing
        """
        print 'Processing chain: ' + self.hchain.name
        print 'CDR set name: ' + cdr_set.name
        closest = {}
        regions = {}
        for i, cdr in enumerate(self.hchain.cdr):
            item_key = 'CDR' + str(i + 1)
            print '\tProcessing ' + item_key
            print '\tSearching closest using Hamming distance...'
            closest['Hamming'][item_key] = self.find_closest(cdr_set, i + 1, hamm)
            print '\tSearching closest using Levenshtein distance...'
            closest['Levenshtein'][item_key] = self.find_closest(cdr_set, i + 1, distance)
            print '\tSearching closest using adjusted distance...'
            closest['Adjusted'][item_key] = self.find_closest(cdr_set, i + 1, adj_distance)
            print '\tCollecting statistics with threshold: %s' % thresh
            regions[item_key] = self.region_stats(cdr_set, i + 1, thresh)
        self.result.add_result(cdr_set.name, ProcessResult(closest=closest, stats=regions))

    def find_closest(self, cdr_set, region_ind, dist_func):
        """ Find closest sequences for given CDR index (1-3)
        according to given distance function.
        :param cdr_set: CDRSet (cdr_set.py)
        :param region_ind: int
        :param dist_func: function
        :return: dict
        """
        result = {}
        cdr = self.hchain.cdr[region_ind - 1]
        for item in cdr_set.from_region(region_ind):
            if len(item) == 0:
                continue
            dist = dist_func(item, cdr)
            if dist:
                if dist in result:
                    result[dist].append(item)
                else:
                    result[dist] = [item]
        min_dist = min(result.keys())
        return {
            REGION: cdr,
            DIST: min_dist,
            COUNT: len(result[min_dist]),
            VLS: np.unique(result[min_dist])
        }

    def region_stats(self, cdr_set, region_ind, thresh):
        """ Collect statistics for all examples from CDR[region_ind]
        according to given threshold value.
        :param cdr_set: CDRSet (cdr_set.py)
        :param region_ind: int
        :param thresh: float
        :return: list of dicts
        """
        assert 0 <= thresh <= 1.0
        assert 1 <= region_ind <= 3
        cdr = self.hchain.cdr[region_ind - 1]
        stats_data = [{} for i in xrange(len(cdr))]
        for item in cdr_set.from_region(region_ind):
            if len(item) == 0:
                continue
            aligns = pairwise2.align.globalxx(cdr, item)[0]
            cdr_aligned = aligns[0]
            item_aligned = aligns[1]
            assert len(cdr_aligned) == len(item_aligned) == aligns[4]
            # dist = adj_distance(cdr, item)
            dist = 1.0 * distance(cdr, item) / aligns[4]  # for faster calculation
            if dist <= thresh:
                i = 0
                for ch1, ch2 in izip(cdr_aligned, item_aligned):
                    if ch1 is not GAP:
                        stats_data[i][ch2] += 1 if ch2 in stats_data[i].keys() else 1
                        i += 1
        return stats_data


class Result(object):
    """ Analysis result keeper. """
    def __init__(self, name):
        self.name = name
        self.entries = {}

    def add_result(self, name, process_result):
        """ Add CDRs collection analysis result.
        :param name: collection name
        :param process_result: ProcessResult
        :return: nothing
        """
        self.entries[name] = process_result

    def __str__(self):
        from json import dumps
        return dumps(self.entries, indent=4)


def entry_point(thresh):
    hchain_analyzers = [CDRAnalyzer(hchain) for hchain in hchains]
    for analyzer in hchain_analyzers:
        analyzer.process(human_set, thresh)
        analyzer.process(llama_set, thresh)
    return hchain_analyzers

if __name__ == '__main__':
    analyzers = entry_point(0.5)
    f_names = ['results1.txt', 'results2.txt']
    for analyzer, f_name in izip(analyzers, f_names):
        with open(f_name, 'w') as f:
            f.write(str(analyzer.result))
