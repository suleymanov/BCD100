__author__ = 'Suleymanov'

import os
import numpy as np
from itertools import izip
from Bio import pairwise2
from input_data import hchains, hchains_names, human_set, llama_set, results_dir
from output_data import record_stats, record_diff, record_closest, record_lengthers
from utils import distance, hamming, adj_distance
from cdr_set import GAP


class CDRAnalyzer(object):
    def __init__(self, hchain):
        self.hchain = hchain

    def find_closest(self, cdr_set, distance_func=hamming):
        """ Find closest CDRs for all 3 regions in cdr_set
        :param cdr_set:
        :param distance_func: strings similarity function
        :return: dict
        """
        result = {}
        for i, cdr in enumerate(self.hchain.cdr):
            k_lengthers = cdr_set.all_k_lengthers(len(cdr), i + 1)
            data = {}
            for k in k_lengthers:
                dist = distance_func(k, cdr)
                if dist in data:
                    data[dist].append(k)
                else:
                    data[dist] = [k]
            item_key = 'CDR' + str(i + 1)
            min_dist = min(data.keys())
            result[item_key] = {
                'region': cdr,
                'dist': min_dist,
                'count': len(data[min_dist]),
                'values': np.unique(data[min_dist])}
        return result

    def find_closest_adj(self, cdr_set):
        """ Find closest CDRs for all 3 regions in cdr_set
        using adjusted distance function
        :param cdr_set:
        :return: dict
        """
        result = {}
        for i, cdr in enumerate(self.hchain.cdr):
            data = {}
            for length in cdr_set.lengths:
                lengthers = cdr_set.all_k_lengthers(length, i + 1)
                for item in lengthers:
                    dist = adj_distance(item, cdr)
                    if dist in data:
                        data[dist].append(item)
                    else:
                        data[dist] = [item]
            item_key = 'CDR' + str(i + 1)
            min_dist = min(data.keys())
            result[item_key] = {
                'region': cdr,
                'dist': min_dist,
                'count': len(data[min_dist]),
                'values': np.unique(data[min_dist])
            }
        return result

    def region_stats(self, cdr_set, thresh, region_ind):
        """ Calculate statistics for all examples from CDR[region_ind]
        with adjusted distance passing threshold
        :param thresh: float
        :param region_ind: int
        :return: list; an element is statistics for corresponding position
        """
        assert 0 < thresh <= 1.0
        assert 1 <= region_ind <= 3
        cdr = self.hchain.cdr[region_ind - 1]
        stats_data = [{} for i in xrange(len(cdr))]
        for length in cdr_set.lengths:  # try all lengths values in queried CDR region
            lengthers = cdr_set.all_k_lengthers(length, region_ind)
            for item in lengthers:
                if len(item) == 0:
                    continue
                aligns = pairwise2.align.globalxx(cdr, item)[0]
                cdr_aligned = aligns[0]
                item_aligned = aligns[1]
                assert len(cdr_aligned) == len(item_aligned) == aligns[4]
                dist = 1.0 * distance(cdr, item) / aligns[4]
                if dist <= thresh:
                    i = 0
                    for ch1, ch2 in izip(cdr_aligned, item_aligned):
                        if ch1 is not GAP:
                            if ch2 in stats_data[i].keys():
                                stats_data[i][ch2] += 1
                            else:
                                stats_data[i][ch2] = 1
                            i += 1
        return stats_data


def closest():
    """ Find and record closest CDRs to given ones.
    :return: nothing
    """
    analyzers = [CDRAnalyzer(hchain) for hchain in hchains]
    # closest_human = [analyzer.find_closest(human_set, distance) for analyzer in analyzers]
    # closest_llama = [analyzer.find_closest(llama_set, distance) for analyzer in analyzers]
    closest_human = [analyzer.find_closest_adj(human_set) for analyzer in analyzers]
    closest_llama = [analyzer.find_closest_adj(llama_set) for analyzer in analyzers]
    f_name_closest = [os.sep.join([results_dir, name + ' (closest, Levenshtein).txt'])
                      for name in hchains_names]
    for f_name, human_res, llama_res in izip(f_name_closest, closest_human, closest_llama):
        if os.path.exists(f_name):
            os.remove(f_name)
        record_closest(human_res, f_name, 'human data')
        record_closest(llama_res, f_name, 'llama data')


def stats(thresh, region_ind):
    analyzers = [CDRAnalyzer(hchain) for hchain in hchains]
    stats_human = [analyzer.region_stats(human_set, thresh, region_ind) for analyzer in analyzers]
    stats_llama = [analyzer.region_stats(llama_set, thresh, region_ind) for analyzer in analyzers]

    f_name = 'human-stats (thresh ' + str(thresh) + ', CDR' + str(region_ind) + ').txt'
    f_name = os.sep.join([results_dir, f_name])
    if os.path.exists(f_name):
        os.remove(f_name)
    for stats_data, hchain in izip(stats_human, hchains):
        record_stats(hchain.cdr[region_ind - 1], stats_data, f_name, hchain.name)

    f_name = 'llama-stats (thresh ' + str(thresh) + ', CDR' + str(region_ind) + ').txt'
    f_name = os.sep.join([results_dir, f_name])
    if os.path.exists(f_name):
        os.remove(f_name)
    for stats_data, hchain in izip(stats_llama, hchains):
        record_stats(hchain.cdr[region_ind - 1], stats_data, f_name, hchain.name)


def diff():
    """ Find and record symbols present in human, but not llama and vice versa.
    :return: nothing
    """
    human_diff = human_set.diff(llama_set)
    f_name = os.sep.join([results_dir, 'human-diff.txt'])
    if os.path.exists(f_name):
        os.remove(f_name)
    record_diff(human_diff, f_name)

    llama_diff = llama_set.diff(human_set)
    f_name = os.sep.join([results_dir, 'llama-diff.txt'])
    if os.path.exists(f_name):
        os.remove(f_name)
    record_diff(llama_diff, f_name)


if __name__ == '__main__':
    # stats(0.9, 3)
    closest()
    # record_lengthers()
    # diff()
