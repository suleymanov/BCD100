__author__ = 'Suleymanov'

import os
import numpy as np
from itertools import izip
from Bio import pairwise2
from input_data import hchains, hchains_names, human_set, llama_set, results_dir
from utils import distance, hamming, GAP, AMINO_ACIDS


class CDRAnalyzer(object):
    def __init__(self, hchain):
        self.hchain = hchain

    def find_closest(self, cdr_set, distance_func=hamming):
        """ Find k closest CDRs for all 3 regions in cdr_set
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


def record_closest(data, f_name, header=None):
    """ Print results of CDRAnalyzer.find_closest()
    :param data: dict
    :param f_name: string
    :param header: string
    :return: nothing
    """
    data_lines = []
    num_lines = max([len(data[key]['values']) for key in data])
    for i, line in enumerate(xrange(num_lines)):
        row_vls = []
        for key in data:
            if i < len(data[key]['values']):
                row_vls.append(data[key]['values'][i])
            else:
                row_vls.append('---')
        data_lines.append('\t'.join(row_vls))

    with open(f_name, 'a') as f:
        if header:
            f.write(header + '\n')
        f.write('\t' + '\t'.join([key for key in data]) + '\n')
        f.write('Region\t' + '\t'.join([data[key]['region'] for key in data]) + '\n')
        f.write('Distance\t' + '\t'.join([str(data[key]['dist']) for key in data]) + '\n')
        f.write('Count\t' + '\t'.join([str(data[key]['count']) for key in data]) + '\n')
        for line in data_lines:
            f.write(line + '\n')
        f.write('\n')


def record_diff(data, f_name, header=None):
    """ Print unique amino acids dependent on position and CDR region of two CDR datasets.
    :param data: dict
    :param f_name: string
    :param header: string
    :return: nothing
    """
    with open(f_name, 'a') as f:
        if header:
            f.write(header + '\n')
        for key in data:  # traverse length values
            f.write(str(key) + '-len words\n')
            for region_ind, region_data in enumerate(data[key]):  # traverse region data
                if region_data:  # omit empty data
                    f.write('CDR' + str(region_ind + 1) + '\n')
                    for i, vls in enumerate(region_data):
                        if vls:  # omit empty data
                            str_vls = [str(val) for val in vls]
                            f.write(str(i) + '\t' + ' '.join(str_vls) + '\n')
                f.write('\n')
            f.write('\n')


def record_stats(cdr, stats_data, f_name, header=None):
    alphabet = sorted([item for sub in [AMINO_ACIDS, ['-']] for item in sub])
    with open(f_name, 'a') as f:
        if header:
            f.write(header + '\n')
        f.write('\t' + '\t'.join(list(cdr)) + '\n')
        for letter in alphabet:
            row = []
            for stat_item in stats_data:
                if letter in stat_item:
                    row.append(stat_item[letter])
                else:
                    row.append(0)
            if any([val != 0 for val in row]):
                f.write(letter + '\t' + '\t'.join([str(val) for val in row]) + '\n')
        f.write('\n')


def closest():
    """ Find and record closest CDRs to given ones.
    :return: nothing
    """
    analyzers = [CDRAnalyzer(hchain) for hchain in hchains]
    closest_human = [analyzer.find_closest(human_set, distance) for analyzer in analyzers]
    closest_llama = [analyzer.find_closest(llama_set, distance) for analyzer in analyzers]
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

    f_name = os.sep.join([results_dir, 'human-stats.txt'])
    if os.path.exists(f_name):
        os.remove(f_name)
    for stats_data, hchain in izip(stats_human, hchains):
        record_stats(hchain.cdr[region_ind - 1], stats_data, f_name, hchain.name)

    f_name = os.sep.join([results_dir, 'llama-stats.txt'])
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


def record_lengthers():
    """ Record CDR-lengthers of interest (lengths 5, 17, 20).
    :return: nothing
    """
    def record_lengthers(data, f_name):
        with open(f_name, 'w') as f:
            for record in data:
                f.write(record + '\n')

    def create_file(short_name):
        name = os.sep.join([results_dir, short_name])
        if os.path.exists(name):
            os.remove(name)
        return name

    lengths = [5, 17, 20]
    for i, item in enumerate(lengths):
        f_name = create_file('human %s-lenghters (CDR%s).txt' % (item, i + 1))
        subset = human_set.all_k_lengthers(item, i + 1)
        record_lengthers(subset, f_name)
        # record_lengthers(human_set.all_k_lengthers(item, i + 1), f_name)


if __name__ == '__main__':
    stats(0.4, 3)
    # closest()
    # record_lengthers()
    # diff()
