__author__ = 'Suleymanov'

import os
import operator
from input_data import hchains, hchains_names, human_set, llama_set, results_dir
from itertools import izip
from utils import hamming


class CDRAnalyzer(object):
    def __init__(self, hchain):
        self.hchain = hchain

    def find_closest(self, cdr_set, k):
        """ Find k closest CDRs for all 3 regions in cdr_set
        :param cdr_set:
        :param k: number of considered closest distances
        :return: dict
        """
        assert k > 0
        result = {}
        for i, cdr in enumerate(self.hchain.cdr):
            k_lengthers = cdr_set.all_k_lengthers(len(cdr), i + 1)
            data = {k: hamming(k, cdr) for k in k_lengthers}
            sorted_data = sorted(data.items(), key=operator.itemgetter(1))
            # result[cdr] = [sorted_data[i] for i in xrange(k)]
            # consider rewriting/optimizing stuff below
            values = [sorted_data[0]]
            count = 1
            for j, item in enumerate(sorted_data[1:]):
                prev = sorted_data[j - 1]
                if prev[1] < item[1]:
                    count += 1
                if count > k:
                    break
                values.append(item)
            item_key = 'CDR' + str(i + 1)
            result[item_key] = {'region': cdr, 'values': values}
        return result


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
                val = data[key]['values'][i]
                row_vls.append(val[0] + ' (' + str(val[1]) + ')')
            else:
                row_vls.append('---')
        data_lines.append('\t'.join(row_vls))

    with open(f_name, 'a') as f:
        if header:
            f.write(header + '\n')
        f.write('\t'.join([key for key in data]) + '\n')
        f.write('\t'.join([data[key]['region'] for key in data]) + '\n')
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


def entry_point():
    # analyzers = [CDRAnalyzer(hchain) for hchain in hchains]
    # num_closest = 1
    # closest_human = [analyzer.find_closest(human_set, num_closest) for analyzer in analyzers]
    # closest_llama = [analyzer.find_closest(llama_set, num_closest) for analyzer in analyzers]
    # f_name_closest = [os.sep.join([results_dir, name + ' (closest).txt'])
    #                   for name in hchains_names]
    # for f_name, human_res, llama_res in izip(f_name_closest, closest_human, closest_llama):
    #     if os.path.exists(f_name):
    #         os.remove(f_name)
    #     record_closest(human_res, f_name, 'human data')
    #     record_closest(llama_res, f_name, 'llama data')

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
    entry_point()
