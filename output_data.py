__author__ = 'Suleymanov'

import os
from cdr_set import AMINO_ACIDS
from input_data import human_set, results_dir


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
