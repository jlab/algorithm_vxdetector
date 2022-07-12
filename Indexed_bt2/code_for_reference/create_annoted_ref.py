#!/usr/bin/python

# This programm takes the variable region reference position
# from variable_region_calc.py and searches for them in the
# aligned gg otus

from os import path as p


def write(region, seq_chrom, annoted_ref):
    with open(annoted_ref, 'a') as t:
        if region['V1_start'] != '':
            for counter, key in enumerate(region):
                if (counter % 2) == 0:
                    wordstart = f'{seq_chrom}\t{str(region[key])}'
                else:
                    var_reg = key[0:2]
                    wordend = f'\t{str(region[key])}\t{var_reg}\n'
                    words = wordstart + wordend
                    t.write(words)


def index(line, region):
    boundary = {'V1_start': 189, 'V1_end': 471,
                'V2_start': 485, 'V2_end': 1867,
                'V3_start': 1915, 'V3_end': 2231,
                'V4_start': 2262, 'V4_end': 4050,
                'V5_start': 4089, 'V5_end': 4521,
                'V6_start': 4652, 'V6_end': 4931,
                'V7_start': 5043, 'V7_end': 5806,
                'V8_start': 5909, 'V8_end': 6426,
                'V9_start': 6449, 'V9_end': 6790}
    for key in boundary:
        list1 = list(line)
        list1[boundary[key]] = '1'
        list1 = ''.join(str(e) for e in list1)
        list1 = list1.replace('-', '')
        region[key] = int(list1.index('1')) + 1
    return region


def main():
    file_ = f'{p.dirname(p.dirname(__file__))}/85_otus_aligned.fasta'
    annoted_ref = f'{p.dirname(p.dirname(__file__))}/annoted_ref.bed'
    with open(file_, 'r') as f:
        for line in f:
            region = {'V1_start': '', 'V1_end': '',
                      'V2_start': '', 'V2_end': '',
                      'V3_start': '', 'V3_end': '',
                      'V4_start': '', 'V4_end': '',
                      'V5_start': '', 'V5_end': '',
                      'V6_start': '', 'V6_end': '',
                      'V7_start': '', 'V7_end': '',
                      'V8_start': '', 'V8_end': '',
                      'V9_start': '', 'V9_end': ''}
            if line.startswith('>'):
                seq_chrom = line[1:]
                seq_chrom = seq_chrom.strip('\n')
            else:
                index(line, region)
                write(region, seq_chrom, annoted_ref)


if __name__ == '__main__':
    main()
