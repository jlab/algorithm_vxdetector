#!/usr/bin/python

# This programm takes the variable region reference position
# from variable_region_calc.py and searches for them in the
# aligned gg otus

from os import path as p


def write_reference_file(seq_chrom, annoted_ref, region):
    r'''Reference file writer

    This function uses the dictionary created by index and
    writes a file containing the start and end position of
    each variable region for each entry in the MSA refrence.

    Parameters
    ----------
    seq_chrom: str
        Name of the entry

    annoted_ref: str
        Filepath where the created Refrence file should
        be saved to.

    region: dict
        Dictionary conataining region start / end as
        keys and the position in the refrence file
        declared in main with "file".

    '''
    with open(annoted_ref, 'a') as t:
        if 'V1_start' in region:
            for counter, key in enumerate(region):
                if (counter % 2) == 0:
                    wordstart = f'{seq_chrom}\t{str(region[key])}'
                else:
                    var_reg = key[0:2]
                    wordend = f'\t{str(region[key])}\t{var_reg}\n'
                    words = wordstart + wordend
                    t.write(words)
    region.clear()


def index(line, region):
    r'''Variable region position finder

    This function searches for the positions of
    the 16S variable regions within a line from
    the MSA refrence file via deleting the gaps.

    Parameters
    ----------
    line: str
        String of a line from the efrence file
        declared in main with "file".

    region: dict
        Empty dictionary.

    Returns
    -------
    region: dict
        Dictionary conataining region start / end as
        keys and the position in the refrence file
        declared in main with "file".

    '''
    boundary = {'V1_start': 189, 'V1_end': 471,
                'V2_start': 485, 'V2_end': 1867,
                'V3_start': 1915, 'V3_end': 2231,
                'V4_start': 2262, 'V4_end': 4050,
                'V5_start': 4089, 'V5_end': 4521,
                'V6_start': 4652, 'V6_end': 4931,
                'V7_start': 5043, 'V7_end': 5806,
                'V8_start': 5909, 'V8_end': 6426,
                'V9_start': 6449, 'V9_end': 6790}
    '''
    The boundary positions were found using code provided by Tony Walters
    (http://qiime.org/home_static/nih_cloud-apr2012/variable_region_position_calculations.pdf).
    It was then modified to allow for the additional search of the
    V1 rev, V2 fwd, V5, V7 and V8 boundary positions.
    Primers used:
    V1: 27f / P2 (original code / Cocolin et al. 2001)
    V2: V2f / 338r (modified V2-V3 fwd primer from "16S V2-V3 Library
                    Preparation Kit for Illumina" / original code)
    V3: 349f / 534r (original code)
    V4: 515f / 806r (original code)
    V5: 806f / 926r (kindly provided by Tony Walters)
    V6: 967f / 1046r (original code)
    V7: 1115f / 1193r (Schneyder et al. 2021 / Bodenhausen et al. 2013)
    V8: 1237f / 1291r (Turner et al. 1999)
    V9: 1391f / 1492r (original code)
    '''
    for key in boundary:
        list1 = list(line)
        list1[boundary[key]] = '1'
        list1 = ''.join(str(e) for e in list1)
        list1 = list1.replace('-', '')
        region[key] = int(list1.index('1')) + 1
    return region


def main():
    r'''Statistical analysis of given directory

    This program is used to create a reference file
    containing the positions of 16S variable regions for
    VXdetector.

    '''
    file_ = f'{p.dirname(p.dirname(__file__))}/85_otus_aligned.fasta'
    annoted_ref = f'{p.dirname(p.dirname(__file__))}/annoted_ref.bed'
    with open(file_, 'r') as f:
        region = dict()
        for line in f:
            if line.startswith('>'):
                seq_chrom = line[1:]
                seq_chrom = seq_chrom.strip('\n')
            else:
                region = index(line, region)
                write_reference_file(seq_chrom, annoted_ref, region)


if __name__ == '__main__':
    main()
