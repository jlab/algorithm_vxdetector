#!/usr/bin/python

from itertools import (takewhile, repeat)
from os.path import exists

regions_list = ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9']
# lists all variable regions for iteration


def rawincount(filename):
    r'''Lines counter

    This function counts the amount of lines within a file.

    Parameters
    ----------
    filename : str
        Filepath to the file which is to be counted.

    Returns
    -------
    lines_count : int
        Amount of lines in the counted file.

    '''
    if exists(filename) is False:
        raise FileNotFoundError(f'It seems {filename} is missing.'
                                ' Check if you deleted the temp_dir')
    with open(filename, 'rb') as f:
        bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024)
                                         for _ in repeat(None)))
        lines_count = sum(buf.count(b'\n') for buf in bufgen)
    return lines_count


def region_count(temp_path, paired, new_row, regions):
    r'''Create dictionary

    This function streamlines already obtained information and
    information generated on its own into a dictionary.

    Parameters
    ----------
    temp_path : str
        Path to the directory containing the noOver.bed
        file in which all reads not intersecting with a
        variable region are listed in and the bed.log which
        is an indicator for how many reads were not properly
        paired.
    paired : Bool
        Wether or not the reads in the directory have paired mates
        or are unpaired.
    new_row : dict
        Dictionary containing all analysed information about current read file.
        keys: ['Number of Reads', 'Unaligned Reads [%]']
    regions : dict
        Dictionary containing the occurences of all variable regions
        within the intersect file.
        keys: ['V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9']

    Returns
    -------
    new_row : dict
        Dictionary containing all analysed information about current read file.
        keys: ['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
               'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
               'V7', 'V8', 'V9', 'Not aligned to a variable region']

    '''
    if paired is True:
        new_row['Not properly paired'] = rawincount(
                                         f'{temp_path}bed.log') \
                                         / new_row['Number of Reads']
        # counts the amount of error messages given by bedtools
        # during bam to bed conversion each error message reprensents
        # one improperly paired read
    else:
        new_row['Not properly paired'] = 'not paired'
    regions['Not aligned to a variable region'] = rawincount(
                                                  f'{temp_path}noOver.bed')
    # counts the number of reads not mapped to any variable region
    aligned_count = 100 - new_row['Unaligned Reads [%]']
    var_re_count = sum(regions.values())
    if aligned_count == 0:
        new_row['Sequenced variable region'] = 'No variable Region'
    else:
        for region in regions:
            regions[region] = (regions[region] / var_re_count) * aligned_count
        # Calculates the percantage of reads mapped to either a specific
        # region or no region
        most_probable_V = [x.replace('V', '') for x in
                           regions if ((regions[x] / aligned_count) * 100) > 20
                           and x.startswith('V')]
        new_row['Sequenced variable region'] = 'V' + ''.join(most_probable_V)
        if new_row['Sequenced variable region'] == 'V':
            new_row['Sequenced variable region'] = 'No variable Region'
    # finds the variable regions with the highest probability
    # any region having a percentage of 20% or higher will be listed.
    # The percentage is calculated within the regions.
    new_row.update(regions)
    return new_row


def create_row(temp_path, paired):
    r'''Base function to create Output

    This function counts the occurences of the variable regions and
    reads the total number of reads and unaligned reads from the bowtie2.log.
    Then it calls other functions which create a dictionary containing
    analised data.

    Parameters
    ----------
    temp_path : str
        Path to the directory containing the bowtie2.log and the .bed
        file in which the intersections between mapped reads and
        variable region reference are listed.
    paired : Bool
        Wether or not the reads in the directory have paired mates
        or are unpaired.

    Returns
    -------
    new_row : dict
        Dictionary containing all analysed information about current read file.
        keys: ['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
               'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
               'V7', 'V8', 'V9', 'Not aligned to a variable region']

    '''
    BED_path = f'{temp_path}BED.bed'
    Log_path = f'{temp_path}bowtie2.log'
    if exists(BED_path) is False:
        raise FileNotFoundError(f'It seems {BED_path} is missing.'
                                ' Check if you deleted the temp_dir')
    if exists(Log_path) is False:
        raise FileNotFoundError(f'It seems {Log_path} is missing.'
                                ' Check if you deleted the temp_dir')
    new_row = dict()
    regions = dict()
    b = open(BED_path, 'r')
    bed_data = b.read()
    b.close()
    for region in regions_list:
        regions[region] = bed_data.count(region)
        # counts all appearing variable Regions
    with open(Log_path, 'r') as log:
        lines = log.readlines()
        new_row['Number of Reads'] = int(lines[0].strip('\n\t ').split()[0])
        unaligned_count = lines[-1].strip('\n\t ')
        # Reads the bowtie2 stdout which is normally seen in the terminal
    unaligned_count = unaligned_count.split()[0].replace('%', '')
    new_row['Unaligned Reads [%]'] = 100 - float(unaligned_count)
    new_row = region_count(temp_path, paired, new_row, regions)

    return new_row
