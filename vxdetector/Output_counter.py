#!/usr/bin/python

import os
import pandas as pd
from itertools import (takewhile, repeat)


new_row = pd.DataFrame({'Read-file': [''], 'Number of Reads': [0],
                        'Unaligned Reads [%]': [0],
                        'Not properly paired': ['unpaired'],
                        'Sequenced variable region': ['V'],
                        'V1': [0], 'V2': [0], 'V3': [0],
                        'V4': [0], 'V5': [0], 'V6': [0],
                        'V7': [0], 'V8': [0], 'V9': [0],
                        'Not aligned to a variable region': [0]})
# declares the dictionary containing all variable regions as
# a global variable. That way it can be easily accessed and modified


def directory_navi(file_name, path, dir_name, dir_path):
    dir_path = dir_name.replace(dir_path, '', 1)
    if dir_name.endswith('/'):
        dir_name = dir_name.rstrip('/')
    # only leaves the part in the directory tree between
    # the file and the original given directory path
    dir_name = os.path.basename(dir_name)
    dir_path = f'{path}Output/{os.path.dirname(dir_path)}'
    # This block should set dir_path in such a way that
    # the directory structure found in the original given
    # directory is mirrored
    os.makedirs(dir_path, exist_ok=True)
    if file_name.endswith('.gz'):
        file_name = file_name.rsplit('.', 2)[0]
    else:
        file_name = file_name.rsplit('.', 1)[0]
    return file_name, dir_name, dir_path


def rawincount(filename):
    f = open(filename, 'rb')
    bufgen = takewhile(lambda x: x, (f.raw.read(1024*1024)
                                     for _ in repeat(None)))
    return sum(buf.count(b'\n') for buf in bufgen)


def region_count(temp_path, mode):
    if mode == 'paired':
        new_row['Not properly paired'] = rawincount(
                                         f'{temp_path}bed.log') \
                                         / new_row['Number of Reads']
        # counts the amount of error messages given by bedtools
        # during bam to bed conversion each error message reprensents
        # one improperly paired read
    new_row['Not aligned to a variable region'] = rawincount(
                                                  f'{temp_path}noOver.bed')
    # counts the number of reads not mapped to any variable region
    aligned_count = 100 - new_row.iat[0, 2]
    count = sum(new_row.iloc[0, 5:])
    for column in new_row.columns[5:]:
        new_row[column] = (new_row[column] / count) * aligned_count
    # Calculates the percantage of reads mapped to either a specific
    # region or no region
    if aligned_count == 0 or (new_row.iat[0, 14]/aligned_count) == 1:
        new_row['Sequenced variable region'] = 'No variable Region'
    else:
        most_probable_V = [x.replace('V', '') for x in
                           new_row.columns[5:14]
                           if ((new_row.at[0, x] / aligned_count) * 100) > 20]
        new_row['Sequenced variable region'] = 'V' + ''.join(most_probable_V)
        if new_row.iat[0, 4] == 'V':
            new_row['Sequenced variable region'] = 'No variable Region'
    # finds the variable regions with the highest probability
    # any region having a percentage of 20% or higher will be listed.
    # The percentage is calculated within the regions.


def create_row(path, file_name, dir_name, dir_path, temp_path, mode):
    file_name, dir_name, dir_path = directory_navi(file_name, path,
                                                   dir_name, dir_path)
    if mode == 'paired':
        new_row['Read-file'] = file_name.replace('_R1_001', '')
    new_file = f'{dir_path}/{dir_name}.csv'
    # The new file is named after the directory containing the reads
    region_count(temp_path, mode)
    for column in new_row:
        if isinstance(new_row.at[0, column], float):
            new_row.at[0, column] = round(new_row.at[0, column], 2)
    return new_file


def print_output(temp_path, mode):
    region_count(temp_path, mode)
    print(f"\n{round(new_row.iat[0, 2], 2)}% were unaligned.")
    print(f"Of all the aligned Reads most were aligned to: \
          {new_row.iat[0, 4]}")
    print('The probabilities of all regions is as follows [%]:')
    print('\n'.join(':'.join((column, str(round(new_row.at[0, column], 2))))
          for column in new_row.columns[5:14]))
    print(f"\nAnd {round(new_row.iat[0, -1], 2)}% were not mapped \
          to any variable region\n")


def count(temp_path, file_name, file_type, path, dir_name, dir_path, mode):
    BED_path = f'{temp_path}BED.bed'
    Log_path = f'{temp_path}bowtie2.log'
    with open(BED_path, 'r') as bed:
        for line in bed:
            line_list = line.split('\t')
            new_row[line_list[-1].strip('\n')] += 1
            # counts all appearing variable Regions
    with open(Log_path, 'r') as log:
        lines = log.readlines()
        new_row['Number of Reads'] = int(lines[0].strip('\n\t ').split()[0])
        unaligned_count = lines[-1].strip('\n\t ')
        # Reads the bowtie2 stdout which is normally seen in the terminal
    unaligned_count = unaligned_count.split()[0].replace('%', '')
    new_row['Unaligned Reads [%]'] = 100 - float(unaligned_count)
    if file_type is None:
        print_output(temp_path, mode)
        # should file_type be None only a single file was given
        # therefore output will be printed in the terminal rather
        # than in a file
    else:
        new_file = create_row(path, file_name, dir_name, dir_path,
                              temp_path, mode)
    if file_type is not None:
        return new_row, new_file
