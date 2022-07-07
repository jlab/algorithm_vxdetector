#!/usr/bin/python

import os
import csv
from itertools import (takewhile, repeat)

dictionary = {'V1': 0, 'V2': 0, 'V3': 0, 'V4': 0, 'V5': 0,
              'V6': 0, 'V7': 0, 'V8': 0, 'V9': 0}
# declares the dictionary containing all variable regions as
# a global variable. That way it can be easily accessed and modified


def directory_navi(file_name, path, dir_name, dir_path):
    dir_path = dir_name.replace(dir_path, '', 1)
    # only leaves the part in the directory tree between
    # the file and the original given directory path
    if dir_name.split('/')[-1] == '':
        dir_name = dir_name.split('/')[-2]
    else:
        dir_name = dir_name.split('/')[-1]
    dir_path = dir_path.replace(dir_name, '')
    dir_path = f'{path}Output/{dir_path}'
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


def region_count(unaligned_count, temp_path, mode, all_reads):
    if mode == 'paired':
        unpaired = round(rawincount(f'{temp_path}bed.log')/all_reads, 2)
        # counts the amount of error messages given by bedtools
        # during bam to bed conversion each error message reprensents
        # one improperly paired read
    elif mode == 'unpaired':
        unpaired = 'unpaired Reads'
    no_V = rawincount(f'{temp_path}noOver.bed')
    # counts the number of reads not mapped to any variable region
    aligned_count = 100 - unaligned_count
    count = sum([dictionary['V1'], dictionary['V2'], dictionary['V3'],
                dictionary['V4'], dictionary['V5'], dictionary['V6'],
                dictionary['V7'], dictionary['V8'], dictionary['V9'], no_V])
    for key in dictionary:
        dictionary[key] = (dictionary[key] / count) * aligned_count
    no_V = (no_V/count)*aligned_count
    # Calculates the percantage of reads mapped to either a specific
    # region or no region
    if aligned_count == 0 or (no_V/aligned_count) == 1:
        most_probable_V = 'No variable Region'
    else:
        most_probable_V = [x[0].replace('V', '') for x in
                           sorted(list(dictionary.items()))
                           if ((x[1] / aligned_count) * 100) > 20]
        most_probable_V = 'V' + ''.join(most_probable_V)
        if most_probable_V == 'V':
            most_probable_V = 'No variable Region'
    # finds the variable regions with the highest probability
    # any region having a percentage of 20% or higher will be listed.
    # The percentage is calculated within the regions.

    return most_probable_V, unpaired, no_V


def create_output(path, file_name, unaligned_count, dir_name,
                  dir_path, temp_path, all_reads, mode):
    file_name, dir_name, dir_path = directory_navi(file_name, path,
                                                   dir_name, dir_path)
    if mode == 'paired':
        file_name = file_name.replace('_R1_001', '')
    new_file = f'{dir_path}{dir_name}.csv'
    # The new file is named after the directory containing the reads
    most_probable_V, unpaired, no_V = region_count(unaligned_count, temp_path,
                                                   mode, all_reads)
    if os.path.exists(new_file):
        header = True
    else:
        header = False
    # The header will only be written if it wasnt written before
    with open(new_file, 'a', newline='') as o:
        fieldnames = ['Read-file', 'Number of Reads', 'Unaligned Reads [%]',
                      'Not properly paired', 'Sequenced variable region', 'V1',
                      'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9',
                      'Not aligned to a variable region']
        writer = csv.DictWriter(o, fieldnames=fieldnames)
        if header is False:
            writer.writeheader()
        for key in dictionary:
            dictionary[key] = round(dictionary[key], 2)
        writer.writerow({'Read-file': file_name, 'Number of Reads': all_reads,
                         'Unaligned Reads [%]': round(unaligned_count, 2),
                         'Not properly paired': unpaired,
                         'Sequenced variable region': most_probable_V,
                         'V1': dictionary['V1'], 'V2': dictionary['V2'],
                         'V3': dictionary['V3'], 'V4': dictionary['V4'],
                         'V5': dictionary['V5'], 'V6': dictionary['V6'],
                         'V7': dictionary['V7'], 'V8': dictionary['V8'],
                         'V9': dictionary['V9'],
                         'Not aligned to a variable region': round(no_V, 2)})


def print_output(unaligned_count, temp_path, mode):
    most_probable_V, unpaired, no_V = region_count(unaligned_count, temp_path,
                                                   mode, all_reads='')
    print(f'\n{str(unaligned_count)}% were unaligned.')
    print(f'Of all the aligned Reads most were aligned to: {most_probable_V}')
    print('The probabilities of all regions is as follows [%]:')
    print('\n'.join(':'.join((key, str(round(val, 2))))
          for (key, val) in dictionary.items()))
    print(f'\nAnd {round(no_V, 2)}% were not mapped to any variable region\n')


def count(temp_path, file_name, file_type, path, dir_name, dir_path, mode):
    BED_path = f'{temp_path}BED.bed'
    Log_path = f'{temp_path}bowtie2.log'
    with open(BED_path, 'r') as bed:
        for line in bed:
            line_list = line.split('\t')
            dictionary[line_list[-1].strip('\n')] += 1
            # counts all appearing variable Regions
    with open(Log_path, 'r') as log:
        lines = log.readlines()
        all_reads = int(lines[0].strip('\n\t ').split()[0])
        unaligned_count = lines[-1].strip('\n\t ')
        # Reads the bowtie2 stdout which is normally seen in the terminal
    unaligned_count = unaligned_count.split()[0].replace('%', '')
    unaligned_count = 100 - float(unaligned_count)
    if file_type is None:
        print_output(unaligned_count, temp_path, mode)
        # should file_type be None only a single file was given
        # therefore output will be printed in the terminal rather
        # than in a file
    else:
        create_output(path, file_name, unaligned_count,
                      dir_name, dir_path, temp_path, all_reads, mode)
