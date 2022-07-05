#!/usr/bin/python

import os
import csv


def directory_navi(file_name, path, dir_name, dir_path):
    dir_path = dir_name.replace(dir_path, '', 1)
    dir_name = dir_path.split('/')[-1]
    dir_path = dir_path.replace(dir_name, '')
    dir_path = path+'Output/'+dir_path
    os.makedirs(dir_path, exist_ok=True)
    if file_name.endswith('.gz'):
        file_name = file_name.rsplit('.', 2)[0]
    else:
        file_name = file_name.rsplit('.', 1)[0]

    return file_name, dir_name, dir_path


def create_output(path, file_name, unaligned_count, most_probable_V,
                  probability, dir_name, dir_path):
    file_name, dir_name, dir_path = directory_navi(file_name, path,
                                                   dir_name, dir_path)
    unaligned_count = unaligned_count.split(' ')[1]
    unaligned_count = float(unaligned_count.strip('()%'))
    new_file = dir_path + dir_name + '.csv'
    if os.path.exists(new_file):
        header = True
    else:
        header = False
    with open(new_file, 'a', newline='') as o:
        fieldnames = ['Read-file', 'Unaligned Reads [%]',
                      'Sequenced variable region', 'Probability [%]']
        writer = csv.DictWriter(o, fieldnames=fieldnames)
        if header is False:
            writer.writeheader()
        writer.writerow({
            'Read-file': file_name,
            'Unaligned Reads [%]': unaligned_count,
            'Sequenced variable region': most_probable_V,
            'Probability [%]': probability})


def count(count, temp_path, file_name, file_type, path, dir_name, dir_path):
    dictionary = {}
    BED_path = temp_path + 'BED.bed'
    Log_path = temp_path + 'bowtie2.log'
    with open(BED_path, 'r') as bed:
        for line in bed:
            line_list = line.split('\t')
            if line_list[-1] not in dictionary:
                dictionary[line_list[-1]] = 1
            else:
                dictionary[line_list[-1]] += 1
                # counts all appearing variable Regions
    # sorts the Dictionary with all variable regions for extracting the most
    # probable one
    dictionary = sorted(dictionary.items(), key=lambda x: x[1], reverse=True)
    most_probable_V = dictionary[0]
    most_probable_V_count = most_probable_V[1]
    most_probable_V = most_probable_V[0].strip('\n')
    probability = round((most_probable_V_count/count)*100, 2)
    with open(Log_path, 'r') as log:
        lines = log.readlines()
        unaligned_count = lines[2].strip('\n\t ')
    if file_type is None:
        print(file_name + ': ')
        print(unaligned_count)
        print('The sequenced variable Region in this file is '
              + most_probable_V + ' with a probability of '
              + str(probability) + '%.')
    else:
        create_output(path, file_name, unaligned_count, most_probable_V,
                      probability, dir_name, dir_path)
