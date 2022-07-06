#!/usr/bin/python

import os
import csv
from itertools import (takewhile,repeat)

def directory_navi(file_name, path, dir_name, dir_path):
    dir_path = dir_name.replace(dir_path, '', 1)
    if dir_name.split('/')[-1] == '':
        dir_name = dir_name.split('/')[-2]
    else:
        dir_name = dir_name.split('/')[-1]
    dir_path = dir_path.replace(dir_name, '')
    dir_path = f'{path}Output/{dir_path}'
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

    return sum( buf.count(b'\n') for buf in bufgen )


def region_count(dictionary, unaligned_count ,temp_path, all_reads, mode):
    if temp_path is not None and mode == 'paired':
        unpaired = round((rawincount(f'{temp_path}bed.log')/all_reads), 2)
    elif mode == 'unpaired':
        unpaired = 'unpaired Reads'
    no_V = rawincount(f'{temp_path}noOver.bed')
    aligned_count = 100 - unaligned_count
    count = sum([dictionary['V1'], dictionary['V2'],dictionary['V3'],
		 dictionary['V4'],dictionary['V5'],dictionary['V6'],dictionary['V7'],
		 dictionary['V8'], dictionary['V9'], no_V])
    dictionary['V1'] = round((dictionary['V1']/count)* aligned_count, 2)
    dictionary['V2'] = round((dictionary['V2']/count)* aligned_count, 2)
    dictionary['V3'] = round((dictionary['V3']/count)* aligned_count, 2)
    dictionary['V4'] = round((dictionary['V4']/count)*aligned_count, 2)
    dictionary['V5'] = round((dictionary['V5']/count)* aligned_count, 2)
    dictionary['V6'] = round((dictionary['V6']/count)*aligned_count, 2)
    dictionary['V7'] = round((dictionary['V7']/count)*aligned_count, 2)
    dictionary['V8'] = round((dictionary['V8']/count)*aligned_count, 2)
    dictionary['V9'] = round((dictionary['V9']/count)*aligned_count, 2)
    no_V = round((no_V/count)*aligned_count, 2)
    if aligned_count == 0 or (no_V/aligned_count) == 1 :
        most_probable_V = 'No variable Region'
    else:
        most_probable_V = [x[0].replace('V', '') for x in sorted(list(dictionary.items())) \
			   if ((x[1]/aligned_count) *100)>20]
        most_probable_V = 'V' + ''.join(most_probable_V)
        if most_probable_V == 'V':
            most_probable_V = 'No variable Region'

    return most_probable_V, dictionary, unpaired, no_V



def create_output(path, file_name, unaligned_count, dictionary, dir_name,
		  dir_path, temp_path, all_reads, mode):
    file_name, dir_name, dir_path = directory_navi(file_name, path,
						   dir_name, dir_path)
    if mode == 'paired':
        file_name = file_name.replace('_R1_001', '')
    new_file = f'{dir_path}{dir_name}.csv'
    most_probable_V, dictionary, unpaired, no_V = \
    region_count(dictionary, unaligned_count, temp_path, all_reads, mode)
    if os.path.exists(new_file):
        header = True
    else:
        header = False
    with open(new_file, 'a', newline='') as o:
        fieldnames = ['Read-file', 'Number of Reads' , 'Unaligned Reads [%]',
		      'Not properly paired', 'Sequenced variable region', 'V1',
		      'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9', 
		      'Not aligned to a variable region']
        writer = csv.DictWriter(o, fieldnames=fieldnames)
        if header is False:
            writer.writeheader()
        writer.writerow({'Read-file' : file_name, 'Number of Reads': all_reads,
			 'Unaligned Reads [%]' : unaligned_count,
			 'Not properly paired': unpaired,
			 'Sequenced variable region' : most_probable_V, 
			 'V1': dictionary['V1'], 'V2': dictionary['V2'],
			 'V3' : dictionary['V3'], 'V4':dictionary['V4'],
			 'V5':dictionary['V5'], 'V6':dictionary['V6'],
			 'V7':dictionary['V7'], 'V8':dictionary['V8'],
			 'V9':dictionary['V9'], 'Not aligned to a variable region':no_V})


def print_output(dictionary, unaligned_count, number_aligned):
    most_probable_V, dictionary = region_count(dictionary, unaligned_count,
					       temp_path=None, all_reads='', mode='')
    dictionary['V1'] = str(dictionary['V1'])
    dictionary['V2'] = str(dictionary['V2'])
    dictionary['V3'] = str(dictionary['V3'])
    dictionary['V4'] = str(dictionary['V4'])
    dictionary['V5'] = str(dictionary['V5'])
    dictionary['V6'] = str(dictionary['V6'])
    dictionary['V7'] = str(dictionary['V7'])
    dictionary['V8'] = str(dictionary['V8'])
    dictionary['V9'] = str(dictionary['V9'])
    print(f'\n{str(unaligned_count)}% were unaligned.')
    print(f'Of all the aligned Reads most were aligned to: {most_probable_V}')
    print('The probabilities of all regions is as follows [%]:')
    print(f'V1: {V1}\nV2: {V2}\nV3: {V3}\nV4: {V4}\nV5: {V5}\n \
	  V6: {V6}\nV7: {V7}\nV8: {V8}\nV9: {V9}\n')



def count(temp_path, file_name, file_type, path, dir_name, dir_path, mode):
    dictionary = {'V1': 0, 'V2': 0, 'V3' : 0, 'V4': 0, 'V5': 0,
		  'V6':0, 'V7':0, 'V8':0, 'V9':0}
    BED_path = f'{temp_path}BED.bed'
    Log_path = f'{temp_path}bowtie2.log'
    with open(BED_path, 'r') as bed:
        for line in bed:
            line_list = line.split('\t')
            dictionary[line_list[-1].strip('\n')] += 1
            #counts all appearing variable Regions
    with open(Log_path, 'r') as log:
        lines = log.readlines()
        all_reads = int(lines[0].strip('\n\t ').split()[0])
        unaligned_count = lines[-1].strip('\n\t ')
    unaligned_count = unaligned_count.split()[0].replace('%', '')
    unaligned_count = round((100 - float(unaligned_count)), 2)
    if file_type == None:
        print_output(dictionary, unaligned_count, number_aligned)
    else:
        create_output(path, file_name, unaligned_count, dictionary, 
		      dir_name, dir_path, temp_path, all_reads, mode)

