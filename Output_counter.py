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


def region_count(dictionary):
    count = sum([dictionary['V1'], dictionary['V2'],dictionary['V3'],dictionary['V4'],dictionary['V5'],dictionary['V6'],dictionary['V7'], dictionary['V8'], dictionary['V9']])
    V1 = round((dictionary['V1']/count)*100, 2)
    V2 = round((dictionary['V2']/count)*100, 2)
    V3 = round((dictionary['V3']/count)*100, 2)
    V4 = round((dictionary['V4']/count)*100, 2)
    V5 = round((dictionary['V5']/count)*100, 2)
    V6 = round((dictionary['V6']/count)*100, 2)
    V7 = round((dictionary['V7']/count)*100, 2)
    V8 = round((dictionary['V8']/count)*100, 2)
    V9 = round((dictionary['V9']/count)*100, 2)
    dictionary = sorted(dictionary.items(), key=lambda x:x[1], reverse=True)
    most_probable_V = dictionary[0][0]
    return most_probable_V, V1, V2, V3, V4, V5, V6, V7, V8, V9




def create_output(path, file_name, unaligned_count, dictionary, dir_name, dir_path):
    file_name, dir_name, dir_path = directory_navi(file_name, path, dir_name, dir_path)
    new_file = dir_path + dir_name +'.csv'
    most_probable_V, V1, V2, V3, V4, V5, V6, V7, V8, V9 = region_count(dictionary)
    if os.path.exists(new_file):
        header = True
    else:
        header = False
    with open(new_file, 'a', newline='') as o:
        fieldnames = ['Read-file', 'Unaligned Reads [%]', 'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'V7', 'V8', 'V9']
        writer = csv.DictWriter(o, fieldnames=fieldnames)
        if header == False:
            writer.writeheader()
        writer.writerow({'Read-file' : file_name, 'Unaligned Reads [%]' : unaligned_count, 'Sequenced variable region' : most_probable_V, 'V1': V1, 'V2': V2, 'V3' : V3, 'V4':V4, 'V5':V5, 'V6':V6, 'V7':V7, 'V8':V8, 'V9':V9})

def print_output(dictionary, unaligned_count):
    most_probable_V, V1, V2, V3, V4, V5, V6, V7, V8, V9 = region_count(dictionary)
    V1 = str(V1)
    V2 = str(V2)
    V3 = str(V3)
    V4 = str(V4)
    V5 = str(V5)
    V6 = str(V6)
    V7 = str(V7)
    V8 = str(V8)
    V9 = str(V9)
    unaligned_count = '\n' + str(unaligned_count) + '% were unaligned.'
    print(unaligned_count)
    print('Of all the aligned Reads most were aligned to: ' + most_probable_V)
    print('The probabilities of all regions is as follows [%]:')
    print('V1: '+ V1 + '\nV2: ' + V2 + '\nV3: ' + V3 + '\nV4: ' + V4 + '\nV5: ' + V5 + '\nV6: ' + V6 + '\nV7: '+ V7 + '\nV8: ' + V8 + '\nV9: ' + V9 + '\n')



def count(temp_path, file_name, file_type, path, dir_name, dir_path):
    dictionary = {'V1': 0, 'V2': 0, 'V3' : 0, 'V4': 0, 'V5': 0, 'V6':0, 'V7':0, 'V8':0, 'V9':0}
    BED_path = temp_path + 'BED.bed'
    Log_path = temp_path + 'bowtie2.log'
    with open(BED_path, 'r') as bed:
        for line in bed:
            line_list = line.split('\t')
            dictionary[line_list[-1].strip('\n')] += 1
            #counts all appearing variable Regions
    with open(Log_path, 'r') as log:
        lines = log.readlines()
        unaligned_count = lines[-1].strip('\n\t ')
    unaligned_count = unaligned_count.split()[0].replace('%', '')
    unaligned_count = round((100 - float(unaligned_count)), 2)
    if file_type == None:
        print_output(dictionary, unaligned_count)
    else:
        create_output(path, file_name, unaligned_count, dictionary, dir_name, dir_path)




