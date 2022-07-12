#!/usr/bin/python

import argparse
import os
import interact_bowtie2
import interact_bedtools
import Output_counter
import files_manager
import pandas as pd

output = pd.DataFrame({'Read-file': [], 'Number of Reads': [],
                       'Unaligned Reads [%]': [],
                       'Not properly paired': [],
                       'Sequenced variable region': [],
                       'V1': [], 'V2': [], 'V3': [],
                       'V4': [], 'V5': [], 'V6': [],
                       'V7': [], 'V8': [], 'V9': [],
                       'Not aligned to a variable region': []})


def do_statistic():
    global output
    avarage = output.mean(numeric_only=True).round(2).to_frame().T
    avarage['Read-file'] = ['Average']
    avarage['Sequenced variable region'] = (output.iloc[:, [4]].
                                            mode().values.tolist())
    std_dev = output.std(numeric_only=True).round(2).to_frame().T
    std_dev['Read-file'] = ['Standard deviation']
    statistic = pd.concat([avarage, std_dev], axis=0)
    statistic = statistic[['Read-file', 'Number of Reads',
                           'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3',
                           'V4', 'V5', 'V6', 'V7', 'V8', 'V9',
                           'Not aligned to a variable region']]
    output = pd.concat([statistic, output], axis=0)


def workflow(path, temp_path, file_path, file_type, old_file_name,
             file_name, dir_name, dir_path, mode, read2_file):
    interact_bowtie2.buildbowtie2(path)
    if file_type is not None:
        aligned_path = interact_bowtie2.mapbowtie2(file_path, read2_file,
                                                   path, temp_path, mode,
                                                   file_type)
    else:
        aligned_path = interact_bowtie2.mapbowtie2(file_path, read2_file,
                                                   path, temp_path, mode,
                                                   file_type=' -q')
    # The Programm bowtie2 is used to align the Reads to a reference
    # 16S database.
    interact_bedtools.overlap(path, temp_path, aligned_path)
    # look which reads intersect with which variable Region
    if file_type is not None:
        global output
        new_row, new_file = Output_counter.count(temp_path, file_name,
                                                 file_type, path, dir_name,
                                                 dir_path, mode)
        if new_file != old_file_name and old_file_name is not None:
            output.sort_values(by=['Read-file'], inplace=True)
            do_statistic()
            output.to_csv(old_file_name, index=False)
            output = output.iloc[0:0]
        output = pd.concat([output, new_row], axis=0)
        old_file_name = new_file
        return new_file, old_file_name
    else:
        Output_counter.count(temp_path, file_name, file_type, path,
                             dir_name, dir_path, mode)
    # counts the Variable Regions that are found with bedtools and prints the
    # highest probable variable Region (single file) or
    # writes a new file (directory)


def main():
    parser = argparse.ArgumentParser(prog='VX detector', description=(
        'This programm tries to find which variable region of the 16S '
        'sequence was sequencend'))
    parser.add_argument('-d', '--directory', dest='dir_path',
                        help=('Directory path of the directory containing '
                              'multiple fastq or fasta files.'))
    parser.add_argument('-sf', '--single-file', dest='fasta_file',
                        default=None,
                        help=('Filepath of the fastq file containing the '
                              'sequencences of which the variable regions '
                              'are to be verified.'))
    args = parser.parse_args()
    # allows terminal input
    path = files_manager.get_lib()
    temp_path = files_manager.tmp_dir(path, temp_path='')
    # sets the paths of the programm itself and a temporary folder
    file_type = None
    fasta_ext = ('.fasta', '.fa', '.ffa', '.ffn', '.faa', '.frn')
    fastq_ext = ('.fq', '.fastq',)
    old_file_name = None
    if args.fasta_file is None:
        # If a directory was given as input the programm will walk through
        # that directory and search for forward reads.
        # If found it will look for its complementary backward read.
        for root, dirs, files in os.walk(args.dir_path, topdown=True):
            for file in files:
                mode = 'unpaired'
                read2_file = ''
                if '_R2_' in file:
                    continue
                if any(elements in file for elements in fastq_ext):
                    file_name = file
                    read2_file = os.path.join(root, file.replace('_R1_',
                                                                 '_R2_'))
                    rev_exists = os.path.exists(read2_file)
                    if '_R1_' in file_name and rev_exists is True:
                        mode = 'paired'
                    dir_name = root
                    file_path = os.path.join(root, file)
                    file_type = ' -q'
                elif any(elements in file for elements in fasta_ext):
                    file_name = file
                    read2_file = os.path.join(root, file.replace('_R1_',
                                                                 '_R2_'))
                    rev_exists = os.path.exists(read2_file)
                    if '_R1_' in file_name and rev_exists is True:
                        mode = 'paired'
                    dir_name = root
                    file_path = os.path.join(root, file)
                    file_type = ' -f'
                else:
                    continue
                new_file, old_file_name = workflow(path, temp_path, file_path,
                                                   file_type, old_file_name,
                                                   file_name, dir_name,
                                                   args.dir_path, mode,
                                                   read2_file)
        if file_type is None:
            print('There were no FASTA or FASTQ files with 16S sequences '
                  'in this directory')
    else:
        workflow(path, temp_path, args.fasta_file, file_type, old_file_name,
                 file_name='Your file', dir_name='', dir_path='',
                 mode='unpaired', read2_file='')
    files_manager.tmp_dir(path, temp_path)
    if old_file_name is not None:
        do_statistic()
        output.to_csv(old_file_name, index=False)


if __name__ == '__main__':
    main()
