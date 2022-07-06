#!/usr/bin/python

import argparse
import os
import interact_bowtie2
import interact_bedtools
import Output_counter
import files_manager


def workflow(path, temp_path, file_path, file_type, file_name, dir_name,
             dir_path, mode, read2_file):
    interact_bowtie2.buildbowtie2(path)
    if file_type is not None:
        # The Programm bowtie2 is used to align the Reads to a reference
        # 16S database.
        aligned_path = interact_bowtie2.mapbowtie2(file_path, read2_file,
                                                   path, temp_path, mode,
                                                   file_type)
    else:
        aligned_path = interact_bowtie2.mapbowtie2(file_path, read2_file,
                                                   path, temp_path, mode,
                                                   file_type=' -q')
    # look which reads intersect with which variable Region
    interact_bedtools.overlap(path, temp_path, aligned_path)
    # counts the Variable Regions that are found with bedtools and prints the
    # highest probable variable Region
    Output_counter.count(temp_path, file_name, file_type, path,
                         dir_name, dir_path, mode)


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
    path = files_manager.get_lib()
    temp_path = files_manager.tmp_dir(path, temp_path='')
    file_type = None
    fasta_ext = ('.fasta', '.fa', '.ffa', '.ffn', '.faa', '.frn')
    fastq_ext = ('.fq', '.fastq',)
    if args.fasta_file is None:
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
                workflow(path, temp_path, file_path, file_type, file_name,
                         dir_name, args.dir_path, mode, read2_file)
        if file_type is None:
            print('There were no FASTA or FASTQ files with 16S sequences '
                  'in this directory')
    else:
        workflow(path, temp_path, args.fasta_file, file_type,
                 file_name='Your file', dir_name='', dir_path='',
                 mode='unpaired', read2_file='')
    # quit() #keeps the tmp folder for troubleshooting
    files_manager.tmp_dir(path, temp_path)


if __name__ == '__main__':
    main()
