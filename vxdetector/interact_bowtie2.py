#!/usr/bin/python

import os
import shutil
from Output_counter import directory_navi
# necessery for the delevopmental codeblock at the end

bowtie2_path = shutil.which('bowtie2')
samtools_path = shutil.which('samtools')
bedtools_path = shutil.which('bedtools')


def buildbowtie2(path):
    input_ref = f'{path}Indexed_bt2/85_otus.fasta'
    # reference genome greengenes is used
    output_path = f'{path}Indexed_bt2/bowtie2'
    if os.path.exists(f'{output_path}.1.bt2'):
        # only builds an index for the alignment if there isn't already one
        pass
    else:
        # if a indexfile is missing a new index is build
        os.system(f'{bowtie2_path}-build -f {input_ref} {output_path}')


def mapbowtie2(fasta_file, read2_file, path, temp_path, mode, file_type):
    index_path = f'{path}Indexed_bt2/bowtie2'
    if os.path.exists(f'{index_path}.1.bt2') is False:
        raise Exception(f'No Index files found under "{index_path}"')
    log_path = f'{temp_path}bowtie2.log'
    bed_logpath = f'{temp_path}bed.log'
    test = f'/homes/jgroos/Documents/{os.path.basename(fasta_file)}'
    # declares various filepaths
    Error = False
    if mode == 'unpaired':
        aligned_path = f'{temp_path}unpaired.bam'
        cmd = f'{bowtie2_path} -x {index_path} {file_type} -U {fasta_file} \
                 --fast 2> {log_path} | {samtools_path} view -b -q 25 -S -F 4 \
                 -o {aligned_path}'
        # Should no backward read be found it will just use the forward
        # read and does an alignment followed by a pipe to convert
        # the bowtie2 output .sam to a .bam file.
    else:
        aligned_path = f'{temp_path}paired.bed'
        cmd = f'{bowtie2_path} -x {index_path} -1 {fasta_file} -2 {read2_file} \
                --fast 2> {log_path} | {samtools_path} view -b -q 25 -S -F 4 | \
                {bedtools_path} bamtobed -bedpe -i stdin > {aligned_path} \
                2> {bed_logpath}'
        # Should a backward read be found both files will be given to bowtie2.
        # After converting .sam to .bam a conversion to .bed is done to
        # properly mate the pairs.
    os.system(cmd)
    with open(log_path, 'r') as log:
        lines = log.readlines()
        if lines[0].startswith('stat: Bad file descriptor') or lines[0] == '':
            wrong_file = lines[1].split('"')[1]
            print(f'{wrong_file} does not look like a compatible file. Skipping.')
            Error = True
    

    return aligned_path, Error

'''
# This Block is a developmental alternative to the mapbowtie2 function above
# which writes a .sam file instead of piping it directly to a .bam converter.
# There is a need to change the given parameters in VXdetector.py
def mapbowtie2(fasta_file, read2_file, path, temp_path, mode,
               file_name, dir_name, dir_path, file_type):
    file_name, dir_name, dir_path = directory_navi(file_name, path,
                                                   dir_name, dir_path)
    SAM_path = f'{dir_path}{file_name}.sam'
    bowtie2_path = '/vol/software/bin/bowtie2'
    index_path = f'{path}Indexed_bt2/bowtie2'
    log_path = f'{temp_path}bowtie2.log'
    bed_logpath = f'{temp_path}bed.log'
    if mode == 'unpaired':
        aligned_path = f'{temp_path}unpaired.bam'
        cmd = f'{bowtie2_path} -x {index_path} {file_type} -U {fasta_file} \
                --fast -S {SAM_path} 2> {log_path}'
        cmd2 = f'{samtools_path} view -b -S -F 4 {SAM_path} -o {aligned_path}'
    else:
        aligned_path = f'{temp_path}paired.bed'
        cmd = f'{bowtie2_path} -x {index_path} -1 {fasta_file} -2 \
              {read2_file} --fast -S {SAM_path} 2> {log_path}'
        cmd2 = f'{samtools_path} view -b -S -F 4  {SAM_path} | \
                 {bedtools_path} bamtobed -bedpe -i \
                 stdin > {aligned_path} 2> {bed_logpath}'
    os.system(cmd)
    os.system(cmd2)
    return aligned_path
'''
