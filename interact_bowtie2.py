#!/usr/bin/python

import os

def buildbowtie2(path):
	bowtie2_path = '/vol/software/bin/bowtie2-build'
	input_ref =  f' {path}Indexed_bt2/85_otus.fasta'	#reference genome greengenes is used
	output_path = f'{path}Indexed_bt2/bowtie2'
	if os.path.exists(f'{output_path}.1.bt2'):
		pass	#only builds an index for the alignment if there isn't already one
	else: 
		os.system(f'{bowtie2_path} -f {input_ref} {output_path}')	#if a indexfile is missing a new index is build

def mapbowtie2(fasta_file, read2_file, path, temp_path, mode, file_type):
    bowtie2_path = '/vol/software/bin/bowtie2'
    samtools_path = '/usr/bin/samtools'
    bedtools_path = '/usr/bin/bedtools'
    index_path = f'{path}Indexed_bt2/bowtie2'
    log_path = f'{temp_path}bowtie2.log'
    bed_logpath = f'{temp_path}bed.log'
    if mode == 'unpaired':
        aligned_path = f'{temp_path}unpaired.bam'
        cmd = f'{bowtie2_path} -x {index_path} {file_type} -U {fasta_file} --fast 2> {log_path} | {samtools_path} view -b -S -F 4 -o {aligned_path}'
    else:
        aligned_path = f'{temp_path}paired.bed'
        cmd = f'{bowtie2_path} -x {index_path} -1 {fasta_file} -2 {read2_file} --fast 2> {log_path} | {samtools_path} view -b -S -F 4 | {bedtools_path} bamtobed -bedpe -i stdin > {aligned_path} 2> {bed_logpath}'
    os.system(cmd)
    return aligned_path
