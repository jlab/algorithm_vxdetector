#!/usr/bin/python

import os
from os.path import exists

def buildbowtie2(path):
	bowtie2_path = '/vol/software/bin/bowtie2-build'
	input_ref = path + 'Indexed_bt2/85_otus.fasta'	#reference genome greengenes is used
	output_path = path + 'Indexed_bt2/bowtie2'
	if exists(output_path + '.1.bt2'):
		pass	#only builds an index for the alignment if there isn't already one
	else:
		cmd = bowtie2_path + ' -f ' + input_ref + ' ' + output_path
		os.system(cmd)	#if a indexfile is missing a new index is build

def mapbowtie2(fasta_file, read2_file, path, temp_path, mode, file_type):
    bowtie2_path = '/vol/software/bin/bowtie2'
    samtools_path = '/usr/bin/samtools'
    index_path = path + 'Indexed_bt2/bowtie2'
    log_path = temp_path + 'bowtie2.log'
    if mode == 'unpaired':
        aligned_path = temp_path + 'unpaired.bam'
        cmd = bowtie2_path + ' -x ' + index_path + file_type + ' -U ' + fasta_file  + ' --fast 2> ' + log_path + ' | ' + samtools_path + ' view -b -S -F 4 -o ' + aligned_path
    else:
        aligned_path = temp_path + 'paired.bed'
        cmd = bowtie2_path + ' -x ' + index_path + ' -1 ' + fasta_file + ' -2 ' + read2_file + ' --fast 2> ' + log_path + ' | ' + samtools_path + ' view -b -S -F 4 | ' + samtools_path + ' sort -n  | /usr/bin/bedtools bamtobed -bedpe -i stdin > ' + aligned_path
    os.system(cmd)
    return aligned_path
