#!/usr/bin/python

import os
import subprocess

def SAMtoBAM(path, temp_path):
    samtools_filepath = '/usr/bin/samtools'
    BAM_filepath = temp_path + 'BAM.bam'
    BAM_sorted_filepath = temp_path + 'bam_sorted.bam'
    cmd = samtools_filepath + ' sort ' + BAM_filepath + ' -o ' + BAM_sorted_filepath #sorts the .bam file
    os.system(cmd)
    cmd = [samtools_filepath,'view','-c',BAM_sorted_filepath]
    count = subprocess.run(cmd, stdout=subprocess.PIPE)
    count = int(count.stdout.decode('utf-8'))
    cmd = samtools_filepath + ' index ' + BAM_sorted_filepath	#indexes the .bam file.
    os.system(cmd)
    
    return count

"""   
def count(path, temp_path):
    samtools_filepath = '/usr/bin/samtools'
    BED_filepath = temp_path + 'BED.bed'
    S_ref = path + 'Indexed_bt2/annoted_ref.bed'
    BAM_filepath = temp_path + 'bam_sorted.bam'
    pileup_path = temp_path + 'pileup.txt'
    pileup_test = temp_path + 'test.txt'
    cmd = samtools_filepath + ' mpileup ' + BAM_filepath + ' -l ' + S_ref + ' -o ' + pileup_path + ' 2> /dev/null'
    os.system(cmd)
    with open(pileup_path, 'r') as pu:
        for line in pu:
            line_list = line.split('\t')
            with open(pileup_test, 'a') as t:
                words = line_list[0] + '\t' + line_list[1] + '\t' + line_list[3] + '\n'
                t.write(words)
"""
