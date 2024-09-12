#!/usr/bin/python

import argparse
import glob
import os
import sys
import warnings
import pandas as pd
import Output_counter as Output_counter
import files_manager as files_manager
from interact_bowtie2 import mapbowtie2, buildbowtie2
from interact_bedtools import overlap
from multiprocessing import pool
import random

def sample_fastq(file_path, sample_size, sampled_indices=None):
    '''Get random parts from the FASTQ file based on shared indices for paired-end reads.'''
    sampled_reads = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        total_reads = len(lines) // 4
        if sampled_indices is None:
            sampled_indices = sorted(random.sample(range(total_reads), sample_size))
        for idx in sampled_indices:
            read = lines[idx*4:(idx+1)*4]
            sampled_reads.extend(read)
    return sampled_reads, sampled_indices

def save_sampled_fastq(sampled_reads, output_path):
    '''Save sampled reads to a temporary FASTQ file'''
    with open(output_path, 'w') as f:
        f.writelines(sampled_reads)

def do_statistic(result):

    average = result.mean(numeric_only=True).to_frame().T
    region = (result['Sequenced variable region'].mode().values)
    region = ' / '.join(str(r) for r in region)
    region = region.replace('\'', '')
    region = region.replace('[', '')
    region = region.replace(']', '')
    average['Sequenced variable region'] = region
    if 'Not properly paired' not in average.columns:
        average['Not properly paired'] = 'not paired'
    std_dev = result.std(numeric_only=True).to_frame().T
    statistic = pd.concat([average, std_dev], axis=0)
    statistic = statistic[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                           'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    statistic['row_descriptor'] = ['Average', 'Standard deviation']
    statistic = statistic.set_index('row_descriptor')
    result = pd.concat([statistic, result], axis=0)
    return result

def do_output(result, new_file, single_file):
    
    warnings.simplefilter(action='ignore', category=FutureWarning)
    result = pd.DataFrame(result).T.sort_index()
    
    print("Available columns in result:", result.columns)
    print(result.head())  # Zeigt die ersten Zeilen des DataFrames an

    
    for column in result:
        result[column] = pd.to_numeric(result[column], errors='ignore')
    if single_file is False:
        result = do_statistic(result)
    else:
        result = result[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                         'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                         'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    result.to_csv(new_file, index=True)

def workflow(file_dir, new_file, write_csv, sample_size):
    path = files_manager.get_lib()
    temp_path = files_manager.tmp_dir(path, temp_path=None)
    paired = False
    single_file = False
    result = dict()
    buildbowtie2(path)

    if not os.path.exists(file_dir):
        raise ValueError(f"The provided path {file_dir} does not exist.")
    
    if glob.glob(f'{file_dir}**/*.fastq*', recursive=True) == [] and os.path.isdir(file_dir):
        files_manager.tmp_dir(path, temp_path)
        raise ValueError('There were no FASTQ files in this directory')
    if os.path.isfile(file_dir):
        single_file = True
        file_name = os.path.basename(file_dir)
        read2_file = os.path.join(os.path.dirname(file_dir), file_name.replace('_R1_', '_R2_'))
        if '_R1_' in file_name and os.path.exists(read2_file):
            paired = True
        # Sample reads from the forward file (R1) and use the same indices for the reverse file (R2)
        sampled_reads_R1, sampled_indices = sample_fastq(file_dir, sample_size)
        sampled_reads_R2, _ = sample_fastq(read2_file, sample_size, sampled_indices)
        
        # Save sampled reads for both R1 and R2
        temp_fastq_R1 = os.path.join(temp_path, 'sampled_R1.fastq')
        temp_fastq_R2 = os.path.join(temp_path, 'sampled_R2.fastq')
        save_sampled_fastq(sampled_reads_R1, temp_fastq_R1)
        save_sampled_fastq(sampled_reads_R2, temp_fastq_R2)

        aligned_path, Error = mapbowtie2(temp_fastq_R1, temp_fastq_R2, path, temp_path, paired)

        print(f"Mapping result: {aligned_path}, Error: {Error}")
        
        if Error is True:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file does not look like a fastq file')
        if paired is True and Output_counter.rawincount(f'{temp_path}paired.bed') == 0:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file has no Reads of the required mapping-quality')
        overlap(path, temp_path, aligned_path)
        if paired is False and Output_counter.rawincount(f'{temp_path}BED.bed') == 0:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file has no Reads of the required mapping-quality')
        file_name = file_name.rsplit('.f', 1)[0]
        file_name = file_name.replace('_R1_001', '')
        result[file_name] = Output_counter.create_row(temp_path, paired)
        
        print(f"Output for {file_name}: {result[file_name]}")
        
    elif os.path.isdir(file_dir):
        single_file = False
        total_files = len(glob.glob(f'{file_dir}**/*.fastq*', recursive=True))
        processed_files = 0
        for fq_file in glob.glob(f'{file_dir}**/*.fastq*', recursive=True):
            processed_files += 1
            paired = False
            if '_R2_' in fq_file:
                continue
            file_name = os.path.basename(fq_file)
            read2_file = os.path.join(os.path.dirname(fq_file), file_name.replace('_R1_', '_R2_'))
            if '_R1_' in file_name and os.path.exists(read2_file):
                paired = True

            # Sample reads from R1 and use the same indices for R2
            sampled_reads_R1, sampled_indices = sample_fastq(fq_file, sample_size)
            sampled_reads_R2, _ = sample_fastq(read2_file, sample_size, sampled_indices)
            
            temp_fastq_R1 = os.path.join(temp_path, 'sampled_R1.fastq')
            temp_fastq_R2 = os.path.join(temp_path, 'sampled_R2.fastq')
            save_sampled_fastq(sampled_reads_R1, temp_fastq_R1)
            save_sampled_fastq(sampled_reads_R2, temp_fastq_R2)

            aligned_path, Error = mapbowtie2(temp_fastq_R1, temp_fastq_R2, path, temp_path, paired)
            if Error is True:
                continue
            overlap(path, temp_path, aligned_path)

            print(f"Overlap completed for {aligned_path}")
            
            file_name = file_name.rsplit('.f', 1)[0]
            file_name = file_name.replace('_R1_001', '')
            result[file_name] = Output_counter.create_row(temp_path, paired)
            
            print(f"Output for {file_name}: {result[file_name]}")

    print("Result content before do_output:", result)
    
    files_manager.tmp_dir(path, temp_path)
    do_output(result, new_file, single_file)
    if write_csv is True:
        new_file = (f'{path}Output/{os.path.basename(os.path.dirname(file_dir))}.csv')
        do_output(result, new_file, single_file)
        

def main():
    parser = argparse.ArgumentParser(prog='VX detector', description=(
        'This programm tries to find which variable region of the 16S sequence was sequencend'))
    parser.add_argument('dir_path', help=('Directory path of the directory containing multiple fastq or fasta files.'))
    parser.add_argument('-o', '--output', dest='output_file', default=sys.stdout, 
                        help='User can specify a file format in which the output is written in the Output folder.')
    parser.add_argument('-c', '--csv', dest='write_csv', action='store_true', 
                        help='If set the output will be written in a .csv file in the Output folder')
    parser.add_argument('-s', '--sample_size', dest='sample_size', type=int, default=1000, 
                        help='Number of reads to sample from each FASTQ file')
    args = parser.parse_args()
    workflow(args.dir_path, args.output_file, args.write_csv, args.sample_size)

if __name__ == '__main__':
    main()
