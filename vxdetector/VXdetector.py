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
import tempfile
from multiprocessing import Pool
import biom
import shutil
import gzip
import random

def biom_to_fastq(biom_file, fastq_output):
    '''Converts a BIOM file to FASTQ format by extracting sequences and ensuring they are compatible with FASTQ standards.'''
    with open(fastq_output, "w") as f:
        table = biom.load_table(biom_file)  # Load the BIOM table
        sequences = table.ids(axis='observation')  # Extract sequences
        for seq_id in sequences:
            # Truncate the sequence ID if it’s too long for Bowtie2's requirements
            truncated_id = seq_id[:50]  # Cut the sequence ID to 50 characters (adjust if needed)
            # Check and sanitize the sequence to avoid unwanted characters
            sanitized_seq = ''.join(filter(str.isalpha, seq_id))  # Ensure only letters are in the sequence
            f.write(f"@{truncated_id}\n{sanitized_seq}\n+\n{'I' * len(sanitized_seq)}\n")  # Write in FASTQ format
    print(f"BIOM file successfully converted to FASTQ format at: {fastq_output}")


def sample_fastq(file_path, sample_size, sampled_indices=None):
    '''Get random parts from the FASTQ file based on shared indices for paired-end reads.'''
    sampled_reads = []  # List to hold sampled reads
    open_func = gzip.open if file_path.endswith('.gz') else open

    # Count total reads in the file by iterating line by line
    with open_func(file_path, 'rt') as f:
        total_reads = sum(1 for _ in f) // 4

    # Adjust sample_size if it's greater than total_reads
    if sample_size > total_reads:
        sample_size = total_reads

    if sampled_indices is None:
        sampled_indices = sorted(random.sample(range(total_reads), sample_size))

    with open_func(file_path, 'rt') as f:
        read_idx = 0
        read_buffer = []
        for line in f:
            read_buffer.append(line)
            if len(read_buffer) == 4:  # Once we have a full read
                if read_idx in sampled_indices:
                    sampled_reads.extend(read_buffer)  # Add the read if it matches the sampled index
                read_buffer = []  # Clear the buffer for the next read
                read_idx += 1  # Move to the next read

                # Stop if we've collected enough reads
                if len(sampled_reads) >= sample_size * 4:
                    break

    return sampled_reads, sampled_indices


def save_sampled_fastq(sampled_reads, output_path):
    '''Save sampled reads to a temporary FASTQ file.'''
    open_func = gzip.open if output_path.endswith('.gz') else open  # Use gzip if output is .gz

    with open_func(output_path, 'wt') as f:  # Write in text mode ('wt')
        f.writelines(sampled_reads)


def do_statistic(result):
    '''Performs statistical analysis on the DataFrame

    Calculates the mean and standard deviation for all numeric columns
    in the DataFrame, and determines the most common sequenced variable 
    region. Adds these statistics to the DataFrame and returns the updated 
    DataFrame.
    '''
    average = result.mean(numeric_only=True).to_frame().T  # Calculate mean for numeric columns
    region = (result['Sequenced variable region'].mode().values)  # Find the most common variable region
    region = ' / '.join(str(r) for r in region)  # Format the result for better readability
    region = region.replace('\'', '').replace('[', '').replace(']', '')  # Clean up formatting
    average['Sequenced variable region'] = region
    if 'Not properly paired' not in average.columns:
        average['Not properly paired'] = 'not paired'  # Handle cases where reads are not paired
    std_dev = result.std(numeric_only=True).to_frame().T  # Calculate standard deviation for numeric columns
    statistic = pd.concat([average, std_dev], axis=0)  # Combine mean and standard deviation dataframes
    # Select relevant columns and reorder them
    statistic = statistic[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                           'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    statistic['row_descriptor'] = ['Average', 'Standard deviation']  # Add descriptors for the statistics
    statistic = statistic.set_index('row_descriptor')  # Set the row descriptors as the index
    result = pd.concat([statistic, result], axis=0)  # Combine the statistical results with the original data
    return result


def do_output(result, new_file, single_file):
    '''Writes the results into a CSV file

    Converts the dictionary of results into a DataFrame, calculates statistics 
    if working with multiple files, and writes the data into a CSV file.
    '''
    warnings.simplefilter(action='ignore', category=FutureWarning)  # Ignore warnings about future behavior
    result = pd.DataFrame(result).T.sort_index()  # Convert the dictionary to a DataFrame and sort by index
    for column in result:
        result[column] = pd.to_numeric(result[column], errors='ignore')  # Convert columns to numeric where possible
    if single_file is False:
        result = do_statistic(result)  # If multiple files, calculate statistics
    else:
        # Select relevant columns for single file processing
        result = result[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                         'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                         'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    result.to_csv(new_file, index=True)  # Write the results to a CSV file


def process_file(fq_file, path, bowtie2_params, sample_size=None):
    """
    Processes a single FASTQ file, including optional read sampling, alignment with Bowtie2, 
    and overlap analysis using Bedtools.
    """
    paired = False
    result = {}
    try:
        with tempfile.TemporaryDirectory() as temp_path:
            file_name = os.path.basename(fq_file)
            read2_file = os.path.join(os.path.dirname(fq_file), file_name.replace('_R1_', '_R2_'))
            if '_R1_' in file_name and os.path.exists(read2_file):
                paired = True  # Paired-end sequencing detected

            if sample_size:  # Perform sampling if required
                sampled_reads_R1, sampled_indices = sample_fastq(fq_file, sample_size)
                temp_fastq_R1 = os.path.join(temp_path, 'sampled_R1.fastq')
                save_sampled_fastq(sampled_reads_R1, temp_fastq_R1)
                fq_file = temp_fastq_R1

                if paired:  # Process the reverse reads for paired-end data
                    sampled_reads_R2, _ = sample_fastq(read2_file, sample_size, sampled_indices)
                    temp_fastq_R2 = os.path.join(temp_path, 'sampled_R2.fastq')
                    save_sampled_fastq(sampled_reads_R2, temp_fastq_R2)
                    read2_file = temp_fastq_R2

            # Perform alignment and overlap analysis
            aligned_path, Error = mapbowtie2(fq_file, read2_file, path, temp_path, paired, bowtie2_params)
            if Error:
                return None
            overlap(path, temp_path, aligned_path)
            file_name = file_name.rsplit('.f', 1)[0].replace('_R1_001', '')
            result[file_name] = Output_counter.create_row(temp_path, paired)
    except Exception as e:
        print(f"Error processing file {fq_file}: {e}")
    return result


def workflow(file_dir, new_file, write_csv, sample_size, bowtie2_params):
    """
    Manages the overall processing workflow, including directory traversal, BIOM file conversion, 
    and multiprocessing for FASTQ file processing.
    """
    path = files_manager.get_lib()
    buildbowtie2(path)  # Prepare Bowtie2 library
    result = {}
    single_file = False

    if file_dir.endswith('.biom'):  # Handle BIOM file input
        temp_fastq = os.path.join(tempfile.gettempdir(), 'temp_sequences.fastq')
        biom_to_fastq(file_dir, temp_fastq)
        result = process_file(temp_fastq, path, bowtie2_params, sample_size)
        single_file = True
    elif os.path.isfile(file_dir):  # Process single FASTQ file
        single_file = True
        result = process_file(file_dir, path, bowtie2_params, sample_size)
    elif os.path.isdir(file_dir):  # Process directory of FASTQ files
        fastq_files = glob.glob(f'{file_dir}/**/*.fastq*', recursive=True)
        with Pool() as pool:
            results = pool.starmap(process_file, [(fq_file, path, bowtie2_params, sample_size)
                                                  for fq_file in fastq_files if '_R2_' not in fq_file])
        for res in results:
            if res:
                result.update(res)
    do_output(result, new_file, single_file)


def main():
    '''Main function to parse user arguments and initiate the workflow.'''
    parser = argparse.ArgumentParser(prog='VX detector', description=(
        'This program tries to find which variable region of the 16S sequence was sequenced'))
    parser.add_argument('dir_path', help=('Directory path of the directory containing multiple fastq or fasta files.'))
    parser.add_argument('-o', '--output', dest='output_file', default=sys.stdout, 
                        help='User can specify a file format in which the output is written in the Output folder.')
    parser.add_argument('-c', '--csv', dest='write_csv', action='store_true', 
                        help='If set the output will be written in a .csv file in the Output folder')
    parser.add_argument('-s', '--sample_size', dest='sample_size', type=int, default=1000, 
                        help='Number of reads to sample from each FASTQ file')
    parser.add_argument('-b','--bowtie2_params', dest='bowtie2_params', default= None,
                        help='Additional parameters to pass to Bowtie2')
    # Parse the user input arguments
    args = parser.parse_args()
    
    # Start the workflow with provided arguments
    workflow(args.dir_path, args.output_file, args.write_csv, args.sample_size, args.bowtie2_params)

if __name__ == '__main__':
    main()  # Entry point of the script