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
import random
import gzip


def sample_fastq(file_path, sample_size, sampled_indices=None, seed=None):
    #'''Get random parts from the FASTQ file based on shared indices for paired-end reads.'''
    """Read a (gzipped) fastQ file and randomly (without replacement) select
       a given number of reads.

       Assumes that each reads comes in exactly four lines. To save memory and
       have reasonable speed, we first pass through the file "just" to count
       the total amount of lines (div by 4 gives the number of reads).
       Then we randomly draw (without replacement) sample_size many numbers
       between 0 and line-numbers (div by 4 to only target header lines) and sort
       this list. Next, we read the file again and keep track of the line
       number (ln). Whenever the line number matches the current random line
       number, the file line and the three following are pushed to the
       sampled_readlines array.

       Parameters
       ----------
       file_path : str
            File path of the fastQ (gzipped) file containing the reads.
        sample_size : int
            Number of random reads taken from the fastQ file.
        sample_indices : ?
            ??? Ask Jeremy
        seed : int
            Random generator seed.

        Returns
        -------
        Tuple of list of random read lines and list of randomly selected line
        numbers.
     """

    def _count_generator(reader):
        b = reader(1024 * 1024)
        while b:
            yield b
            b = reader(1024 * 1024)

    def _get_filehandle(file_path):
        # open file either plain or as gzip compressed file
        FH = None
        if file_path.endswith('.gz'):
            FH = gzip.open(file_path, 'rb')
        else:
            FH = open(file_path, 'rb')
        return FH

    # first pass: count number of lines in the file
    FH = _get_filehandle(file_path)
    c_generator = _count_generator(FH.read)
    number_lines = sum(buffer.count(b'\n') for buffer in c_generator)
    assert number_lines % 4 == 0, "file %s does not contain a multiple of 4 lines! %i" % (file_path, number_lines)
    FH.close()

    # second pass: iterate through file with line counts
    # add next 4 lines to sampled_reads, iff linecount is in random selection
    sampled_readlines = []  # List to hold sampled reads
    if seed is not None:
        random.seed(seed)
    sel_reads = sorted(list(map(lambda x: x*4, random.sample(range(int(number_lines / 4) + 1), sample_size))))
    FH = _get_filehandle(file_path)
    ln = 0  # line number iterator
    j = 0  # iterator through sorted, random header line numbers
    while ln < number_lines:
        line = FH.readline()
        if ln == sel_reads[j]:
            sampled_readlines.append(str(line, encoding='utf-8'))
            for i in range(3):
                line = FH.readline()
                sampled_readlines.append(str(line, encoding='utf-8'))
                ln += 1
            if j + 1 < len(sel_reads):
                 j += 1
            else:
                break
        ln += 1
    FH.close()

    return sampled_readlines, sel_reads


def save_sampled_fastq(sampled_reads, output_path):
    '''Save sampled reads to a temporary FASTQ file.'''
    open_func = gzip.open if output_path.endswith('.gz') else open  # Use gzip if output is .gz

    with open_func(output_path, 'wt') as f:  # Write in text mode ('wt')
        f.writelines(sampled_reads)

def do_statistic(result):
    '''Compute statistics (mean and standard deviation) from alignment results.'''
    average = result.mean(numeric_only=True).to_frame().T  # Compute the mean of numeric columns
    region = (result['Sequenced variable region'].mode().values)  # Get most frequent region
    region = ' / '.join(str(r) for r in region)  # Join the modes into a single string
    region = region.replace('\'', '').replace('[', '').replace(']', '')  # Clean up string format
    average['Sequenced variable region'] = region  # Add region info to the average dataframe
    if 'Not properly paired' not in average.columns:
        average['Not properly paired'] = 'not paired'
    std_dev = result.std(numeric_only=True).to_frame().T  # Compute standard deviation of numeric columns
    statistic = pd.concat([average, std_dev], axis=0)  # Combine mean and standard deviation
    # Select specific columns for the final output
    statistic = statistic[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                           'V7', 'V8', 'V9', 'Not aligned to a variable region']]
    statistic['row_descriptor'] = ['Average', 'Standard deviation']  # Add row descriptors
    statistic = statistic.set_index('row_descriptor')  # Set descriptor as index
    result = pd.concat([statistic, result], axis=0)  # Append stats to original results
    return result

def do_output(result, new_file, single_file):
    '''Process and output results to a CSV file.'''
    warnings.simplefilter(action='ignore', category=FutureWarning)
    result = pd.DataFrame(result).T.sort_index()  # Convert result to DataFrame and sort

    # Convert columns to numeric where applicable
    for column in result:
        result[column] = pd.to_numeric(result[column], errors='ignore')

    if single_file is False:
        result = do_statistic(result)  # Add statistics if processing multiple files
    else:
        result = result[['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
                         'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
                         'V7', 'V8', 'V9', 'Not aligned to a variable region']]

    # Write the result DataFrame to a CSV file
    result.to_csv(new_file, index=True)

def workflow(file_dir, new_file, write_csv, sample_size, bowtie2_params):
    '''Main workflow for processing FASTQ files and performing sequence alignment.'''
    path = files_manager.get_lib()  # Get the main working directory
    temp_path = files_manager.tmp_dir(path, temp_path=None)  # Create a temporary directory
    paired = False  # Initially assume single-end
    single_file = False  # To check if processing a single file
    result = dict()  # Dictionary to store results from different files
    buildbowtie2(path)  # Build the Bowtie2 index

    # Check if the input path exists
    if not os.path.exists(file_dir):
        raise ValueError(f"The provided path {file_dir} does not exist.")

    # Check if there are FASTQ files in the directory
    if glob.glob(f'{file_dir}**/*.fastq*', recursive=True) == [] and os.path.isdir(file_dir):
        files_manager.tmp_dir(path, temp_path)
        raise ValueError('There were no FASTQ files in this directory')

    # Process single file
    if os.path.isfile(file_dir):
        single_file = True
        file_name = os.path.basename(file_dir)

        read2_file = os.path.join(os.path.dirname(file_dir), file_name.replace('_R1_', '_R2_'))  # Expected paired file

        # Check if paired-end data is available
        if '_R1_' in file_name and os.path.exists(read2_file):
            paired = True  # Mark as paired if both R1 and R2 files exist

        # Sample reads from R1 (and R2 if paired) with the same indices
        sampled_reads_R1, sampled_indices = sample_fastq(file_dir, sample_size)
        if paired:
            sampled_reads_R2, _ = sample_fastq(read2_file, sample_size, sampled_indices)
        else:
            sampled_reads_R2 = None  # No R2 file present, set R2 reads as None

        # Save sampled reads to temporary files
        temp_fastq_R1 = os.path.join(temp_path, 'sampled_R1.fastq')
        save_sampled_fastq(sampled_reads_R1, temp_fastq_R1)
        if paired:
            temp_fastq_R2 = os.path.join(temp_path, 'sampled_R2.fastq')
            save_sampled_fastq(sampled_reads_R2, temp_fastq_R2)

        # Run Bowtie2 alignment
        aligned_path, Error = mapbowtie2(temp_fastq_R1, temp_fastq_R2 if paired else None, path, temp_path, paired, bowtie2_params)

        if Error is True:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file does not look like a fastq file')
        if paired is True and Output_counter.rawincount(f'{temp_path}paired.bed') == 0:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file has no Reads of the required mapping-quality')

        # Perform overlap analysis using Bedtools
        overlap(path, temp_path, aligned_path)
        if not paired and Output_counter.rawincount(f'{temp_path}BED.bed') == 0:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file has no Reads of the required mapping-quality')

        # Generate results from output
        file_name = file_name.rsplit('.f', 1)[0]
        file_name = file_name.replace('_R1_001', '')
        result[file_name] = Output_counter.create_row(temp_path, paired)

    # Process directory of files
    elif os.path.isdir(file_dir):
        single_file = False
        total_files = len(glob.glob(f'{file_dir}**/*.fastq*', recursive=True))  # Count number of FASTQ files
        processed_files = 0
        for fq_file in glob.glob(f'{file_dir}**/*.fastq*', recursive=True):
            processed_files += 1
            paired = False
            if '_R2_' in fq_file:
                continue  # Skip R2 files to avoid processing duplicates

            file_name = os.path.basename(fq_file)
            read2_file = os.path.join(os.path.dirname(fq_file), file_name.replace('_R1_', '_R2_'))
            if '_R1_' in file_name and os.path.exists(read2_file):
                paired = True  # Mark as paired if both R1 and R2 files exist

            # Sample reads from R1 (and R2 if paired) with same indices
            sampled_reads_R1, sampled_indices = sample_fastq(fq_file, sample_size)
            sampled_reads_R2 = sample_fastq(read2_file, sample_size, sampled_indices)[0] if paired else None

            # Save sampled reads to temporary files
            temp_fastq_R1 = os.path.join(temp_path, 'sampled_R1.fastq')
            save_sampled_fastq(sampled_reads_R1, temp_fastq_R1)
            if paired:
                temp_fastq_R2 = os.path.join(temp_path, 'sampled_R2.fastq')
                save_sampled_fastq(sampled_reads_R2, temp_fastq_R2)

            # Run Bowtie2 alignment
            aligned_path, Error = mapbowtie2(temp_fastq_R1, temp_fastq_R2 if paired else None, path, temp_path, paired, bowtie2_params)
            if Error is True:
                continue  # Skip if there's an error in processing the file

            # Perform overlap analysis using Bedtools
            overlap(path, temp_path, aligned_path)

            # Generate results from output
            file_name = file_name.rsplit('.f', 1)[0]
            file_name = file_name.replace('_R1_001', '')
            result[file_name] = Output_counter.create_row(temp_path, paired)

    # Clean up temporary directory
    files_manager.tmp_dir(path, temp_path)

    # Output the result as CSV or display
    do_output(result, new_file, single_file)
    if write_csv is True:
        new_file = (f'{path}Output/{os.path.basename(os.path.dirname(file_dir))}.csv')
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
    parser.add_argument('-b','--bowtie2_params', dest='bowtie2_params', default='--threads 4, --fast',
                        help='Additional parameters to pass to Bowtie2')
    # Parse the user input arguments
    args = parser.parse_args()

    # Start the workflow with provided arguments
    workflow(args.dir_path, args.output_file, args.write_csv, args.sample_size, args.bowtie2_params)

if __name__ == '__main__':
    main()  # Entry point of the script
