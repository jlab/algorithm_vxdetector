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
    """
    Converts a BIOM file into a FASTQ file.

    This function reads a BIOM file, which is a binary JSON-like format often used in metagenomics to store 
    sequence counts and metadata. It extracts sequence identifiers and sequences, and formats them into a 
    FASTQ format. Each sequence is given a default quality score to make it compatible with downstream 
    bioinformatics tools that require FASTQ files.

    Parameters
    ----------
    biom_file : str
        Path to the input BIOM file containing sequence data and associated metadata.
    fastq_output : str
        Path to the output FASTQ file where the converted sequences will be saved. Supports gzip compression 
        if the filename ends with '.gz'.

    Returns
    -------
    None
        Writes the converted sequences to the specified FASTQ file.
    
    Notes
    -----
    - Ensure that the input BIOM file is properly formatted and includes sequence data.
    - The FASTQ format generated assumes a fixed quality score, which may need adjustment based on downstream 
      applications.
    """
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


def sample_fastq(file_path, sample_size, sampled_indices=None, seed=None):
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
        sampled_indices : list of int, optional
            Consistent sampling across paired end files
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
    """
    Writes sampled sequencing reads into a FASTQ file.

    The function accepts a list of sampled reads in FASTQ format and saves them to the specified output 
    file. If the output file ends with '.gz', the file is automatically compressed using gzip. This function 
    ensures that the sampled reads are stored efficiently for use in subsequent steps like alignment or 
    analysis.

    Parameters
    ----------
    sampled_reads : list of str
        List of reads in FASTQ format, where each read consists of four lines:
        - Line 1: Sequence identifier (e.g., @SEQ_ID).
        - Line 2: Raw sequence.
        - Line 3: Optional '+' separator.
        - Line 4: Quality scores.
    output_path : str
        Path to the output file. If the file ends with '.gz', it will be written in a compressed format.

    Returns
    -------
    None
        Saves the sampled reads to the output file.
    """
    open_func = gzip.open if output_path.endswith('.gz') else open  # Use gzip if output is .gz

    with open_func(output_path, 'wt') as f:  # Write in text mode ('wt')
        f.writelines(sampled_reads)


def do_statistic(result):
    """
    Performs statistical analysis on sequencing data.

    This function analyzes the results from processing FASTQ files, calculating metrics such as the mean and 
    standard deviation for numeric columns in the data. It also identifies the most common variable region 
    across sequences, appending these statistics as new columns in the DataFrame.

    Parameters
    ----------
    result : pandas.DataFrame
        DataFrame containing processed sequencing data. Expected columns include numeric metrics and 
        categorical data such as variable regions.

    Returns
    -------
    pandas.DataFrame
        Updated DataFrame with additional columns for calculated statistics.
    """
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
    """
    Writes processing results to a CSV file.

    Converts sequencing results stored in a dictionary into a DataFrame, optionally calculates statistics, 
    and writes the data to a CSV file for downstream analysis or record-keeping. If the results are from a 
    single input file, no additional statistics are calculated.

    Parameters
    ----------
    result : dict
        Dictionary containing processing results, keyed by file name.
    new_file : str
        Path to the output CSV file where results will be saved.
    single_file : bool
        If True, skips statistical calculations as the results are from a single input file.

    Returns
    -------
    None
        Saves the results to the specified file.
    """
    
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
    Processes an individual FASTQ file through sampling, alignment, and overlap analysis.

    This function handles all steps for processing a single FASTQ file. It includes optional random sampling 
    of reads to reduce data size, alignment to a reference genome using Bowtie2, and analysis of overlap 
    regions using Bedtools.

    Parameters
    ----------
    fq_file : str
        Path to the input FASTQ file to process.
    path : str
        Path to the Bowtie2 library and index files.
    bowtie2_params : str
        Additional parameters for Bowtie2 alignment (e.g., "--fast --threads 4").
    sample_size : int, optional
        Number of reads to randomly sample. If None, the entire file is processed.

    Returns
    -------
    dict or None
        A dictionary with processing results, including file name and metrics. Returns None if an error 
        occurs during processing.
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
    Executes the full workflow for processing sequencing data.

    The workflow includes directory traversal, conversion of BIOM files to FASTQ format (if needed), and 
    parallelized processing of FASTQ files. Results are aggregated and optionally written to a CSV file.

    Parameters
    ----------
    file_dir : str
        Path to the directory or file containing input sequencing data.
    new_file : str
        Path to the CSV file where results will be saved.
    write_csv : bool
        If True, writes the results to a CSV file.
    sample_size : int
        Number of reads to sample from each FASTQ file. If None, all reads are processed.
    bowtie2_params : str
        Additional parameters for Bowtie2 alignment.

    Returns
    -------
    None
        Executes the workflow and outputs results as specified.
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
    """
    Parses user input and initiates the workflow.

    This function uses argparse to handle command-line arguments, passing them to the workflow function. 
    It provides flexibility for customizing the input directory, Bowtie2 parameters, and output options.

    Parameters
    ----------
    None

    Returns
    -------
    None
        Executes the script based on user-provided inputs.
    
    Notes
    -----
    - Command-line usage should follow the format specified in the argparse help message.
    - Ensure all required dependencies (Bowtie2, Bedtools) are installed and accessible.
    """
    
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