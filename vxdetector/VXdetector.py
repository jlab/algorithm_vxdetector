#!/usr/bin/python

import argparse
import glob
import os
import sys
import time

# SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# sys.path.append(os.path.dirname(SCRIPT_DIR))
# sys.path.append('$CONDA/lib/python3.9/site-packages')

import pandas as pd  # noqa: E402
from vxdetector.interact_bowtie2 import mapbowtie2, buildbowtie2  # noqa: E402
from vxdetector.interact_bedtools import overlap  # noqa: E402
import vxdetector.Output_counter as Output_counter  # noqa: E402
import vxdetector.files_manager as files_manager  # noqa: E402


def do_statistic(result):
    r'''Statistical analysis of given directory

    This function calculates the average (mean) value and
    the standard deviation of all dataFrame collumns containg
    numeric data.
    It also declares which variable region(s) is/are the most
    probable sequenced region over all fastq files.

    Parameters
    ----------
    result : pandas dataFrame
        DataFrame with filenames as indey and the obtained information
        in the columns. Each file is a single row.

    Returns
    -------
    result : pandas dataFrame
        DataFrame is the same as above.
        The first two lines were added which show the average (mean) and
        standard deviation of each column containg numeric data.
        Index of these rows is ['Average', 'Standard deviation']
        Columns are:
        ['Number of Reads', 'Unaligned Reads [%]', 'Not properly paired',
         'Sequenced variable region', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6',
         'V7', 'V8', 'V9', 'Not aligned to a variable region']

    '''
    average = result.mean(numeric_only=True).to_frame().T
    # calculates mean value for every numeric column
    region = (result['Sequenced variable region'].mode().values)
    # finds the most common value in column 'Sequenced variable region'
    region = ' / '.join(str(r) for r in region)
    region = region.replace('\'', '')
    region = region.replace('[', '')
    region = region.replace(']', '')
    # formats the mode() output for easy viewing
    average['Sequenced variable region'] = region
    if 'Not properly paired' not in average.columns:
        average['Not properly paired'] = 'not paired'
        # adds 'Not properly paired' column if reads are unpaired
    std_dev = result.std(numeric_only=True).to_frame().T
    # calculates standard deviation for every numeric column
    statistic = pd.concat([average, std_dev], axis=0)
    statistic = statistic[['Number of Reads',
                           'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3',
                           'V4', 'V5', 'V6', 'V7', 'V8', 'V9',
                           'Not aligned to a variable region']]
    # combines the average and std_dev dataframes and sorts the columns
    # in the correct order
    statistic['row_descriptor'] = ['Average', 'Standard deviation']
    statistic = statistic.set_index('row_descriptor')
    # adds a descriptive index
    result = pd.concat([statistic, result], axis=0)
    return result


def do_output(result, new_file, single_file):
    r'''Writes Output

    This function writes the obtained information in a csv format.

    Parameters
    ---------
    result : dict
        Dictionary with filenames as keys and the obtained information
        as a Dictionary stored in the value.
    new_file : str or _io.TextIOWrapper
        Filepath to the output file.
        If none is given to the program (_io.TextIOWrapper) output
        will go to STDOUT.
    single_file : Bool
        Signifies wether a directory or a single file was given via terminal
        input.

    '''
    result = pd.DataFrame(result).T.sort_index()
    # converts dict to dataframe and sorts the dataframe by index
    # allows for easy navigation within the file
    for column in result:
        result[column] = pd.to_numeric(result[column], errors='ignore')
        # converts all cell values to a numeric data type if applicable
    if single_file is False:
        result = do_statistic(result)
    else:
        result = result[['Number of Reads',
                         'Unaligned Reads [%]', 'Not properly paired',
                         'Sequenced variable region', 'V1', 'V2', 'V3',
                         'V4', 'V5', 'V6', 'V7', 'V8', 'V9',
                         'Not aligned to a variable region']]
    result.to_csv(new_file, index=True)


def workflow(file_dir, new_file, write_csv):
    r'''Worker function

    This function is the center piece of this program.
    It does some preliminary work such as grabbing all fastq files
    in the given directory.
    After that it calls other functions which analyse the found
    fastq-files.

    Parameters
    ----------
    file_dir : str
        Given filepath. Can be a path to single file which leads
        to visual output whichcan easily be understood or a
        path to a directory which produces csv output.
    new_file : str
        Filepath to the Output file. If the flag is not set
        output will be printed via STDOUT.
    write_csv : Bool
        Wether or not a csv file should be written in the
        standard Output folder of this program.

    '''
    path = files_manager.get_lib()
    temp_path = files_manager.tmp_dir(path, temp_path=None)
    # sets the paths of the programm itself and a temporary folder
    paired = False
    result = dict()
    buildbowtie2(path)
    # builds bowtie2 index
    if glob.glob(f'{file_dir}**/*.fastq*', recursive=True) == [] \
       and os.path.isdir(file_dir):
        files_manager.tmp_dir(path, temp_path)
        raise ValueError('There were no FASTQ files '
                         'in this directory')
        # checks if given directory contains fastq files
    if os.path.isfile(file_dir):
        single_file = True
        file_name = os.path.basename(file_dir)
        read2_file = os.path.join(os.path.dirname(file_dir),
                                  file_name.replace('_R1_', '_R2_'))
        if '_R1_' in file_name and os.path.exists(read2_file):
            paired = True
        # searches for a reverse read file
        aligned_path, Error = mapbowtie2(file_dir, read2_file,
                                         path, temp_path, paired)
        # The Programm bowtie2 is used to align the Reads to a reference
        # 16S database.
        if Error is True:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file does not look like a fastq file')
            # raises error if given file is not a fastq file
        if paired is True and Output_counter.rawincount(f'{temp_path}'
                                                        'paired.bed') == 0:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file has no Reads of the required '
                             'mapping-quality')
        overlap(path, temp_path, aligned_path)
        # look which reads intersect with which variable Region
        if paired is False and Output_counter.rawincount(f'{temp_path}'
                                                         'BED.bed') == 0:
            files_manager.tmp_dir(path, temp_path)
            raise ValueError('This file has no Reads of the required '
                             'mapping-quality')
        file_name = file_name.rsplit('.f', 1)[0]
        file_name = file_name.replace('_R1_001', '')
        # reformats filename for easy viewing
        result[file_name] = Output_counter.create_row(temp_path, paired)
        # streamlines generated output to a visual terminal output
    elif os.path.isdir(file_dir):
        single_file = False
        for fq_file in glob.glob(f'{file_dir}**/*.fastq*', recursive=True):
            paired = False
            if '_R2_' in fq_file:
                continue
            file_name = os.path.basename(fq_file)
            read2_file = os.path.join(os.path.dirname(fq_file),
                                      file_name.replace('_R1_', '_R2_'))
            if '_R1_' in file_name and os.path.exists(read2_file):
                paired = True
            # searches for a reverse read file
            aligned_path, Error = mapbowtie2(fq_file, read2_file,
                                             path, temp_path, paired)
            # The Programm bowtie2 is used to align the Reads to a reference
            # 16S database.
            if Error is True:
                continue
            overlap(path, temp_path, aligned_path)
            # look which reads intersect with which variable Region
            file_name = file_name.rsplit('.f', 1)[0]
            file_name = file_name.replace('_R1_001', '')
            # reformats filename for easy viewing
            result[file_name] = Output_counter.create_row(temp_path,
                                                          paired)
            # streamlines generated output to a visual terminal output
    files_manager.tmp_dir(path, temp_path)
    # deletes temporary folder
    do_output(result, new_file, single_file)
    # writes ouput eiher to STDOUT or to a file specified
    # via the -o option
    if write_csv is True:
        new_file = (f'{path}Output/'
                    f'{os.path.basename(os.path.dirname(file_dir))}.csv')
        do_output(result, new_file, single_file)
        # writes csv file to the standard output folder


def main():
    r'''Main function

    This function uses argparse to interpret user
    input and then calls the workflow function which
    does the actual work.

    Iput otions
        -o, --output :
            User specified output path
        -c, --csv :
            If set a csv file will be written into
            the standard Output folder

    '''
    parser = argparse.ArgumentParser(prog='VX detector', description=(
        'This programm tries to find which variable region of the 16S '
        'sequence was sequencend'))
    parser.add_argument('dir_path',
                        help=('Directory path of the directory containing '
                              'multiple fastq or fasta files.'))
    parser.add_argument('-o', '--output', dest='output_file', default=sys.stdout, help='User can specify \
                        a file format in which the output is written in the \
                        Output folder.')
    parser.add_argument('-c', '--csv', dest='write_csv', action='store_true',
                        help='If set the output will be written in a .csv file \
                        in the Output folder')
    args = parser.parse_args()
    # allows terminal input
    workflow(args.dir_path, args.output_file, args.write_csv)


if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
    with open('/homes/jgroos/Desktop/Output/time.txt', 'a') as f:
        f.write("\n --- %s seconds ---" % (time.time() - start_time))
