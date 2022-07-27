#!/usr/bin/python

import os
import shutil

bowtie2_path = shutil.which('bowtie2')
if bowtie2_path is None:
    bowtie2_path = '$CONDA/bin/bowtie2'
samtools_path = shutil.which('samtools')
if samtools_path is None:
    samtools_path = '$CONDA/bin/samtools'
bedtools_path = shutil.which('bedtools')
if bedtools_path is None:
    bedtools_path = '$CONDA/bin/bedtools'
# finds where bowtie2, samtools and bedtools
# are installed


def buildbowtie2(path):
    r'''Builds bowtie2 index

    This function builds an index for bowtie2 it is the equivalent
    of building it manually with the "" command.

    Parameters
    ----------
    path : str
        The program filepath. The directory where it needs to look
        for the reference genome and where the index should be saved.

    '''
    input_ref = f'{path}Indexed_bt2/85_otus.fasta'
    # reference genome greengenes is used
    output_path = f'{path}Indexed_bt2/bowtie2'
    if os.path.exists(f'{output_path}.1.bt2'):
        # only builds an index for the alignment if there isn't already one
        pass
    else:
        os.system(f'{bowtie2_path}-build -f {input_ref} {output_path} '
                  '> /dev/null')
        # if a indexfile is missing a new index is build


def mapbowtie2(fasta_file, read2_file, path, temp_path, paired):
    r'''Maps reads against index

    This function maps read files against a previously build index.
    The bowtie standard output .sam is piped into samtools view and
    converted to .bam. In this step all unmapped reads and those with a
    lower mapquality than requested are discarded.
    If the reads are paired the .bam file is converted to .bed with bedtools
    bamtobed that way they are properly paired and can be intersected.

    Parameters
    ----------
    fasta_file : str
        Path to a .fastq file (or .fastq.gz). This file contains the (forward)
        reads which are mapped against the reference.
    read2_file : str
        Path to a .fastq file (or .fastq.gz) containing backwards reads.
        Is disregarded if paired = False. If not needed can be ''.
    path : str
        The program filepath. The directory where it needs to look
        for the index files.
    temp_path : str
        Path to a (temporary) folder where the resulting .log files and output
        files are saved
    paired : Bool
        Dictates wether bowtie2 alignes paired or unpaired reads.
        If paired is set to False bowtie2 will do an unpaired alignment
        otherwise it will do a paired one.

    Returns
    -------
    aligned_path : str
        Path to the converted bowtie2 Output file in which either paired or
        unpaired reads have been mapped against the reference.
    Error : Bool
        Wether or not bowtie2 ended with an error.
        In the context of the whole program Error = True leads to either
        a raised Exception (if only a single file was given) or to the current
        file being skipped (if a directory was given).

    '''
    index_path = f'{path}Indexed_bt2/bowtie2'
    if os.path.exists(f'{index_path}.1.bt2') is False:
        raise Exception(f'No Index files found under "{index_path}"')
        # raises an Exception if the Index files cannot be found
    log_path = f'{temp_path}bowtie2.log'
    bed_logpath = f'{temp_path}bed.log'
    # declares various filepaths
    Error = False
    if paired is True:
        aligned_path = f'{temp_path}paired.bed'
        cmd = f'{bowtie2_path} -x {index_path} -1 {fasta_file} -2 {read2_file} \
                --fast 2> {log_path} | {samtools_path} view -b -q 25 -S -F 4 \
                | {bedtools_path} bamtobed -bedpe -i stdin > {aligned_path} \
                2> {bed_logpath}'
        # Should a backward read be found both files will be given to bowtie2.
        # After converting .sam to .bam a conversion to .bed is done to
        # properly mate the pairs
    else:
        aligned_path = f'{temp_path}unpaired.bam'
        cmd = f'{bowtie2_path} -x {index_path} -q -U {fasta_file} \
                 --fast 2> {log_path} | {samtools_path} view -b -q 25 -S -F 4 \
                 -o {aligned_path}'
        # Should no backward read be found it will just use the forward
        # read and does an alignment followed by a pipe to convert
        # the bowtie2 output .sam to a .bam file
    os.system(cmd)
    with open(log_path, 'r') as log:
        lines = log.readlines()
        try:
            int(lines[0].split()[0])
        except ValueError:
            Error = True
        # Checks if bowtie2 exited with an error

    return aligned_path, Error
