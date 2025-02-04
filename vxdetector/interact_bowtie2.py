#!/usr/bin/python

import os
import shutil
import subprocess

# Check for the installation paths of Bowtie2, Samtools, and Bedtools
bowtie2_path = shutil.which('bowtie2')
if bowtie2_path is None:
    bowtie2_path = '$CONDA/bin/bowtie2'
samtools_path = shutil.which('samtools')
if samtools_path is None:
    samtools_path = '$CONDA/bin/samtools'
bedtools_path = shutil.which('bedtools')
if bedtools_path is None:
    bedtools_path = '$CONDA/bin/bedtools'


def buildbowtie2(path):
    """Builds Bowtie2 index for the reference genome."""
    input_ref = f'{path}Indexed_bt2/85_otus.fasta'  # Path to reference genome
    output_path = f'{path}Indexed_bt2/bowtie2'      # Output path for the index
    # Check if the index already exists; if not, build it
    if not os.path.exists(f'{output_path}.1.bt2'):
        os.system(
            f'{bowtie2_path}-build -f {input_ref}{output_path} > /dev/null')


def mapbowtie2(
        fasta_file, read2_file, path, temp_path, paired, bowtie2_params=None):
    """Maps reads against a previously built Bowtie2 index."""

    index_path = f'{path}Indexed_bt2/bowtie2'  # Path to the Bowtie2 index
    # Check if the index files exist
    if not os.path.exists(f'{index_path}.1.bt2'):
        raise FileNotFoundError(f'No index files found under "{index_path}"')

    log_path = f'{temp_path}bowtie2.log'       # Log file for Bowtie2
    bed_logpath = f'{temp_path}bed.log'         # Log file for Bedtools
    Error = False

    # Base command for Bowtie2
    bowtie2_base_cmd = [bowtie2_path, '-x', index_path, '--fast']

    # Add additional Bowtie2 parameters if provided
    if bowtie2_params:
        bowtie2_base_cmd += bowtie2_params.split()

    # Handle paired-end reads
    if paired:
        aligned_path = f'{temp_path}paired.bed'  # Output path for paired reads
        bowtie2_base_cmd += ['-1', fasta_file, '-2', read2_file]

        # Commands for Samtools and Bedtools
        cmd_sam = [samtools_path, 'view', '-b', '-q', '30', '-S', '-F', '4']
        cmd_bed = [bedtools_path, 'bamtobed', '-bedpe', '-i', 'stdin']

        # Run the Bowtie2, Samtools, and Bedtools pipeline
        with open(log_path, 'w') as log_file, \
             open(bed_logpath, 'w') as bed_log, \
             open(aligned_path, 'w') as output_file:
            bowtie_process = subprocess.Popen(
                bowtie2_base_cmd, stderr=log_file, stdout=subprocess.PIPE)
            samtools_process = subprocess.Popen(
                cmd_sam, stdin=bowtie_process.stdout, stdout=subprocess.PIPE)
            subprocess.run(cmd_bed, stdin=samtools_process.stdout,
                           stdout=output_file, stderr=bed_log)

    # Handle unpaired reads
    else:
        # Output path for unpaired reads
        aligned_path = f'{temp_path}unpaired.bed'
        bowtie2_base_cmd += ['-U', fasta_file]

        # Commands for Samtools and Bedtools
        cmd_sam = [samtools_path, 'view', '-b', '-q', '30', '-S', '-F', '4']
        cmd_bed = [bedtools_path, 'bamtobed', '-i', 'stdin']

        # Run the Bowtie2 and Samtools pipeline
        with open(log_path, 'w') as log_f, open(aligned_path, 'w') as output_f:
            bowtie_process = subprocess.Popen(
                bowtie2_base_cmd, stderr=log_f, stdout=subprocess.PIPE)
            samtools_process = subprocess.Popen(
                cmd_sam, stdin=bowtie_process.stdout, stdout=subprocess.PIPE)
            subprocess.run(
                cmd_bed, stdin=samtools_process.stdout, stdout=output_f)

    # Check for errors in the Bowtie2 log file
    with open(log_path, 'r') as log:
        if any("error" in line.lower() for line in log):
            Error = True
    # Return the path to the aligned file and error status
    return aligned_path, Error
