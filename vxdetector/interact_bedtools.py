#!/usr/bin/python

import os
import shutil

bedtools_path = shutil.which('bedtools')
if bedtools_path is None:
    bedtools_path = '$CONDA/bin/bedtools'
# finds where bedtools is installed


def no_overlap(path, temp_path, aligned_path):
    r'''Creates a no overlap file

    This function creates a file containing all reads that
    do not overlap with any variable region.
    It uses the intersect -v option to find all reads not
    aligned to any variable region but still mapped to the 16S greengenes
    reference.

    Parameters
    ----------
    path : str
        Program path and the directory where all it needs to
        look for the annoted_ref.bed file.
    temp_path : str
        Path to a temporary directory where the no overlap file
        is saved to.
    aligned_path : str
        Filepath of the converted Bowtie2 output file.
        This can either be a .bam or .bed file containing unpaired and
        paired reads respectivly.

    '''
    S_ref = f'{path}Indexed_bt2/annoted_ref.bed'
    noOver_filepath = f'{temp_path}noOver.bed'
    os.system(f'{bedtools_path} intersect -v -f 0.5 -a \
                {aligned_path} -b {S_ref} -bed > {noOver_filepath}')


def overlap(path, temp_path, aligned_path):
    r'''Creates an overlap file

    This function creates a file containing all reads that
    do overlap with a variable region.

    Parameters
    ----------
    path : str
        Program path and the directory where all it needs to
        look for the annoted_ref.bed file.
    temp_path : str
        Path to a temporary directory where the overlap file
        is saved to.
    aligned_path : str
        Filepath of the converted Bowtie2 output file.
        This can either be a .bam or .bed file containing unpaired and
        paired reads respectivly.

    '''
    S_ref = f'{path}Indexed_bt2/annoted_ref.bed'
    # reference "genome" Created by using a aligned greengenes databank and
    # deleting all "-" while keeping track where the variable Regions are
    BED_filepath = f'{temp_path}BED.bed'
    os.system(f'{bedtools_path} intersect -f 0.5 -a {S_ref} -b \
                {aligned_path} > {BED_filepath}')
    # Intersects the aligned bowtie2 Output file with the reference
    no_overlap(path, temp_path, aligned_path)
