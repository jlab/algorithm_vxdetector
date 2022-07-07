#!/usr/bin/python

import os
import shutil

bedtools_path = shutil.which('bedtools')


def no_overlap(path, temp_path, aligned_path):
    S_ref = f'{path}Indexed_bt2/annoted_ref.bed'
    noOver_filepath = f'{temp_path}noOver.bed'
    os.system(f'{bedtools_path} intersect -v -f 0.5 -a \
                {aligned_path} -b {S_ref} > {noOver_filepath}')
    # This function uses the intersect -v option to find all reads not
    # aligned to any variable region but still mapped to the 16S greengenes
    # reference.


def overlap(path, temp_path, aligned_path):
    S_ref = f'{path}Indexed_bt2/annoted_ref.bed'
    # reference "genome" Created by using a aligned greengenes databank and
    # deleting all "-" while keeping track where the variable Regions are
    BED_filepath = f'{temp_path}BED.bed'
    # Intersects the aligned bowtie2 Output file with the reference
    os.system(f'{bedtools_path} intersect -f 0.5 -a {S_ref} -b \
                {aligned_path} > {BED_filepath}')
    no_overlap(path, temp_path, aligned_path)
