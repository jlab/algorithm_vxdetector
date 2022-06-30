#!/usr/bin/python

import os

def no_overlap(path, temp_path, aligned_path):
	bedtools_filepath = '/usr/bin/bedtools'
	S_ref = f'{path}Indexed_bt2/annoted_ref.bed'  
	noOver_filepath = f'{temp_path}noOver.bed'
	os.system(f'{bedtools_filepath} intersect -v -f 0.5 -a {aligned_path} -b {S_ref} > {noOver_filepath}')


def overlap(path, temp_path, aligned_path):
	bedtools_filepath = '/usr/bin/bedtools'
	S_ref = f'{path}Indexed_bt2/annoted_ref.bed' 
	#reference "genome" Created by using a aligned greengenes databank and deleting all "-" while keeping track where the variable Regions are 
	BED_filepath = f'{temp_path}BED.bed' 
	os.system(f'{bedtools_filepath} intersect -f 0.5 -a {S_ref} -b {aligned_path} > {BED_filepath}') #Intersects the .bam file with the reference
	no_overlap(path, temp_path, aligned_path)


