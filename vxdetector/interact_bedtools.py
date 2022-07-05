#!/usr/bin/python

import os

def overlap(path, temp_path):
	bedtools_filepath = '/usr/bin/bedtools'
	BAM_filepath = temp_path + 'bam_sorted.bam'
	S_ref = path + 'Indexed_bt2/annoted_ref.bed' 
	#reference "genome" Created by using a aligned greengenes databank and deleting all "-" while keeping track where the variable Regions are 
	BED_filepath = temp_path + 'BED.bed'
	cmd = bedtools_filepath + ' intersect -f 0.5 -wb -a ' + BAM_filepath + ' -b ' + S_ref + ' -bed > ' + BED_filepath #Intersects the .bam file with the reference
	os.system(cmd)
