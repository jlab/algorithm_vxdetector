This Programm verifies which 16S variable region was sequenced.

Input:
    .fastq file containing 16S sequencing reads
Output:
    print()
    How many Reads were unaligned
    Which Region was most aligned and how high is the proability of this being the sequenced Region
    
    If no 16S Region was sequenced many reads will be unaligned

Requirements:
    python 3
    bowtie2
    samtools
    bedtools
    Module: os
    Module: subprocess
    Module: argparse
    Module: itertools
    Module: tempfile
    Module: shutil

To Do:
	-Error Message: query [...] is marked as paired, but its mate does not occur next to it in your BAM file
	-Multiple variable regions as a possibility for the most probable region
