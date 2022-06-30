This Programm verifies which 16S variable region was sequenced.

Input:
    .fastq file containing 16S sequencing reads
    or a Directory containing (or subfolders containing) fasta or fastq files
Output:
    print(Probabilities for every region)
    CSV file containing read_designation, Percentage unaligned reads, Probabilities of every region
  
    

Requirements:
    python 3
    bowtie2
    samtools
    bedtools
    Module: os
    Module: subprocess
    Module: argparse
    Module: tempfile
    Module: shutil
    Module: itertools
    Module: csv

To Do:
	-Currently the used programms are refrenced by hardcoded paths
            -> not usable on PCs where these programms are saved somewhere else
        -Adding other functions such as primer verification or making the programm more efficiant
