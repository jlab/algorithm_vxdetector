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
    Module: tempfile
    Module: shutil
    
    To Do List:
        -Currently the used programms are refrenced by hardcoded paths
            -> not usable on PCs where these programms are saved somewhere else
        -Currently only the single most probable region is shown. 
            -> Modify in a way that reads spanning multiple regions are shown correctly
        -Programm filters out all fastq files that contain "_R2_"
            -> Is that the universal way of declaring reverse Reads?
        -Current Output is csv
        -Adding other functions such as primer verification or making the programm more efficiant
