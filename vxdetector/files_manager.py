#!/usr/bin/python

import os
import tempfile
import shutil


def tmp_dir(path, temp_path):
    # creates a temporary folder
    # if called again it will delete the temporary folder
    if os.path.exists(temp_path):
        shutil.rmtree(temp_path)
    else:
        temp_path = f"{tempfile.mkdtemp(prefix='tmp_files_', dir=path)}/"

    return temp_path


def get_lib():
    # looks for the required files (e.g. a greengenes reference file)
    # and reports back if/what file(s) are missing
    programm_path = f'{os.path.dirname(__file__)}'
    programm_path = f"{programm_path.rsplit('/', 1)[0]}/"
    if os.path.exists(f'{programm_path}Output/'):
        pass
    else:
        os.mkdir(f'{programm_path}Output/')
    if os.path.exists(f'{programm_path}Indexed_bt2/'):
        if os.path.exists(f'{programm_path}Indexed_bt2/annoted_ref.bed'):
            pass
        else:
            print('ERROR: It seems the 16S variable region boundary'
                  ' reference file is missing.')
        if os.path.exists(f'{programm_path}Indexed_bt2/85_otus.fasta'):
            pass
        else:
            print('ERROR: It seems the Greengenes "85_otus.fasta" file is'
                  ' missing.')
        if os.path.exists(f'{programm_path}Indexed_bt2/85_otus_aligned.fasta'):
            pass
        else:
            print('ERROR: It seems the Greengenes "85_otus_aligned.fasta" '
                  'file is missing.')
    else:
        print('ERROR: It seems the "Indexed_bt2" directory is missing.')

    return programm_path
