#!/usr/bin/python

import os
import tempfile
import shutil


def tmp_dir(path, temp_path):
    r'''Creates and deletes a temporary folder

    This function creates a temporary folder if none exists
    and deletes the temporary folder if it does exist.

    Parameters
    ----------
    path : str
        The function needs to know where the temporary folder
        should be created. In the context of this program the path
        will always be the path of the program directory.
    temp_path : str or None
        Path of the temporary folder. If None a new temporary directory
        will be created otherwise the directory specified here will be deleted.

    Returns
    -------
    temp_path : str
        Path to the created temporary folder.

    '''
    if temp_path is None or os.path.exists(temp_path) is False:
        temp_path = f"{tempfile.mkdtemp(prefix='tmp_files_', dir=path)}/"
        return temp_path
    else:
        shutil.rmtree(temp_path)


def get_lib(program_path=None):
    r'''Checks necessary files

    This function checks if all necessery files are where they should be.
    It also returns the program path.
    Should a necessary file be missing it needs to be replaced.
    annoted_ref.bed can be newly generated with the code provided
    in vxdetector/Indexed_bt2/code_for_reference/

    Parameter
    ---------
    program_path : None or str
        Enables the user to manually call the function on a specific path.
        Mainly used for testing.

    Returns
    -------
    program_path : str
        Path to the directory in which the programm files are stored.

    '''
    if program_path is None:
        program_path = os.path.dirname(__file__)
        program_path = f"{os.path.dirname(program_path)}/"
    # finds the parent directory of the folder in which the
    # program is saved
    if os.path.exists(f'{program_path}Output/'):
        pass
    else:
        os.mkdir(f'{program_path}Output/')
    # creates the standard output folder
    if os.path.exists(f'{program_path}Indexed_bt2/'):
        if os.path.exists(f'{program_path}Indexed_bt2/annoted_ref.bed'):
            pass
        else:
            raise FileNotFoundError('ERROR: It seems the 16S variable '
                                    'region boundary reference file is '
                                    'missing.')
        if os.path.exists(f'{program_path}Indexed_bt2/85_otus.fasta'):
            pass
        else:
            raise FileNotFoundError('ERROR: It seems the Greengenes '
                                    '"85_otus.fasta" file is missing.')
        if os.path.exists(f'{program_path}Indexed_bt2/85_otus_aligned.fasta'):
            pass
        else:
            raise FileNotFoundError('ERROR: It seems the Greengenes '
                                    '"85_otus_aligned.fasta" file is missing.')
    else:
        raise FileNotFoundError('ERROR: It seems the "Indexed_bt2" '
                                'directory is missing.')

    return program_path
