#!/usr/bin/python

import unittest
import os
import tempfile
import shutil
import sys

import_path = f'{__file__.rsplit("/", 2)[0]}/'
sys.path.append(import_path)

import files_manager as fm  # noqa: E402


class test_tmp_dir(unittest.TestCase):
    def test_basic_function(self):
        path = os.path.dirname(__file__)
        temp_path = fm.tmp_dir(path, temp_path='')
        self.assertTrue(os.path.exists(temp_path))
        fm.tmp_dir(path, temp_path)
        self.assertFalse(os.path.exists(temp_path))


class test_get_lib(unittest.TestCase):
    def setUp(self):
        self.fp_tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.fp_tmpdir)

    def test_program_path(self):
        program_path = fm.get_lib()
        self.assertTrue(os.path.exists(program_path))
        self.assertTrue(os.path.exists(f'{program_path}/vxdetector/'))
        self.assertTrue(os.path.exists(f'{program_path}/vxdetector/'
                                       'VXdetector.py'))

    def test_raise_exceptions_get_lib(self):
        with self.assertRaises(FileNotFoundError) as cm:
            fm.get_lib(program_path=self.fp_tmpdir)
        self.assertEqual('ERROR: It seems the "Indexed_bt2" '
                         'directory is missing.', str(cm.exception))
        os.mkdir(f'{self.fp_tmpdir}Indexed_bt2/')
        with self.assertRaises(FileNotFoundError) as cm:
            fm.get_lib(program_path=self.fp_tmpdir)
        self.assertEqual('ERROR: It seems the 16S variable region '
                         'boundary reference file is missing.',
                         str(cm.exception))
        open(f'{self.fp_tmpdir}Indexed_bt2/annoted_ref.bed', 'w').close()
        with self.assertRaises(FileNotFoundError) as cm:
            fm.get_lib(program_path=self.fp_tmpdir)
        self.assertEqual('ERROR: It seems the Greengenes "85_otus.fasta" '
                         'file is missing.', str(cm.exception))
        open(f'{self.fp_tmpdir}Indexed_bt2/85_otus.fasta', 'w').close()
        with self.assertRaises(FileNotFoundError) as cm:
            fm.get_lib(program_path=self.fp_tmpdir)
        self.assertEqual('ERROR: It seems the Greengenes '
                         '"85_otus_aligned.fasta" file is missing.',
                         str(cm.exception))
        open(f'{self.fp_tmpdir}Indexed_bt2/85_otus_aligned.fasta', 'w').close()
        fm.get_lib(program_path=self.fp_tmpdir)
        self.assertTrue(os.path.exists(f'{self.fp_tmpdir}Output/'))
