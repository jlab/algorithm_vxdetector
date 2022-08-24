#!/usr/bin/python

import unittest
import tempfile
import os
import vxdetector.interact_bedtools as ibe
import shutil

path = f'{__file__.rsplit("/", 3)[0]}/'
test_paired = f'{os.path.dirname(__file__)}/test_data/paired/BED.bed'
test_paired_no = f'{os.path.dirname(__file__)}/test_data/paired/noOver.bed'


class test_no_overlap(unittest.TestCase):
    def setUp(self):
        self.fp_tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.fp_tmpdir)

    def test_basic_function(self):
        temp_path = self.fp_tmpdir
        aligned_path = (f'{os.path.dirname(__file__)}/test_data/'
                        'paired/paired.bed')
        ibe.no_overlap(path, temp_path, aligned_path)
        content = []
        with open(test_paired_no)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(f'{temp_path}noOver.bed')as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)


class test_overlap(unittest.TestCase):
    def setUp(self):
        self.fp_tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.fp_tmpdir)

    def test_basic_function(self):
        temp_path = self.fp_tmpdir
        aligned_path = (f'{os.path.dirname(__file__)}/test_data/'
                        'paired/paired.bed')
        ibe.overlap(path, temp_path, aligned_path)
        content = []
        with open(test_paired)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(f'{temp_path}BED.bed')as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)
