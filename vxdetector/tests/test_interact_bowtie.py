#!/usr/bin/python

import unittest
import tempfile
import os
import sys
from glob import glob

import_path = f'{__file__.rsplit("/", 2)[0]}/'
sys.path.append(import_path)

import interact_bowtie2 as ibo  # noqa: E402
import shutil  # noqa: E402

import_path = f'{__file__.rsplit("/", 3)[0]}/'
sys.path.append(import_path)
path = f'{os.path.dirname(__file__)}/'
fasta_file = f'{path}test_data/5011_S225_L001_R1_001.fastq.gz'
read2_file = f'{path}test_data/5011_S225_L001_R2_001.fastq.gz'
test_unpaired = f'{path}test_data/unpaired/unpaired.sam'
test_paired = f'{path}test_data/paired/paired.bed'
samtools_path = shutil.which('samtools')


class test_mapbowtie2(unittest.TestCase):
    def setUp(self):
        ibo.buildbowtie2(f'{path}test_data/')
        self.fp_tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        file_list = glob(f'{path}test_data/Indexed_bt2/*.bt2')
        for file in file_list:
            os.remove(file)
        shutil.rmtree(self.fp_tmpdir)

    def test_aligned_path(self):
        path = f'{os.path.dirname(__file__)}/test_data/'
        paired = False
        temp_path = self.fp_tmpdir
        aligned_path, Error = ibo.mapbowtie2(fasta_file, read2_file,
                                             path, temp_path, paired)
        self.assertEqual(aligned_path, f'{temp_path}unpaired.bam')
        self.assertEqual(Error, False)
        paired = True
        aligned_path, Error = ibo.mapbowtie2(fasta_file, read2_file,
                                             path, temp_path, paired)
        self.assertEqual(aligned_path, f'{temp_path}paired.bed')
        self.assertEqual(Error, False)

    def test_pipe(self):
        path = f'{os.path.dirname(__file__)}/test_data/'
        paired = False
        temp_path = self.fp_tmpdir
        content = []
        with open(test_unpaired)as f:
            for line in f:
                content.append(line.strip().split())
        ibo.mapbowtie2(fasta_file, read2_file, path, temp_path, paired)
        os.system(f'{samtools_path} view {temp_path}unpaired.bam '
                  f'> {temp_path}unpaired.sam')
        output = []
        with open(f'{temp_path}unpaired.sam')as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)
        paired = True
        content = []
        with open(test_paired)as f:
            for line in f:
                content.append(line.strip().split())
        ibo.mapbowtie2(fasta_file, read2_file, path, temp_path, paired)
        output = []
        with open(f'{temp_path}paired.bed')as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)

    def test_wrong_file_type(self):
        path = f'{os.path.dirname(__file__)}/test_data/'
        read2_file = f'{path}test_data/paired/BED.bed'
        paired = True
        aligned_path, Error = ibo.mapbowtie2(fasta_file, read2_file,
                                             path, self.fp_tmpdir, paired)
        self.assertEqual(Error, True)
        paired = False
        fasta_file_local = f'{path}test_data/paired/BED.bed'
        aligned_path, Error = ibo.mapbowtie2(fasta_file_local, read2_file,
                                             path, self.fp_tmpdir, paired)
        self.assertEqual(Error, True)


class test_buildbowtie2(unittest.TestCase):
    def tearDown(self):
        path = f'{os.path.dirname(__file__)}/'
        file_list = glob(f'{path}test_data/Indexed_bt2/*.bt2')
        for file in file_list:
            os.remove(file)

    def test_index(self):
        path = f'{os.path.dirname(__file__)}/'
        ibo.buildbowtie2(f'{path}test_data/')
        self.assertEqual(os.path.exists(f'{path}test_data/Indexed_bt2'
                                        '/bowtie2.1.bt2'), True)


if __name__ == '__main__':
    unittest.main()