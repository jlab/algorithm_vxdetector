#!/usr/bin/python

import unittest
import tempfile
import io
import os
import sys
from glob import glob
import pandas as pd
import vxdetector.VXdetector as vx
import shutil


class test_do_output(unittest.TestCase):
    def setUp(self):
        self.fp_tmpdir = tempfile.mkdtemp()
        self.path = f'{os.path.dirname(__file__)}/'
        self.output_test = f'{self.path}test_data/Output_test.csv'
        self.result = {'5011_S225_L001': {'Number of Reads': 64421,
                                          'Unaligned Reads [%]':
                                          12.189999999999998,
                                          'Not properly paired':
                                          0.004843141211716676,
                                          'Sequenced variable region': 'V45',
                                          'V1': 0.0, 'V2': 0.0, 'V3': 0.0,
                                          'V4': 45.375854271356786,
                                          'V5': 42.3973743718593,
                                          'V6': 0.0122571189279732,
                                          'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
                                          'Not aligned to a variable region':
                                          0.0245142378559464}}

    def tearDown(self):
        shutil.rmtree(self.fp_tmpdir)

    def test_csv_output(self):
        new_file = f'{self.fp_tmpdir}test.csv'
        single_file = True
        vx.do_output(self.result, new_file, single_file)
        content = []
        with open(self.output_test)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(new_file)as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)

    def test_numeric_conversion(self):
        mixed_result = {'5011_S225_L001': {'Not aligned to a variable region':
                                           0.0245142378559464,
                                           'Number of Reads': '64421',
                                           'Unaligned Reads [%]':
                                           '12.189999999999998',
                                           'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
                                           'Not properly paired':
                                           0.004843141211716676,
                                           'Sequenced variable region':
                                           'V45', 'V1': 0.0, 'V2': 0.0,
                                           'V3': 0.0, 'V4': 45.375854271356786,
                                           'V5': '42.3973743718593',
                                           'V6': 0.0122571189279732,
                                           }}
        single_file = True
        new_file = f'{self.fp_tmpdir}test3.csv'
        vx.do_output(mixed_result, new_file, single_file)
        content = []
        with open(self.output_test)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(new_file)as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)

    def test_output_options(self):
        expected = (',Number of Reads,Unaligned Reads [%],Not properly '
                    'paired,Sequenced variable region,V1,V2,V3,V4,V5,V6,'
                    'V7,V8,V9,Not aligned to a variable region\n'
                    '5011_S225_L001,64421,12.189999999999998,'
                    '0.004843141211716676,V45,0.0,0.0,0.0,45.375854271356786,'
                    '42.3973743718593,0.0122571189279732,0.0,0.0,0.0'
                    ',0.0245142378559464\n')
        capturedOutput = io.StringIO()
        sys.stdout = capturedOutput
        new_file = sys.stdout
        single_file = True
        vx.do_output(self.result, new_file, single_file)
        sys.stdout = sys.__stdout__
        self.assertEqual(capturedOutput.getvalue(), expected)
        new_file = f'{self.fp_tmpdir}test2.csv'
        vx.do_output(self.result, new_file, single_file)
        self.assertTrue(os.path.exists(new_file))


class test_do_statistic(unittest.TestCase):
    def setUp(self):
        self.fp_tmpdir = tempfile.mkdtemp()
        self.path = f'{os.path.dirname(__file__)}/'

    def tearDown(self):
        shutil.rmtree(self.fp_tmpdir)

    def test_order(self):
        result = pd.read_csv(f'{self.path}test_data/mixed.csv', index_col=0)
        columns = ''.join(['Number of Reads',
                           'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3',
                           'V4', 'V5', 'V6', 'V7', 'V8', 'V9',
                           'Not aligned to a variable region'])
        statistic = ''.join(vx.do_statistic(result).columns)
        self.assertEqual(columns, statistic)

    def test_basic_function(self):
        expected = pd.read_csv(f'{self.path}test_data/expected.csv',
                               index_col=0)
        result = pd.read_csv(f'{self.path}test_data/result.csv', index_col=0)
        statistic = vx.do_statistic(result).round(5)
        self.assertTrue(statistic.equals(expected))
        expected = pd.read_csv(f'{self.path}test_data/expected_unpaired.csv',
                               index_col=0)
        result = pd.read_csv(f'{self.path}test_data/result_unpaired.csv',
                             index_col=0)
        statistic = vx.do_statistic(result).round(5)
        statistic.to_csv(f'{self.fp_tmpdir}/statistic1.csv')
        self.assertTrue(statistic.equals(expected))


class test_workflow(unittest.TestCase):
    def setUp(self):
        self.fp_tmpdir = tempfile.mkdtemp()
        self.path = f'{os.path.dirname(__file__)}/'

    def tearDown(self):
        if os.path.exists(f'{__file__.rsplit("/", 3)[0]}/'
                          'Output/test_data.csv'):
            os.remove(f'{__file__.rsplit("/", 3)[0]}/Output/test_data.csv')
        file_list = glob(f'{__file__.rsplit("/", 3)[0]}/Indexed_bt2/*.bt2')
        directory_list = glob(f'{__file__.rsplit("/", 3)[0]}/tmp_files_*')
        if os.path.exists(f'{self.path}test_data/dir_test_actual.csv'):
            file_list.append(f'{self.path}test_data/dir_test_actual.csv')
        for file in file_list:
            os.remove(file)
        shutil.rmtree(self.fp_tmpdir)
        for directory in directory_list:
            shutil.rmtree(directory)

    def test_singleFile(self):
        expected = f'{self.path}test_data/Output_test.csv'
        actual = f'{self.fp_tmpdir}/singleFile_test.csv'
        test_file = f'{self.path}test_data/5011_S225_L001_R1_001.fastq.gz'
        vx.workflow(test_file, actual, False)
        content = []
        with open(expected)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(actual)as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)

    def test_directory(self):
        expected = f'{self.path}test_data/dir_test.csv'
        actual = f'{self.path}/test_data/dir_test_actual.csv'
        test_file = f'{self.path}test_data/test_dir/'
        vx.workflow(test_file, actual, False)
        content = []
        with open(expected)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(actual)as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)

    def test_c_option(self):
        expected = f'{self.path}test_data/Output_test.csv'
        actual = sys.stdout
        test_file = f'{self.path}test_data/5011_S225_L001_R1_001.fastq.gz'
        vx.workflow(test_file, actual, True)
        content = []
        with open(expected)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(f'{__file__.rsplit("/", 3)[0]}/Output/test_data.csv')as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)

    def test_no_fastq(self):
        with self.assertRaises(ValueError) as cm:
            vx.workflow(f'{self.path}test_data/Indexed_bt2', sys.stdout, False)
        self.assertEqual('There were no FASTQ files in this directory',
                         str(cm.exception))

    def test_Error(self):
        with self.assertRaises(ValueError) as cm:
            vx.workflow(f'{self.path}test_data/Output_test.csv',
                        sys.stdout, False)
        self.assertEqual('This file does not look like a fastq file',
                         str(cm.exception))

    def test_raise(self):
        with self.assertRaises(ValueError) as cm:
            vx.workflow(f'{self.path}test_data/test_dir/no_qual_test.fastq',
                        sys.stdout, False)
        self.assertEqual('This file has no Reads of the required '
                         'mapping-quality', str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            vx.workflow(f'{self.path}test_data/test_dir/no_qual_'
                        'paired_R1_001.fastq', sys.stdout, False)
        self.assertEqual('This file has no Reads of the required '
                         'mapping-quality', str(cm.exception))
