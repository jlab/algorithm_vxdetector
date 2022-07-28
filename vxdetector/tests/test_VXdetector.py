#!/usr/bin/python

import unittest
import tempfile
import io
import os
import sys


sys.path.append('$CONDA/lib/python3.9/site-packages')

import pandas as pd  # noqa: E402

sys.path.append(f'{__file__.rsplit("/", 2)[0]}/')

import VXdetector as vx  # noqa: E402
import shutil  # noqa: E402

path = f'{os.path.dirname(__file__)}/'
output_test = f'{path}test_data/Output_test.csv'
result = {'5011_S225_L001': {'Number of Reads': 64421,
                             'Unaligned Reads [%]': 12.189999999999998,
                             'Not properly paired': 0.005588239859673088,
                             'Sequenced variable region': 'V45', 'V1': 0.0,
                             'V2': 0.0, 'V3': 0.0, 'V4': 48.9871668789809,
                             'V5': 38.800461146496815,
                             'V6': 0.007457324840764332,
                             'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
                             'Not aligned to a variable region':
                             0.014914649681528664}}


class test_do_output(unittest.TestCase):
    def setUp(self):
        self.fp_tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.fp_tmpdir)

    def test_csv_output(self):
        new_file = f'{self.fp_tmpdir}test.csv'
        single_file = True
        vx.do_output(result, new_file, single_file)
        content = []
        with open(output_test)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(new_file)as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)

    def test_numeric_conversion(self):
        mixed_result = {'5011_S225_L001': {'Number of Reads': '64421',
                                           'Unaligned Reads [%]':
                                           '12.189999999999998',
                                           'Not properly paired':
                                           0.005588239859673088,
                                           'Sequenced variable region':
                                           'V45', 'V1': 0.0, 'V2': 0.0,
                                           'V3': 0.0, 'V4': '48.9871668789809',
                                           'V5': 38.800461146496815,
                                           'V6': 0.007457324840764332,
                                           'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
                                           'Not aligned to a variable region':
                                           0.014914649681528664}}
        single_file = True
        new_file = f'{self.fp_tmpdir}test3.csv'
        vx.do_output(mixed_result, new_file, single_file)
        content = []
        with open(output_test)as f:
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
                    '0.005588239859673088,V45,0.0,0.0,0.0,48.9871668789809,'
                    '38.800461146496815,0.007457324840764332,0.0,0.0,0.0'
                    ',0.014914649681528664\n')
        capturedOutput = io.StringIO()
        sys.stdout = capturedOutput
        new_file = sys.stdout
        single_file = True
        vx.do_output(result, new_file, single_file)
        sys.stdout = sys.__stdout__
        self.assertEqual(capturedOutput.getvalue(), expected)
        new_file = f'{self.fp_tmpdir}test2.csv'
        vx.do_output(result, new_file, single_file)
        self.assertTrue(os.path.exists(new_file))


class test_do_statistic(unittest.TestCase):
    def test_order(self):
        result = pd.read_csv(f'{path}test_data/mixed.csv', index_col=0)
        columns = ''.join(['Number of Reads',
                           'Unaligned Reads [%]', 'Not properly paired',
                           'Sequenced variable region', 'V1', 'V2', 'V3',
                           'V4', 'V5', 'V6', 'V7', 'V8', 'V9',
                           'Not aligned to a variable region'])
        statistic = ''.join(vx.do_statistic(result).columns)
        self.assertEqual(columns, statistic)

    def test_basic_function(self):
        expected = pd.read_csv(f'{path}test_data/expected.csv', index_col=0)
        result = pd.read_csv(f'{path}test_data/result.csv', index_col=0)
        statistic = vx.do_statistic(result).round(5)
        statistic.to_csv(f'{path}test_data/statistic1.csv')
        self.assertTrue(statistic.equals(expected))


class test_workflow(unittest.TestCase):
    def test_no_fastq(self):
        with self.assertRaises(ValueError) as cm:
            vx.workflow(f'{path}test_data/Indexed_bt2', sys.stdout, False)
        self.assertEqual('There were no FASTQ files in this directory',
                         str(cm.exception))

    def test_Error(self):
        with self.assertRaises(ValueError) as cm:
            vx.workflow(f'{path}test_data/Output_test.csv', sys.stdout, False)
        self.assertEqual('This file does not look like a fastq file',
                         str(cm.exception))


class test_total(unittest.TestCase):
    def setUp(self):
        self.fp_tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.fp_tmpdir)

    def test_singleFile(self):
        expected = f'{path}test_data/Output_test.csv'
        actual = f'{self.fp_tmpdir}/singleFile_test.csv'
        test_file = f'{path}test_data/5011_S225_L001_R1_001.fastq.gz'
        program_path = (f'{os.path.dirname(os.path.dirname(path))}/'
                        'VXdetector.py')
        os.system(f'python {program_path} {test_file} -o {actual}')
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
        expected = f'{path}test_data/dir_test.csv'
        actual = f'{self.fp_tmpdir}/singleFile_test.csv'
        test_file = f'{path}test_data/test_dir/'
        program_path = (f'{os.path.dirname(os.path.dirname(path))}/'
                        'VXdetector.py')
        os.system(f'python {program_path} {test_file} -o {actual}')
        content = []
        with open(expected)as f:
            for line in f:
                content.append(line.strip().split())
        output = []
        with open(actual)as f:
            for line in f:
                output.append(line.strip().split())
        self.assertEqual(output, content)


if __name__ == '__main__':
    unittest.main()
