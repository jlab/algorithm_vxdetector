#!/usr/bin/python

import unittest
import os
import vxdetector.Output_counter as oc


class test_rawincount(unittest.TestCase):
    def setUp(self):
        self.data_path = f'{os.path.dirname(__file__)}/test_data/'

    def test_counting(self):
        count = oc.rawincount(f'{self.data_path}'
                              '5011_S225_L001_R2_001.fastq.gz')
        self.assertEqual(count, 39324)

    def test_raise(self):
        filename = 'notafile'
        with self.assertRaises(FileNotFoundError) as cm:
            oc.rawincount(filename)
        self.assertEqual(f'It seems {filename} is missing. '
                         'Check if you deleted the temp_dir',
                         str(cm.exception))


class test_region_count(unittest.TestCase):
    def setUp(self):
        self.data_path = f'{os.path.dirname(__file__)}/test_data/'

    def test_new_row(self):
        test_data = f'{self.data_path}unpaired/'
        row = {'Number of Reads': 64421,
               'Unaligned Reads [%]': 4.019999999999996,
               'Not properly paired': 'not paired',
               'Sequenced variable region': 'V45', 'V1': 0.0, 'V2': 0.0,
               'V3': 0.0, 'V4': 50.265919045461416, 'V5': 45.714080954538595,
               'V6': 0.0, 'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
               'Not aligned to a variable region': 0.0}
        paired = False
        new_row = {'Number of Reads': 64421,
                   'Unaligned Reads [%]': 4.019999999999996}
        regions = {'V1': 0, 'V2': 0, 'V3': 0, 'V4': 6935, 'V5': 6307,
                   'V6': 0, 'V7': 0, 'V8': 0, 'V9': 0}
        new_row = oc.region_count(test_data, paired, new_row, regions)
        self.assertEqual(new_row, row)
        test_data = f'{self.data_path}paired/'
        row = {'Number of Reads': 64421,
               'Unaligned Reads [%]': 4.019999999999996,
               'Not properly paired': 0.005588239859673088,
               'Sequenced variable region': 'V45', 'V1': 0.0, 'V2': 0.0,
               'V3': 0.0, 'V4': 49.23987456606844, 'V5': 44.780950092025044,
               'V6': 0.0, 'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
               'Not aligned to a variable region': 1.9591753419065112}
        paired = True
        new_row = oc.region_count(test_data, paired, new_row, regions)
        self.assertEqual(new_row, row)

    def test_aligned_to_no_v(self):
        test_data = f'{self.data_path}Indexed_bt2/'
        row = {'Number of Reads': 64421,
               'Unaligned Reads [%]': 4.019999999999996,
               'Not properly paired': 'not paired',
               'Sequenced variable region': 'No variable Region',
               'V1': 0.0, 'V2': 0.0,
               'V3': 0.0, 'V4': 0.0, 'V5': 0.0,
               'V6': 0.0, 'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
               'Not aligned to a variable region': 95.98}
        paired = False
        new_row = {'Number of Reads': 64421,
                   'Unaligned Reads [%]': 4.019999999999996}
        regions = {'V1': 0, 'V2': 0, 'V3': 0, 'V4': 0, 'V5': 0,
                   'V6': 0, 'V7': 0, 'V8': 0, 'V9': 0}
        new_row = oc.region_count(test_data, paired, new_row, regions)
        self.assertEqual(new_row, row)

    def test_no_aligned(self):
        test_data = f'{self.data_path}unpaired/'
        row = {'Number of Reads': 64421,
               'Unaligned Reads [%]': 100,
               'Not properly paired': 'not paired',
               'Sequenced variable region': 'Not 16S',
               'V1': 0.0, 'V2': 0.0,
               'V3': 0.0, 'V4': 0.0, 'V5': 0.0,
               'V6': 0.0, 'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
               'Not aligned to a variable region': 0.0}
        paired = False
        new_row = {'Number of Reads': 64421,
                   'Unaligned Reads [%]': 100}
        regions = {'V1': 0, 'V2': 0, 'V3': 0, 'V4': 0, 'V5': 0,
                   'V6': 0, 'V7': 0, 'V8': 0, 'V9': 0}
        new_row = oc.region_count(test_data, paired, new_row, regions)
        self.assertEqual(new_row, row)


class test_create_row(unittest.TestCase):
    def setUp(self):
        self.data_path = f'{os.path.dirname(__file__)}/test_data/'

    def test_return_value(self):
        test_data = f'{self.data_path}unpaired/'
        row = {'Not aligned to a variable region': 0.0,
               'Not properly paired': 'not paired', 'Number of Reads': 64421,
               'Sequenced variable region': 'V45',
               'Unaligned Reads [%]': 4.019999999999996, 'V1': 0.0,
               'V2': 0.0, 'V3': 0.0, 'V4': 50.265919045461416,
               'V5': 45.714080954538595, 'V6': 0.0, 'V7': 0.0,
               'V8': 0.0, 'V9': 0.0}
        new_row = oc.create_row(test_data, paired=False)
        self.assertEqual(new_row, row)
        test_data = f'{self.data_path}paired/'
        row = {'Number of Reads': 64421,
               'Unaligned Reads [%]': 12.189999999999998,
               'Not properly paired': 'not paired',
               'Sequenced variable region': 'V45', 'V1': 0.0,
               'V2': 0.0, 'V3': 0.0, 'V4': 45.375854271356786,
               'V5': 42.3973743718593, 'V6': 0.0122571189279732,
               'V7': 0.0, 'V8': 0.0, 'V9': 0.0,
               'Not aligned to a variable region': 0.0245142378559464}

        new_row = oc.create_row(test_data, paired=False)
        self.assertEqual(new_row, row)

    def test_raise(self):
        with self.assertRaises(FileNotFoundError) as cm:
            oc.create_row(self.data_path, paired=False)
        self.assertEqual(f'It seems {self.data_path}BED.bed is missing. '
                         'Check if you deleted the temp_dir',
                         str(cm.exception))
        no_log_path = f'{self.data_path}Indexed_bt2/'
        with self.assertRaises(FileNotFoundError) as cm:
            oc.create_row(no_log_path, paired=False)
        self.assertEqual(f'It seems {no_log_path}bowtie2.log is missing. '
                         'Check if you deleted the temp_dir',
                         str(cm.exception))
