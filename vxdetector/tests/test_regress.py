from unittest import TestCase, main
import tempfile
import shutil

from vxdetector.interact_bowtie2 import *
from vxdetector.VXdetector import workflow, output, do_statistic

class AlignmentTests(TestCase):
    def setUp(self):
        self.fp_dubiousV1 = "vxdetector/tests/data/87340_RF42TP0_S10_L001_R1_001.fastq"
        # construct bowtie2 index for Greengenes
        buildbowtie2("")

        self.fp_tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.fp_tmpdir)

    def test_correctV34(self):
        workflow("./", self.fp_tmpdir, self.fp_dubiousV1, None, "output.csv",
                 file_name='Your file', dir_name='', dir_path='',
                 mode='unpaired', read2_file='')
        do_statistic()
        # Todo: I don't understand how I can fill the global output DataFrame?!
        output.to_csv("help.csv", sep="\t")

if __name__ == '__main__':
    main()
