import unittest
from hgtlong.sequence_qc import metrics,reverse_complement


class SequenceQC(unittest.TestCase):
    def test_partial_alignment_not_full_identity(self):
        row=metrics('A\t1000\t0\t500\t+\tB\t1000\t0\t500\t499\t500\t60')
        self.assertEqual(row['block_identity'],.998)
        self.assertEqual(row['sequence_qc'],'partial_local_alignment')

    def test_orientation_conflict(self):
        row=metrics('A\t1000\t0\t1000\t-\tB\t1000\t0\t1000\t999\t1000\t60')
        self.assertEqual(row['sequence_qc'],'orientation_conflict')

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement('ACGTNRY'),'RYNACGT')

    def test_no_alignment(self):
        self.assertEqual(metrics('')['sequence_qc'],'no_local_alignment')
