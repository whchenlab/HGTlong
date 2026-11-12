import unittest
from hgtlong.screening import screen_event


class Screening(unittest.TestCase):
    def event(self,**updates):
        row=dict(status='read_supported',n_a0=3,n_a1=3,n_b0=3,n_b1=3,
                 qc_status='annotation_screen_pass',sequence_qc='both_blocks_covered',
                 block_identity=.995,aligned_fraction_A=.99,aligned_fraction_B=.99)
        row.update(updates);return row

    def test_pass(self):
        self.assertEqual(screen_event(self.event())['screened_candidate'],1)

    def test_lower_identity_not_promoted(self):
        row=screen_event(self.event(block_identity=.973))
        self.assertEqual(row['screened_candidate'],0)
        self.assertIn('block_identity',row['screen_fail_reasons'])

    def test_repeat_unknown_or_low_coverage_fails(self):
        for changes in [dict(qc_status='repeat_associated'),dict(qc_status='not_fully_assessed'),
                        dict(aligned_fraction_A=.5),dict(block_identity=''),dict(n_b1=2)]:
            self.assertEqual(screen_event(self.event(**changes))['screened_candidate'],0)
