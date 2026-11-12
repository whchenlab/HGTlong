import tempfile
from pathlib import Path
import unittest
from hgtlong.qc import annotate_masks


class QC(unittest.TestCase):
    def test_union_and_unassessed(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'mask.bed';p.write_text('a\t10\t30\trRNA\na\t20\t40\trRNA\n')
            events=[dict(ctgA='a',ctgB='b',a0=0,a1=100,b0=100,b1=200),
                    dict(ctgA='c',ctgB='b',a0=0,a1=100,b0=100,b1=200)]
            out=list(annotate_masks(events,p,{'a','b'}))
            self.assertEqual(out[0]['mask_bp_A'],30)
            self.assertEqual(out[0]['qc_status'],'repeat_associated')
            self.assertEqual(out[0]['mask_classes'],'rRNA')
            self.assertEqual(out[1]['repeat_screen'],'incomplete_annotation')
