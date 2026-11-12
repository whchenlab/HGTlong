import tempfile
from pathlib import Path
import unittest
from hgtlong.panel import build_panel


class Panel(unittest.TestCase):
    def test_namespacing(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'in').mkdir()
            (root/'in/A.fa').write_text('>same\nACGT\n')
            (root/'in/B.fa').write_text('>same\nNNNN\n')
            result=build_panel(root/'in',root/'out')
            self.assertEqual(result,dict(mags=2,contigs=2,bases=8))
            self.assertIn('>B__same',(root/'out/panel.fa').read_text())
            with self.assertRaisesRegex(ValueError,'already exist'):build_panel(root/'in',root/'out')

    def test_duplicate_contig_is_atomic(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);(root/'in').mkdir()
            (root/'in/A.fa').write_text('>same\nACGT\n>same\nAAAA\n')
            with self.assertRaisesRegex(ValueError,'Duplicate'):build_panel(root/'in',root/'out')
            self.assertFalse((root/'out/panel.fa').exists())
