import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

class ReleaseDemo(unittest.TestCase):
    def test_installed_workflow_geometry_and_support(self):
        demo=Path(__file__).resolve().parents[1]/'demo'
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            def run(*args):
                p=subprocess.run([sys.executable,'-m','hgtlong',*map(str,args)],capture_output=True,text=True)
                self.assertEqual(p.returncode,0,p.stderr)
            run('extract','--sam',demo/'alignments.sam','--panel',demo/'panel.tsv','--out',out/'evidence.tsv','--workers','1')
            for mode in ['call','call-blocks']:
                extra=['--blocks',demo/'blocks.tsv'] if mode=='call-blocks' else []
                run(mode,'--evidence',out/'evidence.tsv','--outdir',out/mode,*extra)
                with open(out/mode/'events.tsv') as f: rows=list(csv.DictReader(f,delimiter='\t'))
                self.assertEqual(len(rows),1)
                self.assertEqual(rows[0]['status'],'read_supported')
                self.assertEqual([int(rows[0][k]) for k in ['a0','a1','b0','b1']], [1500,2500,3500,4500])
                self.assertEqual([int(rows[0]['n_'+k]) for k in ['a0','a1','b0','b1']], [4,4,4,4])
