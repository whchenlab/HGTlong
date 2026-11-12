import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from hgtlong.cli import filter_lengths
from hgtlong.reporting import write_report

class ReleaseUpdates(unittest.TestCase):
 def test_each_block_cutoff_and_default(self):
  events=[dict(event_id='short',a0=0,a1=499,b0=0,b1=500),dict(event_id='equal',a0=0,a1=500,b0=0,b1=500)]
  members=[dict(event_id=x['event_id']) for x in events]
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)
   self.assertEqual(len(filter_lengths(events,members,p,None)[0]),2)
   retained, evidence=filter_lengths(events,members,p,500)
   self.assertEqual([x['event_id'] for x in retained],['equal'])
   self.assertEqual(evidence,[dict(event_id='equal')])
   self.assertEqual(json.loads((p/'length_filter.json').read_text())['excluded'],1)
 def test_censored_reads_and_cli_cutoff(self):
  demo=Path(__file__).resolve().parents[1]/'demo'
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)
   def run(*args):
    r=subprocess.run([sys.executable,'-m','hgtlong',*map(str,args)],capture_output=True,text=True)
    self.assertEqual(r.returncode,0,r.stderr)
   run('extract','--sam',demo/'alignments.sam','--panel',demo/'panel.tsv','--out',p/'e.tsv','--workers',1)
   with open(p/'e.tsv') as f:rows=list(csv.DictReader(f,delimiter='\t'))
   singles=[r for r in rows if '_only' in r['read']]
   self.assertEqual(len(singles),4)
   for r in singles:
    side=r['resident'].lower()
    self.assertEqual(sum(int(r['anchor_'+side+k]) for k in ['0','1']),1)
   for mode in ['call','call-blocks']:
    extra=['--blocks',demo/'blocks.tsv'] if mode=='call-blocks' else []
    for limit,count in [(1000,1),(1001,0)]:
     out=p/(mode+str(limit))
     run(mode,'--evidence',p/'e.tsv','--outdir',out,'--min-hgt-length',limit,*extra)
     with open(out/'events.tsv') as f:self.assertEqual(len(list(csv.DictReader(f,delimiter='\t'))),count)
 def test_tooltip_single_source(self):
  with tempfile.TemporaryDirectory() as d:
   out=Path(d)/'report.html';write_report([],[],out)
   s=out.read_text()
   self.assertIn("return parent;}const x=document.createElementNS",s)
   self.assertIn(r'\nIdentity to A:',s)
