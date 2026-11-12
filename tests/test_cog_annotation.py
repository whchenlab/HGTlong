import csv,json,subprocess,sys,tempfile,unittest
from pathlib import Path
class COGAnnotation(unittest.TestCase):
 def test_consensus_conflict_and_footprint(self):
  script=Path(__file__).resolve().parents[1]/'scripts/annotate_cog_v2.py'
  if not script.exists():self.skipTest('Source utility not included in installed wheel')
  with tempfile.TemporaryDirectory() as d:
   p=Path(d);(p/'cog-24.def.tab').write_text('COG1\tX\tTransposase\ttnp\nCOG2\tJ\tRibosomal protein\trpl\nCOG3\tG\tHydrolase\thyd\nCOG4\tL\tSite-specific recombinase XerC\txerC\n')
   (p/'cog-24.cog.csv').write_text('g,a,P1,100,1-100,100,COG1,COG1,0\ng,a,P2,100,1-100,100,COG2,COG2,0\ng,a,P3,100,1-100,100,COG3,COG3,0\ng,a,P4,200,160-200,41,COG1,COG1,1\ng,a,P5,100,1-100,100,COG4,COG4,0\n')
   (p/'genes.tsv').write_text('gene_id\tmag\tcontig\tstart\tend\tstrand\nq1\tA\ta\t0\t300\t+\nq2\tA\ta\t300\t600\t+\nq3\tA\ta\t600\t900\t+\nq4\tA\ta\t900\t1200\t+\n')
   (p/'hits.tsv').write_text('q1\tP1\t90\t100\t1\t100\t1\t100\t1e-30\t100\t100\t100\nq2\tP2\t90\t100\t1\t100\t1\t100\t1e-30\t100\t100\t100\nq2\tP3\t90\t100\t1\t100\t1\t100\t1e-30\t99\t100\t100\nq3\tP4\t90\t100\t1\t100\t1\t140\t1e-30\t100\t100\t200\nq4\tP5\t90\t100\t1\t100\t1\t100\t1e-30\t100\t100\t100\n')
   result=subprocess.run([sys.executable,str(script),'--genes',str(p/'genes.tsv'),'--hits',str(p/'hits.tsv'),'--db',str(p),'--out',str(p/'out.tsv')],capture_output=True,text=True)
   self.assertEqual(result.returncode,0,result.stderr)
   with open(p/'out.tsv') as f:rows=list(csv.DictReader(f,delimiter='\t'))
   self.assertEqual(rows[0]['cog_categories'],'X');self.assertEqual(rows[0]['mobile_types'],'transposition')
   self.assertEqual(rows[1]['annotation_status'],'conflicting_near_best_COGs');self.assertEqual(rows[1]['cog_categories'],'')
   self.assertEqual(rows[2]['annotation_status'],'no_qualified_COG_hit')

   self.assertEqual(rows[3]['cog_categories'],'L');self.assertEqual(rows[3]['mobile_types'],'')
