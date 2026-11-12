"""Run the installed CLI against a synthetic alignment fixture and verify outputs."""
import argparse, csv, json, shutil, subprocess
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--outdir',type=Path,default=Path('demo-output'));p.add_argument('--min-hgt-length',type=int);args=p.parse_args()
here=Path(__file__).resolve().parent
exe=shutil.which('hgtlong')
if not exe: raise SystemExit('Install HGTlong in the active environment first.')
out=args.outdir.resolve();out.mkdir(parents=True,exist_ok=True)
def run(*args): subprocess.run([exe,*map(str,args)],check=True)
run('extract','--sam',here/'alignments.sam','--panel',here/'panel.tsv','--out',out/'evidence.tsv.gz','--workers','1')
extra=['--min-hgt-length',str(args.min_hgt_length)] if args.min_hgt_length is not None else []
run('call',*extra,'--evidence',out/'evidence.tsv.gz','--outdir',out/'read-first')
run('call-blocks',*extra,'--blocks',here/'blocks.tsv','--evidence',out/'evidence.tsv.gz','--outdir',out/'fixed-block')
truth=json.loads((here/'truth.json').read_text())
for method in ['read-first','fixed-block']:
 report=out/method/'report.html'
 report.write_text(report.read_text().replace('<h1>HGTlong</h1>','<h1>HGTlong · hgtlong_demo</h1>'))
 with open(out/method/'events.tsv') as f: rows=list(csv.DictReader(f,delimiter='\t'))
 if args.min_hgt_length is not None and args.min_hgt_length>1000:
  assert not rows
  print(method+': short demo block excluded as expected')
  continue
 assert len(rows)==1,(method,rows)
 for k in ['a0','a1','b0','b1','n_a0','n_a1','n_b0','n_b1']: assert int(rows[0][k])==truth[k],(method,k,rows[0][k])
 assert rows[0]['status']==truth['status']
 print(method+': expected geometry and four boundary counts verified; '+str(out/method/'report.html'))
