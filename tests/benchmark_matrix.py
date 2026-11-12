"""Independent engineered-sequence sensitivity matrix, not an empirical ONT benchmark.

Fixed before inspecting the matrix: errors 0/0.5/1/2%, nominal depths 3/5/10/20,
read lengths 2/5/10 kb, three seeds, forward and inverted 1.2-kb shared blocks.
All cases contain an engineered insert; this matrix does not estimate specificity.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import itertools
import json
import math
from pathlib import Path
import random
import subprocess
import sys
import time

from hgtlong.sequence_qc import reverse_complement


def mutate(sequence,error,rng):
    out=[]
    for base in sequence:
        p=rng.random()
        if p<error*.2:continue
        if p<error*.4:out.append(rng.choice('ACGT'))
        out.append(rng.choice([b for b in 'ACGT' if b!=base]) if error*.4<=p<error else base)
    return ''.join(out)


def run_case(case,root,minimap2):
    error,depth,length,seed,inverted=case
    name=f'e{error}_d{depth}_l{length}_s{seed}_i{int(inverted)}'
    folder=root/name;folder.mkdir()
    rng=random.Random(3917+seed)
    genome_length=100000
    a=''.join(rng.choices('ACGT',k=genome_length))
    b=''.join(rng.choices('ACGT',k=genome_length))
    block=''.join(rng.choices('ACGT',k=1200))
    a=a[:40000]+block+a[41200:]
    b=b[:60000]+(reverse_complement(block) if inverted else block)+b[61200:]
    (folder/'panel.fa').write_text('>a\n'+a+'\n>b\n'+b+'\n')
    (folder/'panel.tsv').write_text('contig\tmag\na\tA\nb\tB\n')
    n=int(depth*genome_length/length)
    q=40 if error==0 else max(2,min(40,round(-10*math.log10(error))))
    with open(folder/'reads.fq','w') as f:
        for side,sequence in [('a',a),('b',b)]:
            starts=random.Random(8191+seed+(0 if side=='a' else 100))
            for i in range(n):
                start=starts.randrange(genome_length-length+1)
                read=mutate(sequence[start:start+length],error,random.Random(109+seed*100000+i+(0 if side=='a' else 50000)))
                if i%2:read=reverse_complement(read)
                f.write(f'@{side}{i}\n{read}\n+\n'+chr(33+q)*len(read)+'\n')
    start=time.monotonic()
    with open(folder/'reads.sam','w') as out,open(folder/'mapping.log','w') as err:
        subprocess.run([minimap2,'-t','1','-ax','map-ont','-p','0','-N','100','-Y',
                        str(folder/'panel.fa'),str(folder/'reads.fq')],stdout=out,stderr=err,check=True)
    with open(folder/'caller.log','w') as log:
        subprocess.run([sys.executable,'-m','hgtlong','extract','--sam',str(folder/'reads.sam'),
                        '--panel',str(folder/'panel.tsv'),'--out',str(folder/'evidence.tsv.gz'),
                        '--workers','1'],stdout=log,stderr=log,check=True)
        subprocess.run([sys.executable,'-m','hgtlong','call','--evidence',str(folder/'evidence.tsv.gz'),
                        '--outdir',str(folder/'calls')],stdout=log,stderr=log,check=True)
    with open(folder/'calls/events.tsv') as f:events=list(csv.DictReader(f,delimiter='\t'))
    supported=[e for e in events if e['status']=='read_supported']
    matched=[e for e in supported if e['orientation']==('-' if inverted else '+') and
             max(abs(int(e[k])-v) for k,v in zip(['a0','a1','b0','b1'],[40000,41200,60000,61200]))<=50]
    return dict(case=name,error_rate=error,nominal_depth=depth,read_length=length,seed=seed,
                inverted=int(inverted),reads=2*n,detected=int(bool(matched)),
                supported_candidates=len(supported),off_target_supported=len(supported)-len(matched),
                partial_candidates=sum(e['status']=='partial' for e in events),
                wall_seconds=round(time.monotonic()-start,3))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--outdir',required=True)
    parser.add_argument('--minimap2',required=True)
    args=parser.parse_args()
    root=Path(args.outdir).resolve();root.mkdir(parents=True)
    cases=list(itertools.product([0,.005,.01,.02],[3,5,10,20],[2000,5000,10000],range(3),[False,True]))
    (root/'design.json').write_text(json.dumps(dict(cases=len(cases),genome_length=100000,
        insert_length=1200,model='independent substitutions/insertions/deletions; not empirically calibrated',
        seed_design='genomes and read starts paired across conditions',
        interpretation='engineered-coordinate detection; not biological HGT sensitivity or specificity'),indent=2)+'\n')
    results=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for result in pool.map(lambda case:run_case(case,root,args.minimap2),cases):
            results.append(result)
            if len(results)%24==0:print(f'{len(results)}/{len(cases)} cases completed',flush=True)
    with open(root/'matrix.tsv','w') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(results)
    print('Complete:',len(results),'cases; detections',sum(r['detected'] for r in results),flush=True)


if __name__=='__main__':main()
