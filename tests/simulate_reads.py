"""Deterministic sequence-level smoke tests, not an ONT performance benchmark.

Generate engineered shared inserts, independent random genomes, and ligation chimeras.
Mapping is performed independently by minimap2 in the queued workflow.
"""
import json
from pathlib import Path
import random
import sys


def rc(s):
    return s.translate(str.maketrans('ACGT','TGCA'))[::-1]


def mutate(s,rng):
    out=[]
    for base in s:
        p=rng.random()
        if p<.001: continue
        if p<.002: out.append(rng.choice('ACGT'))
        out.append(rng.choice([b for b in 'ACGT' if b!=base]) if .002<=p<.005 else base)
    return ''.join(out)


def generate(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    cases=[]
    for kind in ['insert','inverted_insert','null','chimera','duplicate_molecule','long_insert']:
        for rep in range(3):
            rng=random.Random(9235+rep)
            dna=lambda n: ''.join(rng.choices('ACGT',k=n))
            a,b=dna(30000),dna(30000)
            size=8000 if kind=='long_insert' else 1200
            shared=dna(size)
            if kind in ['insert','inverted_insert','duplicate_molecule','long_insert']:
                a=a[:8000]+shared+a[8000+size:]
                b=b[:14000]+(rc(shared) if kind=='inverted_insert' else shared)+b[14000+size:]
            name=f'{kind}_{rep}'
            folder=out/name;folder.mkdir(exist_ok=True)
            (folder/'panel.fa').write_text('>a\n'+a+'\n>b\n'+b+'\n')
            (folder/'panel.tsv').write_text('contig\tmag\na\tA\nb\tB\n')
            reads=[]
            for ctg,seq,position in [('a',a,8000),('b',b,14000)]:
                for j in range(16):
                    start=position-rng.randint(400,1600)
                    if kind=='long_insert' and j%2:
                        start=position+size-rng.randint(400,1600)
                    read=seq[start:start+4000]
                    if kind=='chimera':
                        read=a[start:start+2000]+b[14000:16000]
                    read=mutate(read,rng)
                    if j%2: read=rc(read)
                    reads.append((ctg+str(j),read))
            if kind=='duplicate_molecule':
                reads=[('a0',reads[0][1])]*16+[('b0',reads[16][1])]*16
            with open(folder/'reads.fq','w') as f:
                for read,seq in reads: f.write('@'+read+'\n'+seq+'\n+\n'+'I'*len(seq)+'\n')
            cases.append(dict(name=name,kind=kind,replicate=rep,
                              expect_supported=kind in ['insert','inverted_insert'],
                              coordinates=[8000,8000+size,14000,14000+size]))
    (out/'truth.json').write_text(json.dumps(cases,indent=2)+'\n')


if __name__=='__main__': generate(sys.argv[1])
