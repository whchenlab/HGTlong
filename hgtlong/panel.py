"""Build an explicitly namespaced, checksummed MAG reference panel."""
import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import os


def build_panel(genomes,outdir):
    genomes=Path(genomes).resolve();outdir=Path(outdir)
    files=sorted(p for p in genomes.iterdir() if p.is_file() and
                 any(p.name.endswith(ext) for ext in ('.fa','.fna','.fasta','.fa.gz','.fna.gz','.fasta.gz')))
    if not files:raise ValueError('No MAG FASTA files found')
    outdir.mkdir(parents=True,exist_ok=True)
    names=['panel.fa','panel.tsv','panel_manifest.json']
    if any((outdir/name).exists() for name in names):raise ValueError('Panel outputs already exist; use a new output directory')
    seen=set();mags=set();metadata=[];count=0;bases=0
    with tempfile.TemporaryDirectory(prefix='.panel-',dir=outdir) as temp:
        temp=Path(temp)
        with open(temp/'panel.fa','w') as fasta,open(temp/'panel.tsv','w') as table:
            table.write('contig\tmag\toriginal_contig\tlength\n')
            for path in files:
                name=path.name[:-3] if path.name.endswith('.gz') else path.name
                mag=Path(name).stem
                if mag in mags or any(c.isspace() for c in mag):raise ValueError('Invalid/duplicate MAG ID: '+mag)
                mags.add(mag);length=0;contig=None;original=None
                opener=gzip.open if path.name.endswith('.gz') else open
                with opener(path,'rt') as source:
                    for line in source:
                        if line.startswith('>'):
                            if contig is not None:
                                if length==0:raise ValueError('Empty contig: '+contig)
                                table.write(f'{contig}\t{mag}\t{original}\t{length}\n')
                            tokens=line[1:].split()
                            if not tokens:raise ValueError('Empty FASTA header')
                            original=tokens[0];contig=mag+'__'+original;length=0
                            if contig in seen:raise ValueError('Duplicate contig: '+contig)
                            seen.add(contig);count+=1;fasta.write('>'+contig+'\n')
                        elif line.strip():
                            if contig is None:raise ValueError('Sequence before FASTA header')
                            sequence=line.strip().upper()
                            if set(sequence)-set('ACGTRYSWKMBDHVN'):raise ValueError('Non-DNA sequence in '+str(path))
                            length+=len(sequence);bases+=len(sequence);fasta.write(sequence+'\n')
                    if contig is None or length==0:raise ValueError('Empty MAG/contig: '+str(path))
                    table.write(f'{contig}\t{mag}\t{original}\t{length}\n')
                digest=hashlib.sha256()
                with open(path,'rb') as source:
                    for block in iter(lambda:source.read(1024*1024),b''):digest.update(block)
                metadata.append(dict(mag=mag,path=str(path),sha256=digest.hexdigest()))
        (temp/'panel_manifest.json').write_text(json.dumps(dict(mags=len(mags),contigs=count,
            bases=bases,inputs=metadata),indent=2)+'\n')
        for name in names:os.replace(temp/name,outdir/name)
    return dict(mags=len(mags),contigs=count,bases=bases)
