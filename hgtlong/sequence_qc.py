"""Align ONLY the two already called blocks; never search for or replace partners."""
from pathlib import Path
import subprocess
import tempfile


def read_fasta(path):
    sequences={};name=None;parts=[]
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:sequences[name]=''.join(parts).upper()
                name=line[1:].split()[0];parts=[]
                if name in sequences:raise ValueError('Duplicate FASTA contig: '+name)
            elif name is not None:parts.append(line.strip())
        if name is not None:sequences[name]=''.join(parts).upper()
    return sequences


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN'))[::-1]


def metrics(paf):
    hits=[]
    for line in paf.splitlines():
        f=line.split('\t')
        if len(f)<12:raise ValueError('Malformed minimap2 PAF')
        hits.append((int(f[9]),int(f[10]),int(f[3])-int(f[2]),int(f[1]),
                     int(f[8])-int(f[7]),int(f[6]),f[4]))
    if not hits:return dict(block_identity='',aligned_fraction_A='',aligned_fraction_B='',sequence_qc='no_local_alignment')
    matches,length,qa,ql,ta,tl,strand=max(hits)
    return dict(block_identity=round(matches/length,6),aligned_fraction_A=round(qa/ql,6),
                aligned_fraction_B=round(ta/tl,6),
                sequence_qc='orientation_conflict' if strand!='+' else
                            'both_blocks_covered' if min(qa/ql,ta/tl)>=.95 else 'partial_local_alignment')


def annotate_sequences(events,reference,minimap2):
    sequences=read_fasta(reference)
    with tempfile.TemporaryDirectory(prefix='hgtlong-sequence-qc-') as d:
        query,target=Path(d)/'A.fa',Path(d)/'B.fa'
        for event in events:
            result=dict(block_identity='',aligned_fraction_A='',aligned_fraction_B='',sequence_qc='not_assessed')
            if event['status']=='read_supported':
                blocks=[]
                for side in ['A','B']:
                    seq=sequences[event['ctg'+side]]
                    start,end=int(event[side.lower()+'0']),int(event[side.lower()+'1'])
                    if not 0<=start<end<=len(seq):raise ValueError('Called coordinates outside reference')
                    blocks.append(seq[start:end])
                if event['orientation']=='-':blocks[1]=reverse_complement(blocks[1])
                query.write_text('>A\n'+blocks[0]+'\n');target.write_text('>B\n'+blocks[1]+'\n')
                command=[minimap2,'-c','-x','asm5','-k','15','-w','5','-m','40','-s','80',
                         '--secondary=no','-t','1',str(target),str(query)]
                p=subprocess.run(command,capture_output=True,text=True,check=True,timeout=120)
                result=metrics(p.stdout)
            yield dict(event,**result)
