import argparse
from collections import Counter
import csv
import gzip
import hashlib
import html
import itertools
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from datetime import datetime, timezone

from . import __version__
from .core import Alignment, classify
from .events import call_events


def open_text(path, mode='rt'):
    return gzip.open(path,mode) if str(path).endswith('.gz') else open(path,mode)


def write_tsv(path, rows, fields=None):
    rows=iter(rows)
    first=next(rows,None)
    path=Path(path)
    suffix='.gz' if str(path).endswith('.gz') else '.tmp'
    descriptor,temporary=tempfile.mkstemp(prefix='.'+path.name+'.',suffix=suffix,dir=path.parent)
    os.close(descriptor)
    try:
        with open_text(temporary,'wt') as f:
            writer=csv.DictWriter(f,fieldnames=fields or (list(first) if first else ['event_id']),delimiter='\t',lineterminator='\n')
            writer.writeheader()
            if first is not None:
                writer.writerow(first)
            writer.writerows(rows)
        os.replace(temporary,path)
    finally:
        if os.path.exists(temporary):os.unlink(temporary)


def load_rows(path):
    with open_text(path) as f:
        for row in csv.DictReader(f,delimiter='\t'):
            for k in ('a0','a1','b0','b1','anchor_a0','anchor_a1','anchor_b0','anchor_b1',
                      'resident_start','resident_end','ambiguous_fragment'):
                if k in row: row[k]=int(row[k])
            yield row


def panel_map(path):
    with open(path) as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows or not {'contig','mag'} <= rows[0].keys():
        raise ValueError('Panel must have contig and mag columns')
    result={}
    for r in rows:
        if r['contig'] in result:
            raise ValueError('Duplicate panel contig: '+r['contig'])
        result[r['contig']]=r['mag']
    return result


def extract(args):
    panel=panel_map(args.panel)
    stats=Counter()
    headers=[]
    if args.sam == '-':
        source=sys.stdin
    else:
        source=open(args.sam)
    from .extraction import stream_evidence
    write_tsv(args.out,stream_evidence(source,panel,args,stats,headers),fields=[
        'read','ctgA','magA','a0','a1','ctgB','magB','b0','b1','orientation','resident',
        'anchor_a0','anchor_a1','anchor_b0','anchor_b1','resident_start','resident_end',
        'resident_identity','fragment_identity','ambiguous_fragment','resident_blocks'])
    if source is not sys.stdin: source.close()
    Path(args.out+'.stats.json').write_text(json.dumps(dict(stats),indent=2)+'\n')
    manifest=dict(version=__version__,evidence_schema='1.1',parameters={k:v for k,v in vars(args).items() if k!='func'},
                  sam_header=headers,panel_sha256=hashlib.sha256(Path(args.panel).read_bytes()).hexdigest(),
                  source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted(Path(__file__).parent.glob('*.py'))},
                  created_utc=datetime.now(timezone.utc).isoformat())
    Path(args.out+'.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


def report(events, path):
    counts=Counter(r['status'] for r in events)
    columns=['event_id','status','ctgA','a0','a1','ctgB','b0','b1','orientation',
             'molecules','n_a0','n_a1','n_b0','n_b1','length_ratio']
    table='<tr>'+''.join('<th>'+c+'</th>' for c in columns)+'</tr>'
    for r in events:
        table+='<tr>'+''.join('<td>'+html.escape(str(r[c]))+'</td>' for c in columns)+'</tr>'
    Path(path).write_text('''<!doctype html><meta charset="utf-8"><title>HGTlong evidence</title>
<style>body{font:15px system-ui;margin:32px;color:#182433}table{border-collapse:collapse;font-size:12px}td,th{padding:7px;border:1px solid #ddd}th{position:sticky;top:0;background:#edf3f8}input{padding:9px;width:400px}</style>
<h1>HGTlong</h1><p>Read-resolved candidate loci. Coordinates: 0-based, half-open.</p>
<p>Read support is not proof of horizontal transfer. Transfer direction, repeats, MAG contamination,
and evolutionary alternatives require further assessment. HDMI is an external comparator only.</p>
<p>'''+html.escape(json.dumps(counts))+'''</p><input id="filter" placeholder="Filter by contig, event, or status">
<table id="events">'''+table+'''</table><script>document.getElementById('filter').oninput=function(){for(const r of document.querySelectorAll('#events tr')){if(r.querySelector('td'))r.hidden=!r.textContent.toLowerCase().includes(this.value.toLowerCase());}};</script>''')


def filter_lengths(events, members, out, minimum):
    excluded=[e for e in events if minimum is not None and min(e['a1']-e['a0'],e['b1']-e['b0'])<minimum]
    rejected={e['event_id'] for e in excluded}
    write_tsv(out/'length_excluded.tsv',excluded,fields=list(events[0]) if events else ['event_id'])
    (out/'length_filter.json').write_text(json.dumps(dict(min_hgt_length=minimum,rule='both_blocks_at_least_cutoff',excluded=len(excluded),retained=len(events)-len(excluded)),indent=2)+'\n')
    return [e for e in events if e['event_id'] not in rejected],[m for m in members if m['event_id'] not in rejected]


def call(args):
    events,members=call_events(load_rows(args.evidence),args)
    out=Path(args.outdir)
    out.mkdir(parents=True,exist_ok=True)
    events,members=filter_lengths(events,members,out,args.min_hgt_length)
    event_fields=['event_id','status','magA','ctgA','a0','a1','magB','ctgB','b0','b1',
                  'orientation','molecules','resident_A','resident_B','n_a0','n_a1','n_b0','n_b1',
                  'length_ratio','unresolved','ambiguous_molecules','unassigned_molecules','repeat_screen',
                  'evolutionary_status','transfer_direction']
    write_tsv(out/'events.tsv',events,fields=event_fields)
    write_tsv(out/'event_evidence.tsv.gz',members)
    from .reporting import write_report
    write_report(events,members,out/'report.html')
    params={k:v for k,v in vars(args).items() if k!='func'}
    manifest=dict(version=__version__,created_utc=datetime.now(timezone.utc).isoformat(),
                  parameters=params,status_counts=dict(Counter(e['status'] for e in events)),
                  coordinate_system='0-based-half-open',comparator_used_for_calling=False,
                  evidence_sha256=hashlib.sha256(Path(args.evidence).read_bytes()).hexdigest())
    manifest['source_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(Path(__file__).parent.glob('*.py'))}
    (out/'run.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest['status_counts']))


def compare(args):
    calls=list(load_rows(args.events))
    with open(args.hdmi) as f:
        refs=list(csv.DictReader(f,delimiter='\t'))
    matches=[]
    for ref in refs:
        a,b=ref['ONT_query_Reference_ID'],ref['ONT_subject_Reference_ID']
        x=sorted([int(ref['ONT_query_Start_1based']),int(ref['ONT_query_End_1based'])])
        y=sorted([int(ref['ONT_subject_Start_1based']),int(ref['ONT_subject_End_1based'])])
        x[0]-=1; y[0]-=1
        expected_orientation=None
        if all(ref.get(k) for k in ('q_start','q_end','s_start','s_end')):
            expected_orientation='-' if ((int(ref['q_end'])-int(ref['q_start']))*
                                         (int(ref['s_end'])-int(ref['s_start'])))<0 else '+'
        for call in calls:
            if (call['ctgA'],call['ctgB'])==(a,b): xx,yy=x,y
            elif (call['ctgA'],call['ctgB'])==(b,a): xx,yy=y,x
            else: continue
            da=max(abs(call['a0']-xx[0]),abs(call['a1']-xx[1]))
            db=max(abs(call['b0']-yy[0]),abs(call['b1']-yy[1]))
            overlaps=min(call['a1'],xx[1])>max(call['a0'],xx[0]) and min(call['b1'],yy[1])>max(call['b0'],yy[0])
            tier='exact' if max(da,db)==0 else 'within100' if max(da,db)<=100 else 'within500' if max(da,db)<=500 else 'bilateral_overlap' if overlaps else 'same_contig_pair'
            geometry_tier=tier
            orientation_match='unknown'
            if expected_orientation and call.get('orientation') in ('+','-'):
                orientation_match=str(call['orientation']==expected_orientation).lower()
                if orientation_match=='false':tier='orientation_conflict'
            matches.append(dict(hdmi_id=ref['HGT_ID'],event_id=call['event_id'],status=call['status'],tier=tier,
                                geometry_tier=geometry_tier,orientation_match=orientation_match,max_boundary_distance=max(da,db)))
    write_tsv(args.out,matches,fields=['hdmi_id','event_id','status','tier','geometry_tier','orientation_match','max_boundary_distance'])
    summary={'hdmi_events':len(refs),'hgtlong_candidates':len(calls),'interpretation':'concordance, not sensitivity or precision','tiers':{}}
    for tier in ['exact','within100','within500','bilateral_overlap','same_contig_pair','orientation_conflict']:
        m=[r for r in matches if r['tier']==tier]
        summary['tiers'][tier]={'hdmi_events':len({r['hdmi_id'] for r in m}),'hgtlong_candidates':len({r['event_id'] for r in m})}
    summary['cumulative_boundary_tolerance']={}
    for tolerance in (0,100,500):
        m=[r for r in matches if r['max_boundary_distance']<=tolerance and r['orientation_match']!='false']
        summary['cumulative_boundary_tolerance'][str(tolerance)]={
            'hdmi_events':len({r['hdmi_id'] for r in m}),
            'hgtlong_candidates':len({r['event_id'] for r in m}),
            'read_supported_candidates':len({r['event_id'] for r in m if r['status']=='read_supported'})}
    overlap=[r for r in matches if r['tier'] not in ('same_contig_pair','orientation_conflict')]
    summary['bilateral_overlap_including_boundary_matches']={
        'hdmi_events':len({r['hdmi_id'] for r in overlap}),
        'hgtlong_candidates':len({r['event_id'] for r in overlap})}
    summary['orientation_policy']='Require agreement when both orientations are available; unknown is reported explicitly'
    Path(args.out+'.summary.json').write_text(json.dumps(summary,indent=2)+'\n')


def qc(args):
    from .qc import annotate_masks
    with open(args.screened_contigs) as f:
        screened={line.strip() for line in f if line.strip()}
    rows=list(load_rows(args.events))
    if not rows:
        with open(args.events) as f:
            fields=next(csv.reader(f,delimiter='\t'))
        write_tsv(args.out,[],fields=list(dict.fromkeys(fields+['mask_bp_A','mask_bp_B','qc_status','mask_classes'])))
    else:
        write_tsv(args.out,annotate_masks(rows,args.mask_bed,screened))


def sequence_qc(args):
    from .sequence_qc import annotate_sequences
    rows=list(load_rows(args.events))
    if rows:
        write_tsv(args.out,annotate_sequences(rows,args.reference,args.minimap2))
    else:
        with open(args.events) as f:fields=next(csv.reader(f,delimiter='\t'))
        write_tsv(args.out,[],fields=list(dict.fromkeys(fields+['block_identity','aligned_fraction_A','aligned_fraction_B','sequence_qc'])))


def build_panel(args):
    from .panel import build_panel as build
    print(json.dumps(build(args.genomes,args.outdir)))


def screen(args):
    from .screening import screen_event
    rows=[screen_event(r,args.min_block_identity,args.min_block_coverage,args.min_boundary_reads)
          for r in load_rows(args.events)]
    for row in rows:
        if args.min_hgt_length is not None and min(row['a1']-row['a0'],row['b1']-row['b0'])<args.min_hgt_length:
            row['screened_candidate']=0
            row['screen_fail_reasons']=';'.join(filter(None,[row['screen_fail_reasons'],'short_hgt_block']))
    with open(args.events) as f:fields=next(csv.reader(f,delimiter='\t'))
    fields=list(dict.fromkeys(fields+['screened_candidate','screen_fail_reasons']))
    write_tsv(args.out,rows,fields=fields)
    candidates=[r for r in rows if r['screened_candidate']]
    write_tsv(args.candidates,candidates,fields=fields)
    Path(args.out+'.criteria.json').write_text(json.dumps(dict(
        min_hgt_length=args.min_hgt_length,min_block_identity=args.min_block_identity,min_block_coverage=args.min_block_coverage,
        min_boundary_reads=args.min_boundary_reads,requires='read_supported + annotation_screen_pass + both_blocks_covered',
        comparator_used=False,interpretation='screened candidates, not confirmed evolutionary transfers',
        input_candidates=len(rows),screened_candidates=len(candidates)),indent=2)+'\n')


def call_blocks(args):
    from .blast_seeded import validate
    from .reporting import write_report
    candidates=list(load_rows(args.blocks))
    ids=[c['event_id'] for c in candidates]
    if len(ids)!=len(set(ids)):
        raise ValueError('Candidate event_id values must be unique')
    for c in candidates:
        if c['a0']<0 or c['b0']<0 or c['a1']<=c['a0'] or c['b1']<=c['b0'] or c['orientation'] not in ('+','-'):
            raise ValueError('Invalid block coordinates or orientation')
    events,members=validate(candidates,load_rows(args.evidence),args.min_support,args.flank,args.junction_tolerance)
    out=Path(args.outdir);out.mkdir(parents=True,exist_ok=True)
    events,members=filter_lengths(events,members,out,args.min_hgt_length)
    write_tsv(out/'events.tsv',events,fields=list(candidates[0])+['status','molecules','ambiguous_molecules','unassigned_molecules','length_ratio','transfer_direction','n_a0','n_a1','n_b0','n_b1'] if candidates else ['event_id','status'])
    write_tsv(out/'event_evidence.tsv.gz',members)
    write_report(events,members,out/'report.html')
    (out/'run.json').write_text(json.dumps(dict(version=__version__,mode='fixed-block-validation',
        parameters={k:v for k,v in vars(args).items() if k!='func'},
        blocks_sha256=hashlib.sha256(Path(args.blocks).read_bytes()).hexdigest(),
        evidence_sha256=hashlib.sha256(Path(args.evidence).read_bytes()).hexdigest(),
        status_counts=dict(Counter(e['status'] for e in events))),indent=2)+'\n')


def main():
    parser=argparse.ArgumentParser(prog='hgtlong')
    parser.add_argument('--version',action='version',version=__version__)
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('build-panel',help='Namespace and checksum one FASTA per MAG')
    p.add_argument('--genomes',required=True)
    p.add_argument('--outdir',required=True)
    p.set_defaults(func=build_panel)
    p=sub.add_parser('extract',help='Read evidence from query-name collated SAM')
    p.add_argument('--sam',required=True)
    p.add_argument('--panel',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--workers',type=int,default=max(1,min(8,int(os.environ.get('SLURM_CPUS_PER_TASK','3'))-2)),
                   help='Extraction workers; defaults to min(8, Slurm CPUs minus 2), otherwise 1')
    for name,default in [('full-cov',.90),('full-identity',.99),('frag-cov',.75),('frag-identity',.95)]:
        p.add_argument('--'+name,type=float,default=default)
    for name,default in [('min-mapq',20),('min-block',200),('max-indel',500),('flank',100)]:
        p.add_argument('--'+name,type=int,default=default)
    p.set_defaults(func=extract)
    p=sub.add_parser('call')
    p.add_argument('--evidence',required=True)
    p.add_argument('--outdir',required=True)
    p.add_argument('--junction-tolerance',type=int,default=50)
    p.add_argument('--min-support',type=int,default=3)
    p.add_argument('--flank',type=int,default=100)
    p.add_argument('--min-length-ratio',type=float,default=.9)
    p.set_defaults(func=call)
    p=sub.add_parser('call-blocks',help='Validate supplied fixed blocks with ONT evidence')
    p.add_argument('--blocks',required=True)
    p.add_argument('--evidence',required=True)
    p.add_argument('--outdir',required=True)
    p.add_argument('--min-support',type=int,default=3)
    p.add_argument('--flank',type=int,default=100)
    p.add_argument('--junction-tolerance',type=int,default=50)
    p.set_defaults(func=call_blocks)
    p=sub.add_parser('compare',help='HDMI concordance; never a ground-truth benchmark')
    p.add_argument('--events',required=True)
    p.add_argument('--hdmi',required=True)
    p.add_argument('--out',required=True)
    p.set_defaults(func=compare)
    p=sub.add_parser('qc',help='Annotate supplied repeat masks after discovery')
    p.add_argument('--events',required=True)
    p.add_argument('--mask-bed',required=True)
    p.add_argument('--screened-contigs',required=True)
    p.add_argument('--out',required=True)
    p.set_defaults(func=qc)
    p=sub.add_parser('sequence-qc',help='Compare already called block sequences; never changes partners or calls')
    p.add_argument('--events',required=True)
    p.add_argument('--reference',required=True)
    p.add_argument('--minimap2',default='minimap2')
    p.add_argument('--out',required=True)
    p.set_defaults(func=sequence_qc)
    p=sub.add_parser('screen',help='Separate explicitly screened candidates from all read-derived loci')
    p.add_argument('--events',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--candidates',required=True)
    p.add_argument('--min-block-identity',type=float,default=.99)
    p.add_argument('--min-block-coverage',type=float,default=.95)
    p.add_argument('--min-boundary-reads',type=int,default=3)
    p.set_defaults(func=screen)
    for command in ('call','call-blocks','screen'):
        sub.choices[command].add_argument('--min-hgt-length',type=int,default=None,help='Minimum length in bp of EACH HGT block; default: no cutoff')
    args=parser.parse_args()
    if getattr(args,'min_hgt_length',None) is not None and args.min_hgt_length<1:
        parser.error('min_hgt_length must be positive')
    for key in ('full_cov','full_identity','frag_cov','frag_identity','min_length_ratio','min_block_identity','min_block_coverage'):
        if hasattr(args,key) and not 0<getattr(args,key)<=1:
            parser.error(key+' must be in (0,1]')
    for key in ('min_block','flank','min_support','workers','min_boundary_reads'):
        if hasattr(args,key) and getattr(args,key)<1:
            parser.error(key+' must be positive')
    for key in ('max_indel','junction_tolerance','min_mapq'):
        if hasattr(args,key) and getattr(args,key)<0:
            parser.error(key+' must be nonnegative')
    if hasattr(args,'min_mapq') and args.min_mapq>=255:
        parser.error('min_mapq must be below 255 (unknown MAPQ)')
    args.func(args)


if __name__=='__main__': main()
