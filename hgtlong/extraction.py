"""Bounded, deterministic parallel extraction from query-grouped SAM."""
from collections import Counter,deque
from concurrent.futures import ProcessPoolExecutor
import itertools
import multiprocessing
import sys
from types import SimpleNamespace
from .core import Alignment,classify

_PANEL=None
_CONFIG=None


def initialize(panel,config):
    global _PANEL,_CONFIG
    _PANEL=panel;_CONFIG=config


def process_batch(batch):
    stats=Counter();evidence=[]
    for name,lines in batch:
        records=[]
        for line in lines:
            record=Alignment.from_sam(line)
            if record is None:stats['unmapped_records']+=1
            else:records.append(record)
        if records:
            stats['reads']+=1
            evidence.extend(classify(records,_PANEL,_CONFIG,stats))
    return evidence,stats


def stream_evidence(source,panel,args,stats,headers):
    config=SimpleNamespace(**{k:getattr(args,k) for k in (
        'full_cov','full_identity','frag_cov','frag_identity','min_mapq','min_block','max_indel','flank')})
    progress=0
    def absorb(counts):
        nonlocal progress
        stats.update(counts)
        bucket=stats.get('reads',0)//100000
        if bucket>progress:
            progress=bucket
            print('HGTlong: %s mapped read groups; %s evidence rows' %
                  (format(stats['reads'],','),format(stats.get('evidence_rows',0),',')),
                  file=sys.stderr,flush=True)
    def lines():
        for line in source:
            if line.startswith('@'):headers.append(line.rstrip())
            else:yield line
    def batches():
        seen=set();batch=[];size=0
        for name,group in itertools.groupby(lines(),lambda line:line.split('\t',1)[0]):
            if name in seen:
                raise ValueError('SAM is not query-name collated: repeated read group '+name)
            seen.add(name)
            records=list(group);batch.append((name,records));size+=sum(map(len,records))
            if len(batch)>=128 or size>=8*1024*1024:
                yield batch;batch=[];size=0
        if batch:yield batch
    if args.workers==1:
        initialize(panel,config)
        for batch in batches():
            rows,counts=process_batch(batch);absorb(counts)
            yield from rows
    else:
        with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,
                                 initargs=(panel,config),mp_context=multiprocessing.get_context('spawn')) as pool:
            pending=deque()
            for batch in batches():
                pending.append(pool.submit(process_batch,batch))
                if len(pending)>=2*args.workers:
                    rows,counts=pending.popleft().result();absorb(counts)
                    yield from rows
            while pending:
                rows,counts=pending.popleft().result();absorb(counts)
                yield from rows
