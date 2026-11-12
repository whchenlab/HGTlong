"""Conservative junction clustering with explicit censored boundaries."""
from collections import defaultdict
import hashlib
import statistics
from .core import union_length

COORDS = ('a0','a1','b0','b1')


def compatible(row, bounds, tolerance):
    for k in COORDS:
        if row['anchor_'+k] and abs(row[k]-bounds[k]) > tolerance:
            return False
    # Censored alignments must still overlap both inferred loci.
    return all(min(row[end],bounds[end])-max(row[start],bounds[start]) > 0
               for start,end in [('a0','a1'),('b0','b1')])


def call_events(rows, config):
    grouped=defaultdict(list)
    for row in rows:
        if 'resident_blocks' not in row:
            raise ValueError('Evidence lacks resident_blocks; re-extract with the current HGTlong version')
        grouped[(row['magA'],row['ctgA'],row['magB'],row['ctgB'],row['orientation'])].append(row)
    events, memberships = [], []
    for key, group in sorted(grouped.items()):
        # Fully bracketed fragments define seeds. Never merge by interval extent.
        seeds=[]
        complete = [r for r in group if all(r['anchor_'+k] for k in COORDS)]
        for row in sorted(complete,key=lambda r:tuple(r[k] for k in COORDS)+(r['read'],)):
            matches=[s for s in seeds if all(
                max(s['max'][k],row[k])-min(s['min'][k],row[k]) <= config.junction_tolerance
                for k in COORDS)]
            if len(matches)==1:
                seed=matches[0]
                seed['rows'].append(row)
                for k in COORDS:
                    seed['min'][k]=min(seed['min'][k],row[k])
                    seed['max'][k]=max(seed['max'][k],row[k])
            elif not matches:
                seeds.append(dict(rows=[row],min={k:row[k] for k in COORDS},
                                  max={k:row[k] for k in COORDS}))
        for seed in seeds:
            seed['bounds']={k:int(statistics.median(r[k] for r in seed['rows'])) for k in COORDS}
            seed['assigned']=[]
        unresolved=[]
        for row in group:
            matches=[s for s in seeds if compatible(row,s['bounds'],config.junction_tolerance)]
            if len(matches)==1:
                matches[0]['assigned'].append(row)
            else:
                unresolved.append(row)
        # Retain exact unresolved geometries as partial candidates; no fabricated endpoints.
        ungroup=defaultdict(list)
        for row in unresolved:
            ungroup[tuple(row[k] for k in COORDS)].append(row)
        candidates=[(s['bounds'],s['assigned'],False) for s in seeds]
        candidates += [(dict(zip(COORDS,k)),v,True) for k,v in sorted(ungroup.items())]
        # An unresolved geometry can coincide with a resolved seed while being
        # compatible with several seeds. Keep one canonical locus, but never
        # let those ambiguously assigned molecules contribute boundary support.
        canonical=defaultdict(list)
        for bounds, evidence, unresolved_flag in candidates:
            canonical[tuple(bounds[k] for k in COORDS)].append((evidence,unresolved_flag))
        candidates=[]
        for coordinates, parts in sorted(canonical.items()):
            resolved=any(not flag for _,flag in parts)
            evidence=[dict(row,assignment_ambiguous=int(flag and resolved))
                      for rows,flag in parts for row in rows]
            candidates.append((dict(zip(COORDS,coordinates)),evidence,not resolved))
        for bounds, evidence, unresolved_flag in candidates:
            unique={r['read'] for r in evidence}
            if len(unique)<config.min_support:
                continue
            support={k:set() for k in COORDS}
            usable=[r for r in evidence if not r['ambiguous_fragment'] and not r['assignment_ambiguous']]
            for row in usable:
                prefix='a' if row['resident']=='A' else 'b'
                for suffix in ['0','1']:
                    k=prefix+suffix
                    position=bounds[k]
                    blocks=[tuple(map(int,t.split(':'))) for t in row.get('resident_blocks','').split(';') if t]
                    aligned=(min(union_length([(max(s,lo),min(e,hi)) for s,e in blocks
                                                if min(e,hi)>max(s,lo)])/config.flank
                                 for lo,hi in [(position-config.flank,position),
                                               (position,position+config.flank)])>=.9) if blocks else False
                    if (row['anchor_'+k] and abs(row[k]-position)<=config.junction_tolerance
                        and aligned
                        and row['resident_start'] <= position-config.flank
                        and row['resident_end'] >= position+config.flank):
                        support[k].add(row['read'])
            ratio=min(bounds['a1']-bounds['a0'],bounds['b1']-bounds['b0'])/max(
                bounds['a1']-bounds['a0'],bounds['b1']-bounds['b0'])
            status='read_supported' if (not unresolved_flag and
                min(map(len,support.values()))>=config.min_support and ratio>=config.min_length_ratio) else 'partial'
            if not usable:
                status='ambiguous'
            identity='|'.join(map(str,key+tuple(bounds[k] for k in COORDS)))
            eid='HGL_'+hashlib.sha256(identity.encode()).hexdigest()[:16]
            event=dict(zip(('magA','ctgA','magB','ctgB','orientation'),key))
            event.update(bounds)
            event.update(event_id=eid,status=status,molecules=len(unique),
                         resident_A=len({r['read'] for r in usable if r['resident']=='A'}),
                         resident_B=len({r['read'] for r in usable if r['resident']=='B'}),
                         length_ratio=round(ratio,6), unresolved=int(unresolved_flag),
                         ambiguous_molecules=len({r['read'] for r in evidence if r['ambiguous_fragment']}),
                         unassigned_molecules=len({r['read'] for r in evidence if r['assignment_ambiguous']}),
                         repeat_screen='not_assessed',evolutionary_status='candidate_only',
                         transfer_direction='unresolved')
            event.update({'n_'+k:len(support[k]) for k in COORDS})
            events.append(event)
            for row in evidence:
                memberships.append(dict(event_id=eid,**row))
    return sorted(events,key=lambda e:e['event_id']),memberships
