"""Experimental fixed BLAST-block validation using unclustered ONT evidence.

BLAST discovery is an explicit input, not an evolutionary truth label. Read
assignment and boundary requirements match the read-first caller, but no
fully bracketed read seed is required to propose an event.
"""
from collections import defaultdict
from .core import union_length
from .events import compatible, COORDS


def key(row):
    return tuple(row[k] for k in ('magA', 'ctgA', 'magB', 'ctgB', 'orientation'))


def boundary_supported(row, candidate, coordinate, flank=100, tolerance=50):
    position = candidate[coordinate]
    if not row['anchor_' + coordinate] or abs(row[coordinate] - position) > tolerance:
        return False
    if row['resident_start'] > position-flank or row['resident_end'] < position+flank:
        return False
    blocks = [tuple(map(int, b.split(':'))) for b in row['resident_blocks'].split(';') if b]
    return all(union_length([(max(s, lo), min(e, hi)) for s, e in blocks
                             if min(e, hi) > max(s, lo)]) >= .9*flank
               for lo, hi in [(position-flank, position), (position, position+flank)])


def validate(candidates, evidence, min_support=3, flank=100, tolerance=50):
    by_pair = defaultdict(list)
    states = {}
    for c in candidates:
        by_pair[key(c)].append(c)
        states[c['event_id']] = dict(reads=set(), ambiguous=set(), multi=set(),
                                     usable=set(), support={k: set() for k in COORDS})
    memberships = []
    for row in evidence:
        matches = [c for c in by_pair[key(row)] if compatible(row, c, tolerance)]
        for c in matches:
            state = states[c['event_id']]
            state['reads'].add(row['read'])
            if row['ambiguous_fragment']:
                state['ambiguous'].add(row['read'])
            if len(matches) != 1:
                state['multi'].add(row['read'])
            supported = []
            if not row['ambiguous_fragment'] and len(matches) == 1:
                state['usable'].add(row['read'])
                prefix = 'a' if row['resident'] == 'A' else 'b'
                for suffix in ('0', '1'):
                    coordinate = prefix+suffix
                    if boundary_supported(row, c, coordinate, flank, tolerance):
                        state['support'][coordinate].add(row['read'])
                        supported.append(coordinate)
            memberships.append(dict(event_id=c['event_id'], **row,
                                    assignment_ambiguous=int(len(matches) != 1),
                                    supported_boundaries=';'.join(supported)))
    events = []
    for c in candidates:
        state = states[c['event_id']]
        counts = {k: len(v) for k, v in state['support'].items()}
        ratio = min(c['a1']-c['a0'], c['b1']-c['b0'])/max(c['a1']-c['a0'], c['b1']-c['b0'])
        status = ('read_supported' if min(counts.values()) >= min_support and ratio >= .9
                  else 'no_qualifying_evidence' if not state['reads']
                  else 'ambiguous' if not state['usable'] else 'partial')
        events.append(dict(c, status=status, molecules=len(state['reads']),
                           ambiguous_molecules=len(state['ambiguous']),
                           unassigned_molecules=len(state['multi']), length_ratio=ratio,
                           transfer_direction='unresolved',
                           **{'n_'+k: v for k, v in counts.items()}))
    return events, memberships
