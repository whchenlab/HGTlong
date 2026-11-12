"""Alignment geometry. All coordinates are zero-based, half-open.

An alignment's query coordinates are always in the original molecule orientation.
No assembly-to-assembly comparison or comparator table enters this module.
"""
import re
from dataclasses import dataclass
from functools import cached_property

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def union_length(intervals):
    end, total = -1, 0
    for s, e in sorted(intervals):
        total += max(0, e - max(s, end))
        end = max(end, e)
    return total


@dataclass(frozen=True)
class Alignment:
    read: str
    flag: int
    contig: str
    start: int
    mapq: int
    cigar: tuple
    nm: int

    @classmethod
    def from_sam(cls, line):
        f = line.rstrip().split("\t")
        if len(f) < 11:
            raise ValueError("Malformed SAM record")
        if int(f[1]) & 4:
            return None
        ops = tuple((int(n), op) for n, op in CIGAR.findall(f[5]))
        if not ops or ''.join(str(n) + op for n, op in ops) != f[5]:
            raise ValueError("Malformed CIGAR: " + f[5])
        nm = next((int(t[5:]) for t in f[11:] if t.startswith('NM:i:')), -1)
        return cls(f[0], int(f[1]), f[2], int(f[3])-1, int(f[4]), ops, nm)

    @cached_property
    def length(self):
        return sum(n for n, op in self.cigar if op in 'MISH=X')

    @cached_property
    def reverse(self):
        return bool(self.flag & 16)

    @cached_property
    def end(self):
        return self.start + sum(n for n, op in self.cigar if op in 'MDN=X')

    @cached_property
    def identity(self):
        denominator = sum(n for n, op in self.cigar if op in 'MID=X')
        return max(0., 1-self.nm/denominator) if self.nm >= 0 and denominator else 0.

    @cached_property
    def query_interval(self):
        left = 0
        for n, op in self.cigar:
            if op not in 'SH':
                break
            left += n
        size = sum(n for n, op in self.cigar if op in 'MI=X')
        right = left + size
        return (self.length-right, self.length-left) if self.reverse else (left, right)

    @cached_property
    def max_gap(self):
        return max((n for n, op in self.cigar if op in 'IDN'), default=0)

    def project(self, original_query_position):
        """Project a query boundary through CIGAR; insertion interiors are unknown.

        Exact deletion boundaries use the left reference coordinate. Such uncertainty
        is bounded by the configured maximum indel and reported as an approximation.
        """
        target = self.length-original_query_position if self.reverse else original_query_position
        q, r = 0, self.start
        for n, op in self.cigar:
            if op in 'M=X':
                if q <= target <= q+n:
                    return r+target-q
                q += n
                r += n
            elif op in 'ISH':
                if q < target < q+n:
                    return None
                q += n
            elif op in 'DN':
                r += n
        return r if target == q else None

    @cached_property
    def reference_blocks(self):
        blocks=[]
        r=self.start
        for n,op in self.cigar:
            if op in 'M=X':
                if blocks and blocks[-1][1]==r:
                    blocks[-1]=(blocks[-1][0],r+n)
                else:
                    blocks.append((r,r+n))
                r+=n
            elif op in 'DN':
                r+=n
        return blocks

    def anchor_coverage(self, position, flank):
        """Require actual aligned reference bases on both sides, not just a span."""
        return min(sum(max(0,min(e,hi)-max(s,lo)) for s,e in self.reference_blocks)/flank
                   for lo,hi in [(position-flank,position),(position,position+flank)])


def classify(records, panel, config, stats):
    """Conservative MVP: one intact resident alignment, no stitched full hits."""
    records = list(dict.fromkeys(r for r in records if not r.flag & (512 | 1024)))
    if not records:
        return []
    if len({r.length for r in records}) != 1:
        stats['inconsistent_read_lengths'] += 1
        return []
    length = records[0].length
    stats['missing_NM_records'] += sum(r.nm<0 for r in records)
    stats['reads_with_100_or_more_secondary_records'] += sum(bool(r.flag&256) for r in records)>=100
    primary=[r for r in records if not r.flag&(256|2048)
             and (r.query_interval[1]-r.query_interval[0])/length>=config.full_cov
             and config.min_mapq<=r.mapq<255]
    if primary:
        stats['primary_identity_percent_'+str(int(max(r.identity for r in primary)*100))]+=1
    full = [r for r in records
            if (r.query_interval[1]-r.query_interval[0])/length >= config.full_cov
            and r.identity >= config.full_identity and r.max_gap <= config.max_indel]
    # Count alternative intact loci even if their MAPQ is zero.
    loci = {(r.contig, r.start, r.end, r.reverse) for r in full}
    if len(loci) != 1:
        stats['ambiguous_resident' if loci else 'no_resident'] += 1
        return []
    resident = max(full, key=lambda r:(r.mapq,r.identity,-r.flag,r.cigar))
    if not config.min_mapq <= resident.mapq < 255:
        stats['resident_mapq_fail'] += 1
        return []
    by_mag = {}
    for r in records:
        if r.contig not in panel:
            raise ValueError('Contig absent from panel: '+r.contig)
        by_mag.setdefault(panel[r.contig], []).append(r)
    out = []
    for mag, hits in sorted(by_mag.items()):
        if mag == panel[resident.contig]:
            continue
        if union_length([r.query_interval for r in hits])/length > config.frag_cov:
            continue
        eligible = [r for r in hits if r.identity >= config.frag_identity
                    and r.query_interval[1]-r.query_interval[0] >= config.min_block
                    and r.max_gap <= config.max_indel]
        placements={}
        for r in eligible:
            key=(r.contig,r.start,r.end,r.query_interval,r.reverse)
            previous=placements.get(key)
            if previous is None or (r.identity,r.mapq,-r.flag)>(previous.identity,previous.mapq,-previous.flag):
                placements[key]=r
        eligible=[placements[k] for k in sorted(placements)]
        # Keep all eligible partners; paralog alternatives are flagged, not ranked away.
        for frag in eligible:
            qs, qe = frag.query_interval
            p, q = resident.project(qs), resident.project(qe)
            if p is None or q is None or p == q:
                stats['unprojectable_fragment'] += 1
                continue
            a0, a1 = sorted((p, q))
            fqs, fqe = resident.query_interval
            left, right = qs-fqs >= config.flank, fqe-qe >= config.flank
            if resident.reverse:
                left, right = right, left
            left = left and resident.anchor_coverage(a0,config.flank)>=.9
            right = right and resident.anchor_coverage(a1,config.flank)>=.9
            orientation = '-' if resident.reverse != frag.reverse else '+'
            bleft, bright = (right, left) if orientation == '-' else (left, right)
            row = dict(read=resident.read, ctgA=resident.contig, magA=panel[resident.contig],
                       a0=a0, a1=a1, ctgB=frag.contig, magB=mag,
                       b0=frag.start, b1=frag.end, orientation=orientation, resident='A',
                       anchor_a0=int(left), anchor_a1=int(right),
                       anchor_b0=int(bleft), anchor_b1=int(bright),
                       resident_start=resident.start, resident_end=resident.end,
                       resident_identity=resident.identity, fragment_identity=frag.identity,
                       ambiguous_fragment=int(len(eligible)>1))
            # Preserve matched bases around both proposed junctions so support can
            # be recomputed at the final clustered boundary, not only at each read's
            # original endpoint. Wider requested flanks/tolerances fail conservatively.
            windows=[(a0-2*config.flank,a0+2*config.flank),
                     (a1-2*config.flank,a1+2*config.flank)]
            blocks=set()
            for s,e in resident.reference_blocks:
                for lo,hi in windows:
                    if min(e,hi)>max(s,lo):blocks.add((max(s,lo),min(e,hi)))
            row['resident_blocks']=';'.join(str(s)+':'+str(e) for s,e in sorted(blocks))
            if (row['magA'],row['ctgA']) > (row['magB'],row['ctgB']):
                for x,y in [('magA','magB'),('ctgA','ctgB'),('a0','b0'),('a1','b1'),
                            ('anchor_a0','anchor_b0'),('anchor_a1','anchor_b1')]:
                    row[x],row[y]=row[y],row[x]
                row['resident']='B'
            out.append(row)
    stats['resident_reads'] += 1
    resident_mag=panel[resident.contig]
    stats['resident_reads::'+resident_mag]+=1
    stats['resident_aligned_bases::'+resident_mag]+=sum(n for n,op in resident.cigar if op in 'M=X')
    stats['evidence_rows'] += len(out)
    return out
