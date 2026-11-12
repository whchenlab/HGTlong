import unittest
from collections import Counter
from types import SimpleNamespace
from hgtlong.core import Alignment, classify, union_length
from hgtlong.events import call_events

C=SimpleNamespace(full_cov=.9,full_identity=.99,frag_cov=.75,frag_identity=.95,
                  min_mapq=20,min_block=200,max_indel=500,flank=100,
                  min_support=3,junction_tolerance=50,min_length_ratio=.9)


def aln(name='r',ctg='a',pos=1001,cigar='2000M',flag=0,nm=0,mapq=60):
    return Alignment.from_sam(f'{name}\t{flag}\t{ctg}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t*\t*\tNM:i:{nm}')


def evidence(name='r',reverse=False):
    return classify([aln(name),aln(name,'b',5001,'500S1000M500S',272 if reverse else 256)],{'a':'A','b':'B'},C,Counter())[0]


class Geometry(unittest.TestCase):
    def test_indels_identity(self):
        r=aln(cigar='500M10I500M10D990M',nm=20)
        self.assertAlmostEqual(r.identity,1-20/2010)
        self.assertEqual(r.query_interval,(0,2000))

    def test_hardclip(self):
        r=aln(cigar='300H500M200S',flag=16)
        self.assertEqual(r.length,1000)
        self.assertEqual(r.query_interval,(200,700))
        self.assertEqual(r.project(200),1500)
        self.assertEqual(r.project(700),1000)

    def test_insertion_projection(self):
        r=aln(cigar='500M10I1490M')
        self.assertIsNone(r.project(505))
        self.assertEqual(r.project(600),1590)

    def test_deletion_cannot_support_anchor(self):
        r=aln(cigar='500M100D1500M')
        self.assertEqual(r.anchor_coverage(1550,100),.5)

    def test_real_anchor(self):
        self.assertEqual(aln().anchor_coverage(1500,100),1)

    def test_union(self):
        self.assertEqual(union_length([(0,100),(10,50),(90,110)]),110)

    def test_pair(self):
        r=evidence()
        self.assertEqual((r['a0'],r['a1'],r['b0'],r['b1']),(1500,2500,5000,6000))

    def test_reverse_pair(self):
        self.assertEqual(evidence(reverse=True)['orientation'],'-')

    def test_no_quality_borrowing(self):
        rows=classify([aln(nm=100),aln(ctg='a',cigar='1800S200M'),aln(ctg='b',cigar='500S1000M500S',flag=256)],{'a':'A','b':'B'},C,Counter())
        self.assertEqual(rows,[])

    def test_duplicate_records(self):
        r=aln(); f=aln(ctg='b',cigar='500S1000M500S',flag=256)
        self.assertEqual(len(classify([r,r,f,f],{'a':'A','b':'B'},C,Counter())),1)

    def test_duplicate_placement_flags_not_paralogs(self):
        r=aln();f=aln(ctg='b',cigar='500S1000M500S',flag=256)
        g=aln(ctg='b',cigar='500S1000M500S',flag=2048)
        rows=classify([r,f,g],{'a':'A','b':'B'},C,Counter())
        self.assertEqual(len(rows),1)
        self.assertEqual(rows[0]['ambiguous_fragment'],0)

    def test_alignment_order_invariance(self):
        rows=[aln(),aln(ctg='b',pos=5001,cigar='500S1000M500S',flag=256),
              aln(ctg='b',pos=9001,cigar='500S1000M500S',flag=256)]
        self.assertEqual(classify(rows,{'a':'A','b':'B'},C,Counter()),
                         classify(reversed(rows),{'a':'A','b':'B'},C,Counter()))

    def test_alternative_full_locus(self):
        self.assertEqual(classify([aln(),aln(ctg='b',mapq=0,flag=256)],{'a':'A','b':'B'},C,Counter()),[])

    def test_unknown_mapq(self):
        self.assertEqual(classify([aln(mapq=255)],{'a':'A'},C,Counter()),[])

    def test_resident_depth_counts_include_reads_without_fragments(self):
        stats=Counter()
        self.assertEqual(classify([aln()],{'a':'A'},C,stats),[])
        self.assertEqual(stats['resident_reads::A'],1)
        self.assertEqual(stats['resident_aligned_bases::A'],2000)

    def test_missing_nm(self):
        self.assertEqual(aln(nm=-1).identity,0)

    def test_four_boundaries(self):
        rows=[evidence('a'+str(i)) for i in range(3)]
        for i in range(3):
            r=evidence('b'+str(i)); r.update(resident='B',resident_start=4500,resident_end=6500,
                                          resident_blocks='4500:6500')
            rows.append(r)
        events,_=call_events(rows,C)
        self.assertEqual(len(events),1)
        self.assertEqual(events[0]['status'],'read_supported')
        self.assertEqual(events[0]['n_b0'],3)

    def test_canonical_locus_does_not_borrow_unassigned_support(self):
        rows=[evidence('a'+str(i)) for i in range(3)]
        for i in range(3):
            row=evidence('neighbor'+str(i))
            for k in ('a0','a1','b0','b1'):row[k]+=60
            rows.append(row)
        for i in range(10):
            row=evidence('uncertain'+str(i))
            row.update(resident='B',resident_start=4500,resident_end=6500,
                       resident_blocks='4500:6500')
            for k in ('a0','a1','b0','b1'):row['anchor_'+k]=0
            rows.append(row)
        events,members=call_events(rows,C)
        self.assertEqual(len(events),2)
        self.assertEqual(len({e['event_id'] for e in events}),2)
        event=next(e for e in events if e['a0']==1500)
        self.assertEqual(event['unassigned_molecules'],10)
        self.assertEqual(event['n_b0'],0)
        self.assertEqual(event['status'],'partial')
        self.assertEqual(sum(r['assignment_ambiguous'] for r in members),10)

    def test_no_alignment_count_inflation(self):
        events,_=call_events([evidence()]*20,C)
        self.assertEqual(events,[])

    def test_one_sided_is_partial(self):
        events,_=call_events([evidence(str(i)) for i in range(3)],C)
        self.assertEqual(events[0]['status'],'partial')

    def test_final_boundary_checks_aligned_bases(self):
        rows=[dict(evidence(str(i)),resident_blocks='1000:1480;1520:3000') for i in range(3)]
        events,_=call_events(rows,C)
        self.assertEqual(events[0]['n_a0'],0)
        self.assertEqual(events[0]['n_a1'],3)

    def test_old_evidence_fails_explicitly(self):
        row=evidence();del row['resident_blocks']
        with self.assertRaisesRegex(ValueError,'re-extract'):
            call_events([row],C)

    def test_paralog_not_supported(self):
        rows=[dict(evidence(str(i)),ambiguous_fragment=1) for i in range(3)]
        self.assertEqual(call_events(rows,C)[0][0]['status'],'ambiguous')

    def test_order_stability(self):
        rows=[evidence(str(i)) for i in range(3)]
        self.assertEqual(call_events(rows,C)[0],call_events(reversed(rows),C)[0])


if __name__=='__main__': unittest.main()
