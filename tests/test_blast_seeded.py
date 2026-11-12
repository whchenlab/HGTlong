import unittest
from hgtlong.blast_seeded import validate, boundary_supported


class FixedBlockTests(unittest.TestCase):
    def candidate(self, **updates):
        return dict(dict(event_id='B1', magA='A', magB='B', ctgA='a', ctgB='b',
                         a0=1000, a1=2000, b0=3000, b1=4000, orientation='+'), **updates)

    def row(self, side, n):
        r=self.candidate();r.pop('event_id')
        r.update(read=side+str(n), resident=side, ambiguous_fragment=0,
                 resident_start=0, resident_end=5000,
                 resident_blocks='0:5000', **{'anchor_'+k:1 for k in ('a0','a1','b0','b1')})
        return r

    def test_two_resident_sides_required(self):
        a=[self.row('A',i) for i in range(3)]
        events,_=validate([self.candidate()],a)
        self.assertEqual(events[0]['status'],'partial')
        events,_=validate([self.candidate()],a+[self.row('B',i) for i in range(3)])
        self.assertEqual(events[0]['status'],'read_supported')

    def test_duplicate_names_not_extra_support(self):
        events,_=validate([self.candidate()],[self.row('A',0)]*4)
        self.assertEqual(events[0]['n_a0'],1)

    def test_multiple_candidate_assignment_excluded(self):
        events,_=validate([self.candidate(),self.candidate(event_id='B2',a0=1010)],
                          [self.row('A',i) for i in range(3)])
        self.assertTrue(all(e['n_a0']==0 and e['unassigned_molecules']==3 for e in events))

    def test_gap_is_not_anchor(self):
        row=self.row('A',0);row['resident_blocks']='0:950;1050:5000'
        self.assertFalse(boundary_supported(row,self.candidate(),'a0'))

    def test_uncensored_fixed_boundary_without_full_read_seed(self):
        rows=[]
        for side in ('A','B'):
            for endpoint in ('0','1'):
                for i in range(3):
                    r=self.row(side,endpoint+str(i))
                    for k in ('a0','a1','b0','b1'):r['anchor_'+k]=int(k.endswith(endpoint))
                    rows.append(r)
        events,_=validate([self.candidate()],rows)
        self.assertEqual(events[0]['status'],'read_supported')

    def test_reverse_and_ambiguous_partner(self):
        r=self.row('B',0);r['orientation']='-';r['ambiguous_fragment']=1
        events,_=validate([self.candidate(orientation='-')],[r])
        self.assertEqual(events[0]['status'],'ambiguous')
        self.assertEqual(events[0]['n_b0'],0)


if __name__=='__main__':unittest.main()
