import unittest
from hgtlong.context_analysis import merge_intervals,shift_intervals,bh_adjust,boundary_clusters,parse_cog_categories
class ContextAnalysis(unittest.TestCase):
    def test_circular_shift_preserves_covered_bases(self):
        original=merge_intervals([(10,30),(20,40),(70,90)])
        shifted=shift_intervals(original,25,100)
        self.assertEqual(sum(b-a for a,b in original),sum(b-a for a,b in shifted))
        self.assertEqual(shifted,[(0,15),(35,65),(95,100)])
    def test_bh_adjustment(self):
        self.assertEqual(bh_adjust([.01,.04,.03]),[.03,.04,.04])
    def test_no_chained_hotspot(self):
        groups=boundary_clusters([dict(position=p) for p in [0,900,1800]],1000)
        self.assertEqual([len(g) for g in groups],[2,1])
    def test_bad_coordinates_rejected(self):
        with self.assertRaises(ValueError):shift_intervals([(1,101)],2,100)
    def test_cog_category_missing_color_and_group_headers(self):
        cats=parse_cog_categories(['2\tCELLULAR PROCESSES\n','J\t1\tFCCCFC\tTranslation\n','X\t2\tMobilome: prophages, transposons\n'])
        self.assertEqual(set(cats),{'J','X'})
        self.assertEqual(cats['X']['name'],'Mobilome: prophages, transposons')
        self.assertTrue(cats['X']['color'].startswith('#'))
