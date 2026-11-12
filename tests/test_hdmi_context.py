import unittest
from hgtlong.hdmi_context import compare_geometry,overlap_category
class HDMIContext(unittest.TestCase):
    def setUp(self):self.call=dict(ctgA='a',ctgB='b',a0=100,a1=200,b0=300,b1=400,orientation='-')
    def test_swapped_references_are_exact(self):
        other=dict(ctgA='b',ctgB='a',a0=300,a1=400,b0=100,b1=200,orientation='-')
        self.assertTrue(compare_geometry(self.call,other)['exact'])
    def test_one_sided_and_opposite_orientation_do_not_match(self):
        self.assertFalse(compare_geometry(self.call,dict(self.call,b0=400,b1=500))['bilateral'])
        self.assertFalse(compare_geometry(self.call,dict(self.call,orientation='+'))['bilateral'])
    def test_final_overlap_precedes_raw_exact(self):
        raw=dict(compare_geometry(self.call,self.call),final_retained=False)
        final=dict(compare_geometry(self.call,dict(self.call,a0=150)),final_retained=True)
        self.assertEqual(overlap_category([raw,final]),'partial_overlap_with_HDMI_final')
        final['exact']=True
        self.assertEqual(overlap_category([raw,final]),'identical_to_HDMI_final')
    def test_raw_only_and_no_overlap(self):
        self.assertEqual(overlap_category([dict(compare_geometry(self.call,self.call),final_retained=False)]),'overlap_HDMI_raw_not_final')
        self.assertEqual(overlap_category([]),'no_bilateral_HDMI_raw_overlap')
