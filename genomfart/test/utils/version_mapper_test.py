import unittest
from genomfart.utils.version_mapper import version_mapper
from genomfart.data.data_constants import VERSION_TEST_FILE

debug = False

class version_mapper_test(unittest.TestCase):
    """ Tests for version_mapper.py """
    @classmethod
    def setUpClass(cls):
        cls.mapper = version_mapper(VERSION_TEST_FILE)
    def test_v1_to_v2_map(self):
        if debug: print("Testing v1_to_v2_map")
        # Cis direction
        chrom,pos,orient = self.mapper.v1_to_v2_map(1, 13857655)
        self.assertEqual((chrom,pos,orient),(1,13856155,1))
        # Trans direction
        chrom,pos,orient = self.mapper.v1_to_v2_map(10,139836890)
        self.assertEquals((chrom,pos,orient),(2,16760272,'-'))
    def test_v2_to_v1_map(self):
        if debug: print("Testing v2_to_v1_map")
        # Cis direction
        chrom,pos,orient = self.mapper.v2_to_v1_map(1,13856155)
        self.assertEqual((chrom,pos,orient),(1, 13857655,1))
        # Trans direction
        chrom,pos,orient = self.mapper.v2_to_v1_map(2,16760272)
        self.assertEquals((chrom,pos,orient),(10,139836890,'-'))
    def test_v1_to_v2_seg_map(self):
        if debug: print("Testing v1_to_v2_seg_map")
        # Cis direction    
        seg_map = self.mapper.v1_to_v2_seg_map(1, 13857655, 13857659)
        self.assertEqual(list(seg_map.values())[0],
            (1, 13856155, 13856159, 1))
        # Trans direction
        seg_map = self.mapper.v1_to_v2_seg_map(10,139836890,139836895)
        self.assertEqual(list(seg_map.values())[0],
            (2,16760272,16760267,'-')
        )
    def test_v2_to_v1_seg_map(self):
        if debug: print("Testing v2_to_v1_seg_map")
        # Cis direction    
        # seg_map = self.mapper.v1_to_v2_seg_map(1, 13857655, 13857659)
        # self.assertEqual(list(seg_map.values())[0],
        #     (1, 13856155, 13856159, 1))
        seg_map = self.mapper.v2_to_v1_seg_map(1, 13856155, 13856159)
        self.assertEqual(list(seg_map.values())[0],
                         (1, 13857655, 13857659, 1))
        # Trans direction
        # seg_map = self.mapper.v1_to_v2_seg_map(10,139836890,139836895)
        # self.assertEqual(list(seg_map.values())[0],
        #     (2,16760272,16760267,'-')
        # )
        seg_map = self.mapper.v2_to_v1_seg_map(2,16760267,16760272)
        self.assertEqual(list(seg_map.values())[0],
                         (10,139836895,139836890,'-'))
if __name__ == '__main__':
    debug = True
    unittest.main(exit = False)
