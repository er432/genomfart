import unittest
from genomfart.utils.bigDataFrame import BigDataFrame
from genomfart.data.data_constants import FRAME_TEST_FILE
from Ranger import RangeSet, Range
import numpy as np

debug = False

class bigDataFrameTest(unittest.TestCase):
    """ Tests for bigDataFrame.py """
    @classmethod
    def setUpClass(cls):
        cls.frame = BigDataFrame(FRAME_TEST_FILE)
    def test_getitem(self):
        if debug: print("Testing __getitem__")
        self.assertEqual(self.frame[1,2],27140818)
        self.assertEqual(self.frame[1,'pos'],27140818)
        self.assertEqual(self.frame[6,'pos'],146678420)
        self.assertEqual(self.frame[0,'pos'],146650283)
        self.assertAlmostEqual(self.frame[1,'fval'],36.76025729538268)
        row = self.frame[1]
        self.assertAlmostEqual(row['fval'],36.76025729538268)
    def test_iter(self):
        if debug: print("Testing __iter__")
        row_iter = iter(self.frame)
        row = next(row_iter)
        self.assertEqual(row['pos'],146650283)
        row = next(row_iter)
        self.assertEqual(row['pos'],27140818)
        for row in row_iter:
            try:
                pass
            except StopIteration:
                break
    def test_iterrow(self):
        if debug: print("Testing iterrows")
        # Regular row iterator
        row_iter = self.frame.iterrows(cache=False)
        row = next(row_iter)
        self.assertEqual(row['pos'],146650283)
        row = next(row_iter)
        self.assertEqual(row['pos'],27140818)
        for row in row_iter:
            try:
                pass
            except StopIteration:
                break
        # Column specific row iterator
        row_iter = self.frame.iterrows(colspec=('chr','pos','allele'))
        row = next(row_iter)
        self.assertEqual(len(row),3)
        self.assertEqual(row['pos'],146650283)
        self.assertEqual(row['chr'],10)
        self.assertEqual(row['allele'],'C/A_hm2')
        row = next(row_iter)
        self.assertEqual(len(row),3)
        self.assertEqual(row['pos'],27140818)
        self.assertEqual(row['chr'],10)
        self.assertEqual(row['allele'],'A/T_hm2')
    def test_make_numpy_array(self):
        if debug: print("Testing make_numpy_array")
        arr = self.frame.make_numpy_array(rows=[1,2,4,5],cols=['pos','cm'])
        self.assertEqual(arr.shape,(4,2))
        self.assertEqual(arr[0,0],27140818)
        self.assertAlmostEqual(arr[2,1],43.47540355575215)
        arr2 = self.frame.make_numpy_array(rows=[1,2,4,5],cols=[2,3])
        self.assertTrue(np.allclose(arr,arr2))
        arr2 = self.frame.make_numpy_array(rows=[1,2,4,5],cols=RangeSet(
            ranges=[Range.closed(2,3)]
        ))
        self.assertTrue(np.allclose(arr,arr2))
        arr2 = self.frame.make_numpy_array(rows=[1,2,4,5],cols=RangeSet(
            ranges=[Range.closedOpen(2,3), Range.closedOpen(3,4)]
        ))
        self.assertTrue(np.allclose(arr,arr2))
        arr2 = self.frame.make_numpy_array(rows=RangeSet(
            ranges = [Range.closed(1,2), Range.closed(4,5)]
        ),cols=[2,3])
        self.assertTrue(np.allclose(arr,arr2))
if __name__ == '__main__':
    debug = True
    unittest.main(exit = False)    
