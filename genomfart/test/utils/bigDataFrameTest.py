import unittest
from genomfart.utils.bigDataFrame import BigDataFrame
from genomfart.data.data_constants import FRAME_TEST_FILE

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
if __name__ == '__main__':
    debug = False
    unittest.main(exit = False)    
