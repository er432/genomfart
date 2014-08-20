import unittest
from genomfart.parsers.vcf import VCF_parser
from genomfart.data.data_constants import VCF_TEST_FILE

debug = False

class vcfTest(unittest.TestCase):
    """ Unit tests for vcf.py """
    @classmethod
    def setUpClass(cls):
        cls.parser = VCF_parser(VCF_TEST_FILE)
    def test_parse_geno_depths(self):
        if debug: print("Testing parse_geno_depths")
        depth_parser = self.parser.parse_geno_depths()
        chrom,pos,allels,depths = next(depth_parser)
        self.assertEqual(sum(depths['COSTICH_2014']),1)
        self.assertEquals(sum(depths['FL_05_15_1']),0)
        chrom,pos,allels,depths = next(depth_parser)
        chrom,pos,allels,depths = next(depth_parser)
        self.assertEqual(depths['COSTICH_2014'],(2,1))
        self.assertEquals(sum(depths['FL_05_15_1']),0)
    def test_parse_select_geno_depths(self):
        if debug: print("Testing parse_select_geno_depths")
        depth_parser = self.parser.parse_select_geno_depths(['COSTICH_2014','FL_34798'],
                                                            info_dict=True, start=100,
                                                            end=130)
        chrom,pos,alleles,depths,info = next(depth_parser)
        self.assertEqual(len(depths), 2)
        self.assertEqual(pos, 105)
        self.assertEqual(depths['COSTICH_2014'],(9,2))
        self.assertEqual(depths['FL_34798'],(1,0))
        chrom,pos,alleles,depths,info = next(depth_parser)
        self.assertEqual(pos, 124)
        chrom,pos,alleles,depths,info = next(depth_parser)
        with self.assertRaises(StopIteration):
            next(depth_parser)
    def test_parse_select_geno_generic(self):
        if debug: print("Testing parse_select_geno_generic")
        geno_parser = self.parser.parse_select_geno_generic(['COSTICH_2014','FL_05_15_1','FL_34798'],
                                                            info_dict=True, start=100,
                                                            end=130)
        chrom,pos,alleles,genos,info = next(geno_parser)
        self.assertEqual(genos['COSTICH_2014']['GT'], '0/1')
        self.assertEqual(genos['COSTICH_2014']['AD'], '9,2')
        self.assertEqual(genos['FL_05_15_1']['GT'], './.')
        with self.assertRaises(KeyError):
            genos['FL_05_15_1']['AD']
        chrom,pos,alleles,gemps,info = next(geno_parser)
        self.assertEqual(len(genos), 3)
        self.assertEqual(pos, 124)
        chrom,pos,alleles,genos,info = next(geno_parser)
        with self.assertRaises(StopIteration):
            next(geno_parser)
    def test_parse_site_infos(self):
        if debug: print("Testing parse_site_infos")
        info_parser = self.parser.parse_site_infos()
        chrom,pos,alleles,info = next(info_parser)
        self.assertEqual(pos, 13)
        self.assertEqual(info['GN'],(2,1,0))
        self.assertEqual(info['AQ'],(36,33))
    def test_get_affected_ref_bases(self):
        if debug: print("Testing get_affected_ref_bases")
        self.assertEqual(VCF_parser.get_affected_ref_bases(20,'C','T'),set([20]))
        self.assertEqual(VCF_parser.get_affected_ref_bases(20,'C','CTAG'), set([]))
        self.assertEqual(VCF_parser.get_affected_ref_bases(20,'TCG','T'), set([21,22]))
        self.assertEqual(VCF_parser.get_affected_ref_bases(20,'TCGCG','TCG'),set([24,23]))
        self.assertEqual(VCF_parser.get_affected_ref_bases(20,'TCGCG','TCGCGCG'),set([]))
if __name__ == "__main__":
    debug = True
    unittest.main(exit = False)
