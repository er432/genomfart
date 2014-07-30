import unittest
from Ranger import Range
from genomfart.parsers.gff import gff_parser
from genomfart.data.data_constants import GFF_TEST_FILE

debug = False

class gff_parserTest(unittest.TestCase):
    """ Tests for gff.py """
    @classmethod
    def setUpClass(cls):
        cls.parser = gff_parser(GFF_TEST_FILE)
    def test_get_overlapping_element_ids(self):
        if debug: print("Testing get_overlapping_element_ids")
        elements = self.parser.get_overlapping_element_ids('Pt',100,4000)
        self.assertEqual(len(elements),12)
        self.assertTrue('repeat_region:Pt_320_1262:+' in elements)
        self.assertTrue('repeat_region:Pt_3550_3560:?' in elements)
        self.assertTrue('repeat_region:Pt_3683_3696:?' in elements)
        self.assertTrue('repeat_region:Pt_3764_3775:?' in elements)
        self.assertTrue('exon:Pt_1674_3308:-' in elements)
        self.assertTrue('gene:GRMZM5G836994' in elements)
        self.assertTrue('transcript:GRMZM5G836994_T01' in elements)
        self.assertFalse('CDS:GRMZM5G811749_P01' in elements)
        self.assertTrue('CDS:GRMZM5G836994_P01' in elements)
    def test_get_element_info(self):
        if debug: print("Testing get_element_info")
        gene_info = self.parser.get_element_info('gene:GRMZM5G811749')
        self.assertEqual(gene_info['type'],'gene')
        self.assertEqual(len(gene_info['Ranges']),1)
        self.assertEqual(gene_info['Ranges'][0],Range.closed(3363,5604))
        self.assertEqual(len(gene_info['attributes']),1)
        self.assertEqual(gene_info['attributes'][0]['external_name'],'RPS16')
        self.assertEqual(gene_info['attributes'][0]['biotype'],'protein_coding')
    def test_get_element_ids_of_type(self):
        if debug: print("Testing get_element_ids_of_type")
        iterator = self.parser.get_element_ids_of_type('Pt','gene',start=100,
                                                       end=4000)
        self.assertEqual(next(iterator), 'gene:GRMZM5G836994')
        self.assertEqual(next(iterator), 'gene:GRMZM5G811749')
        with self.assertRaises(StopIteration):
            next(iterator)
if __name__ == "__main__":
    debug = True
    unittest.main(exit = False)
