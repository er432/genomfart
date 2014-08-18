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
        self.assertEqual(gene_info['seqid'], 'Pt')
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
    def test_get_element_children_ids(self):
        if debug: print("Testing get_element_children_ids")
        children = self.parser.get_element_children_ids('gene:GRMZM5G811749')
        self.assertEqual(len(children),1)
        self.assertEqual(children[0],'transcript:GRMZM5G811749_T01')
        children = self.parser.get_element_children_ids('transcript:GRMZM5G811749_T01')
        self.assertEqual(len(children), 3)
        self.assertFalse('gene:GRMZM5G811749' in children)
        self.assertTrue('CDS:GRMZM5G811749_P01' in children)
    def test_get_element_parent_ids(self):
        if debug: print("Testing get_element_parent_ids")
        parents = self.parser.get_element_parent_ids('gene:GRMZM5G811749')
        self.assertEqual(len(parents),0)
        parents = self.parser.get_element_parent_ids('transcript:GRMZM5G811749_T01')
        self.assertEqual(len(parents),1)
        self.assertEqual(parents[0],'gene:GRMZM5G811749')
        parents = self.parser.get_element_parent_ids('CDS:GRMZM5G811749_P01')
        self.assertEqual(len(parents),1)
        self.assertEqual(parents[0],'transcript:GRMZM5G811749_T01')
    def test_get_closest_element_id(self):
        if debug: print("Testing get_closest_element_id")
        closest_ids = self.parser.get_closest_element_id('Pt',2000,2200)
        self.assertEqual(closest_ids, set(['gene:GRMZM5G836994','transcript:GRMZM5G836994_T01',
                                           'CDS:GRMZM5G836994_P01','exon:Pt_1674_3308:-']))
        closest_ids = self.parser.get_closest_element_id('Pt',15838,15840)
        self.assertEqual(closest_ids, set(['repeat_region:Pt_15787_15836:?']))
        closest_ids = self.parser.get_closest_element_id('Pt',15880,15882)
        self.assertEqual(closest_ids, set(['repeat_region:Pt_15885_15961:+']))
    def test_get_closest_element_id_of_type(self):
        if debug: print("Testing get_closest_element_id")
        closest_ids = self.parser.get_closest_element_id_of_type('Pt',2000,2200,'gene')
        self.assertEqual(closest_ids, set(['gene:GRMZM5G836994']))
        closest_ids = self.parser.get_closest_element_id_of_type('Pt',15838,15840,'repeat_region')
        self.assertEqual(closest_ids, set(['repeat_region:Pt_15787_15836:?']))
        closest_ids = self.parser.get_closest_element_id_of_type('Pt',15880,15882,'gene')
        self.assertEqual(closest_ids, set(['gene:GRMZM5G813608']))        
if __name__ == "__main__":
    debug = True
    unittest.main(exit = False)
