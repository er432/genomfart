import unittest
from genomfart.parsers.gff import gff_parser
from genomfart.data.data_constants import GFF_TEST_FILE

class gff_parserTest(unittest.TestCase):
    """ Tests for gff.py """
    @classmethod
    def setUpClass(cls):
        cls.parser = gff_parser(GFF_TEST_FILE)
