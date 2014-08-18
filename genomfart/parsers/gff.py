from genomfart.utils.genomeAnnotationGraph import genomeAnnotationGraph
import gzip

class gff_parser(genomeAnnotationGraph):
    """ Class used to parse and analyze GFF (version 3) files.
    The class represents any hierarchical structure as a directed graph
    and puts the individual pieces into a RangeBucketMap

    All coordinates are 1-based

    Examples
    --------

    >>> from genomfart.parsers.gff import gff_parser
    >>> from genomfart.data.data_constants import GFF_TEST_FILE
    >>> parser = gff_parser(GFF_TEST_FILE)
    >>> parser.get_overlapping_element_ids('Pt',100,4000)
    set(['three_prime_UTR:Pt_3363_4490:-', 'exon:Pt_1674_3308:-',
    'CDS:GRMZM5G836994_P01', 'exon:Pt_3363_4708:-',
    'transcript:GRMZM5G811749_T01', 'transcript:GRMZM5G836994_T01',
    'repeat_region:Pt_3550_3560:?', 'repeat_region:Pt_3683_3696:?',
    'gene:GRMZM5G811749', 'gene:GRMZM5G836994', 'repeat_region:Pt_320_1262:+',
    'repeat_region:Pt_3764_3775:?'])
    >>> gene_info = parser.get_element_info('gene:GRMZM5G811749')
    >>> gene_info
    {'Ranges': [[3363 , 5604]], 'attributes': [{'logic_name': 'genebuilder',
    'external_name': 'RPS16',
    'description': '30S ribosomal protein S16%2C chloroplastic  [Source:UniProtKB/Swiss-Prot%3BAcc:P27723]', 'ID': 'gene:GRMZM5G811749',
    'biotype': 'protein_coding'}], 'seqid': 'Pt', 'type': 'gene', 'strand': '-'}
    >>> [x for x in parser.get_element_ids_of_type('Pt','gene',start=100,end=4000)]
    ['gene:GRMZM5G836994', 'gene:GRMZM5G811749']
    """
    def __init__(self, gff_file, exclude_types = None):
        """ Instantiates the gff file

        Parameters
        ----------
        gff_file : str
            The filename of a gff file
        exclude_types : set
            The names of types (e.g. 'repeat') that you don't want to store

        Raises
        ------
        IOError
            If the file isn't correctly formatted
        """
        super(gff_parser, self).__init__()
        if exclude_types is None: exclude_types = set()
        if gff_file.endswith('.gz'):
            gff_handle = gzip.open(gff_file)
        else:
            gff_handle = open(gff_file)
        ## Parse the file
        with gff_handle:
            for line in gff_handle:
                if line.startswith('#'): continue
                elif len(line) < 2: continue
                line = line.strip().split('\t')
                if len(line) != 9:
                    raise IOError("Line(s) do not conform to gff v. 3 format")
                if line[2] in exclude_types: continue
                # Make a new RangeBucketMap for the seqid if necessary
                seqid = line[0]
                # Parse the attributes into a dictionary of key->val
                attr_dict = dict((k,v) for k,v in map(lambda x: x.split('='),
                                                      line[8].split(';')))
                # Check if the element has an ID. If not, give it one
                if 'ID' in attr_dict:
                    element_id = attr_dict['ID']
                else:
                    element_id = '%s:%s_%s_%s:%s' % (line[2],line[0],line[3],line[4],
                                                     line[6])
                parents = set()
                if 'Parent' in attr_dict:
                    parents = parents.union(attr_dict['Parent'].split(','))
                self.add_annotation(element_id, line[0], int(line[3]), int(line[4]),
                                    line[2], strand=line[6], parents=parents,
                                    **attr_dict)
