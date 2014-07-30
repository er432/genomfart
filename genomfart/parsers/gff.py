import networkx as nx
from Ranger import RangeBucketMap, Range

class gff_parser(object):
    """ Class used to parse and analyze GFF (version 3) files.
    The class represents any hierarchical structure as a directed graph
    and puts the individual pieces into a RangeBucketMap

    All coordinates are 1-based
    """
    def __init__(self, gff_file):
        """ Instantiates the gff file

        Parameters
        ----------
        gff_file : str
            The filename of a gff file

        Raises
        ------
        IOError
            If the file isn't correctly formatted
        """
        ## Set up a directed graph where nodes are ids and a directed edge
        # is placed going from parent to child. Each node has a "Ranges" attribute
        # that lists the Ranges corresponding to the ID, a "seqid" attribute
        # that gives the coordinate system where the node is located, a
        # "type" attribute that gives the element type, a "strand" attribute
        # that gives the element strand (not always applicable), and an "attributes"
        # attribute that is a list of dictionaries of key-> val for any other attributes
        self.graph = nx.DiGraph()
        ## Dictionary of seq_id -> RangeBucketMap for each coordinate system. Each Range
        # in a RangeBucketMap maps to an ID that corresponds to a node in the graph
        self.bucketmaps = {}
        ## Parse the file
        with open(gff_file) as gff_handle:
            for line in gff_handle:
                if line.startswith('#'): continue
                elif len(line) < 2: continue
                line = line.strip().split('\t')
                if len(line) != 9:
                    raise IOError("Line(s) do not conform to gff v. 3 format")
                # Make a new RangeBucketMap for the seqid if necessary
                seqid = line[0]
                if seqid not in self.bucketmaps:
                    self.bucketmaps[seqid] = RangeBucketMap()
                # Parse the attributes into a dictionary of key->val
                attr_dict = dict((k,v) for k,v in map(lambda x: x.split('='),
                                                      line[8].split(';')))
                # Check if the element has an ID. If not, give it one
                if 'ID' in attr_dict:
                    element_id = attr_dict['ID']
                else:
                    element_id = '%s:%s_%s_%s:%s' % (line[2],line[0],line[3],line[4],
                                                     line[6])
                # Make the Range for this element
                element_range = Range.closed(int(line[3]),int(line[4]))
                # Put element in RangeBucketMap
                self.bucketmaps[seqid].put(element_range, element_id)
                # Put node in the graph
                if element_id not in self.graph:
                    self.graph.add_node(element_id)
                    self.graph.node[element_id]['Ranges'] = []
                    self.graph.node[element_id]['attributes'] = []
                # Add/update parameters in graph
                self.graph.node[element_id]['seqid'] = line[0]
                self.graph.node[element_id]['Ranges'].append(element_range)
                self.graph.node[element_id]['attributes'].append(attr_dict)
                self.graph.node[element_id]['type'] = line[2]
                self.graph.node[element_id]['strand'] = line[6]
                # Make edges from parent to child if necessary
                if 'Parent' in attr_dict:
                    parents = attr_dict['Parent'].split(',')
                    for parent in parents:
                        self.graph.add_edge(parent, element_id)
    def get_overlapping_element_ids(self, seqid, start, end):
        """ Gets the ids for any elements that overlap a given range

        Parameters
        ----------
        seqid : str
            The name of the coordinate system to check
        start : int
            The start of the range to check (inclusive, 1-based)
        end : int
            The end of the range to check (inclusive, 1-based)

        Raises
        ------
        KeyError
            If the seqid is not present
        
        Returns
        -------
        Set of ids for elements overlapping the range
        """
        checkRange = Range.closed(start, end)
        return self.bucketmaps[seqid].get(checkRange)
    def get_element_info(self, element_id):
        """ Gets information on a particular element

        Parameters
        ----------
        element_id : str
            The id of the element

        Returns
        -------
        Dictionary of {'seqid' -> seqid, type'->type, 'strand'->'strand', 'Ranges'->[Ranges],
        'attributes'->[attribute_dicts]}
        """
        return {'type':self.graph.node[element_id]['type'],
                'strand':self.graph.node[element_id]['strand'],
                'seqid':self.graph.node[element_id]['seqid'],
                'Ranges':self.graph.node[element_id]['Ranges'],
                'attributes':self.graph.node[element_id]['attributes']}
    def get_element_ids_of_type(self, seqid, element_type, start = None, end = None):
        """ Gets element ids of some type along a coordinate system

        Parameters
        ----------
        seqid : str
            The name of the coordinate system to check
        element_type : str
            The type of the elements you want (e.g. 'gene' or 'mRNA')
        start : int, optional
            The start point for getting the elements (inclusive, 1-based)
        end : int, optional
            The end point for getting the elements (inclusive, 1-based)

        Raises
        ------
        KeyError
            If the seqid isn't included
        
        Returns
        -------
        Generator of element_ids
        """
        added = set()
        iterator = self.bucketmaps[seqid].iteritems(start=start,end=end)
        for theRange, element_id in iterator:
            if element_id in added:
                continue
            elif self.graph.node[element_id]['type'] == element_type:
                yield element_id
                added.add(element_id)
    def get_element_children_ids(self, element_id):
        """ Gets the ids of the children of an element

        Parameters
        ----------
        element_id : str
            The element for which you want the children

        Raises
        ------
        KeyError
            If the element is not present

        Returns
        -------
        List of the children IDs of the element
        """
        return self.graph.successors(element_id)
    def get_element_parent_ids(self, element_id):
        """ Gets the ids of the parents of an element

        Parameters
        ----------
        element_id : str
            The element for which you want the parents

        Raises
        ------
        KeyError
            If the element is not present

        Returns
        -------
        List of the parent IDs of the element
        """
        return self.graph.predecessors(element_id)
