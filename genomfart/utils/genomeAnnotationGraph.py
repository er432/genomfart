import networkx as nx
from Ranger import RangeBucketMap, Range

class genomeAnnotationGraph(object):
    """ Representation of a genome, in which annotations in the genome are Ranges
    in a RangeBucketMap, and annotations can be hierarchical (e.g. transcripts are the children
    of genes). The hierarchy is stored as a directed graph
    """
    def __init__(self):
        """ Instantiates the genomeAnnotationGraph
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
    def add_annotation(self, element_id, seqid, start, end, element_type,
                       strand='?', parents = None, children = None,  **attr):
        """ Adds an annotation to the genome

        Note that you can add more than one annotation under an element id,
        so this can be preexisting. If it is preexisting, the Range and attribute dictionaries
        will be appended

        Parameters
        ----------
        element_id : hashable
            The id of the element you're adding. 
        seqid : hashable
            The name of the coordinate system (e.g. the chromosome)
        start : int
            The start of the element (inclusive)
        end : int
            The end of the element (inclusive)
        element_type : str
            The type of the element (e.g. 'gene' or 'mRNA')
        strand : str, optional
            The element strand ('+','-','?')
        parents : iterable of strings, optional
            The parent element_ids of this element_id
        children : iterable of strings, optional
            The child element_ids of this element_id
        attr : Further keyword arguments that will be stored as attributes
               of the element
        """
        if not parents: parents = set()
        if not children: children = set()
        # Put the seqid in the bucketmaps dictionary if necessary
        if seqid not in self.bucketmaps:
            self.bucketmaps[seqid] = RangeBucketMap()
        # Make the Range for this element
        element_range = Range.closed(start, end)
        # Put element in RangeBucketMap
        self.bucketmaps[seqid].put(element_range, element_id)
        # Put node in graph if necessary
        if element_id not in self.graph:
            self.graph.add_node(element_id)
            self.graph.node[element_id]['Ranges'] = []
            self.graph.node[element_id]['attributes'] = []
        # Add/update parameters in graph
        self.graph.node[element_id]['seqid'] = seqid
        self.graph.node[element_id]['Ranges'].append(element_range)
        self.graph.node[element_id]['attributes'].append(attr)
        self.graph.node[element_id]['type'] = element_type
        self.graph.node[element_id]['strand'] = strand
        # Make edges if necessary
        for parent in parents:
            self.graph.add_edge(parent, element_id)
        for child in children:
            self.graph.add_edge(element_id, child)
    def add_node_annotations(self, element_id, **annots):
        """ Adds annotations to all annotation dictionaries under a node.
        Note that if the key for an annotation is the same as a previously
        existing one, that annotation will be overwritten

        Parameters
        ----------
        element_id : hashable
            The id of the element you're adding.
        annots : keyword arguments
            The elements you want to add

        Examples
        --------

        >>> my_genome.add_node_annotations('myNode1', expression1 = 0.2, expression2 = 10.)
        """
        for attr_dict in self.graph.node[element_id]['attributes']:
            attr_dict.update(annots)
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
        if seqid not in self.bucketmaps:
            raise KeyError("%s not present" % seqid)
        checkRange = Range.closed(start, end)
        try:
            return self.bucketmaps[seqid].get(checkRange)
        except KeyError:
            return set()
    def get_overlapping_element_ids_of_type(self, seqid, start, end, element_type):
        """ Gets the ids for any elements that overlap a given range and are
        of a given type

        Parameters
        ----------
        seqid : str
            The name of the coordinate system to check
        start : int
            The start of the range to check (inclusive, 1-based)
        end : int
            The end of the range to check (inclusive, 1-based)
        element_type : str
            The element type that you want (e.g. 'gene', 'mRNA')

        Raises
        ------
        KeyError
            If the seqid is not present

        Returns
        -------
        Set of ids for elements overlapping the range that are of the
        element type
        """
        all_overlaps = self.get_overlapping_element_ids(seqid, start, end)
        return set([x for x in all_overlaps if self.graph.node[x]['type'] == element_type])
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
    def get_closest_element_id(self, seqid, rangeStart, rangeEnd, radius = 10000):
        """ Gets the element id(s) of the whatever element is closest to a range

        Parameters
        ----------
        seqid : str
            The name of the coordinate system to check
        rangeStart : int
            The position (inclusive) beginning the range for which you want
            the closest element
        rangeEnd : int
            The position (inclusive) ending the range for which you want
            the closest element
        radius : int
            How far on either side of the search range you want to search
            for the closest element

        Returns
        -------
        Set that will contain all element ids overlapping the search range
        if there are any. Otherwise, it will contain the element(s) with the
        shortest distance to the search range (only more than 1 if some elements
        are equidistant)
        """
        checkRange = Range.closed(rangeStart, rangeEnd)
        # Go through items overlapping the Range
        iterator = self.bucketmaps[seqid].iteritems(start=(rangeStart-radius),
                                                    end=(rangeEnd+radius))
        # Store the distance and element id of closest element that is not
        # overlapping the search range and is to the left of the search range
        less_than_closest_dist = radius
        less_than_closest_element = set()
        # Store the distance and element id of the closest element that is
        # not overlapping the search range and is to the right of the search range
        greater_than_closest_dist = radius
        greater_than_closest_element = set()
        # Store any overlapping elements
        overlapping = set()
        for theRange, element_id in iterator:
            # Close the Range if necessary
            if not all((theRange.isLowerBoundClosed(), theRange.isUpperBoundClosed())):
                newLowerBd = theRange.lowerEndpoint() if theRange.isLowerBoundClosed() else \
                  (theRange.lowerEndpoint()+1)
                newUpperBd = theRange.upperEndpoint() if theRange.isUpperBoundClosed() else \
                  (theRange.upperEndpoint()-1)
                theRange = Range.closed(newLowerBd, newUpperBd)
            if checkRange.isConnected(theRange):
                # Add element to overlapping set if overlapping
                overlapping.add(element_id)
            elif theRange.lowerEndpoint() < checkRange.lowerEndpoint():
                # If the element is to the left of the search range
                range_dist = checkRange.getDistanceFromRange(theRange)
                if range_dist < less_than_closest_dist:
                    # Add the element to the less than set after clearing current
                    # contents 
                    less_than_closest_element.clear()
                    less_than_closest_element.add(element_id)
                    less_than_closest_dist = range_dist
                elif range_dist == less_than_closest_dist:
                    # Add the element to the current set
                    less_than_closest_element.add(element_id)
            elif theRange.lowerEndpoint() > checkRange.upperEndpoint():
                # If the element is to the right of the search range
                if len(overlapping) > 0:
                    # Return the overlapping set if there are any overlapping elements
                    return overlapping
                range_dist = checkRange.getDistanceFromRange(theRange)
                if range_dist > less_than_closest_dist:
                    # Return the closest elements if the current range distance is
                    # greater than the left distances
                    if less_than_closest_dist < greater_than_closest_dist:
                        return less_than_closest_element
                    elif greater_than_closest_dist < less_than_closest_dist:
                        return greater_than_closest_element
                    else:
                        return less_than_closest_element.union(greater_than_closest_element)
                elif range_dist < greater_than_closest_dist:
                    # Add the element to the greater than set after clearing
                    # current contents
                    greater_than_closest_element.clear()
                    greater_than_closest_element.add(element_id)
                    greater_than_closest_dist = range_dist
                elif range_dist == greater_than_closest_dist:
                    # Add the element to the current set
                    greater_than_closest_element.add(element_id)
                elif range_dist > greater_than_closest_dist:
                    # Return the closest elements if the current range distance
                    # is greater than the minimum right distance
                    if less_than_closest_dist < greater_than_closest_dist:
                        return less_than_closest_element
                    elif greater_than_closest_dist < less_than_closest_dist:
                        return greater_than_closest_element
                    else:
                        return less_than_closest_element.union(greater_than_closest_element)
        # Return the closest if finished going through all elements
        if len(overlapping) > 0:
            return overlapping
        elif less_than_closest_dist < greater_than_closest_dist:
            return less_than_closest_element
        elif greater_than_closest_dist < less_than_closest_dist:
            return greater_than_closest_element
        elif less_than_closest_dist == greater_than_closest_dist:
            return less_than_closest_element.union(greater_than_closest_element)

    def get_closest_element_id_of_type(self, seqid, rangeStart, rangeEnd, element_type,
                                       radius = 10000):
        """ Gets the element id(s) of the whatever element is closest to a range

        Parameters
        ----------
        seqid : str
            The name of the coordinate system to check
        rangeStart : int
            The position (inclusive) beginning the range for which you want
            the closest element
        rangeEnd : int
            The position (inclusive) ending the range for which you want
            the closest element
        radius : int
            How far on either side of the search range you want to search
            for the closest element
        element_type : str
            The type of the elements you want (e.g. 'gene' or 'mRNA')

        Returns
        -------
        Set that will contain all element ids overlapping the search range
        if there are any. Otherwise, it will contain the element(s) with the
        shortest distance to the search range (only more than 1 if some elements
        are equidistant)
        """
        checkRange = Range.closed(rangeStart, rangeEnd)
        # Go through items overlapping the Range
        iterator = self.bucketmaps[seqid].iteritems(start=(rangeStart-radius),
                                                    end=(rangeEnd+radius))
        # Store the distance and element id of closest element that is not
        # overlapping the search range and is to the left of the search range
        less_than_closest_dist = radius
        less_than_closest_element = set()
        # Store the distance and element id of the closest element that is
        # not overlapping the search range and is to the right of the search range
        greater_than_closest_dist = radius
        greater_than_closest_element = set()
        # Store any overlapping elements
        overlapping = set()
        for theRange, element_id in iterator:
            # Skip if element_id not of the right type
            if self.graph.node[element_id]['type'] != element_type: continue
            # Close the Range if necessary
            if not all((theRange.isLowerBoundClosed(), theRange.isUpperBoundClosed())):
                newLowerBd = theRange.lowerEndpoint() if theRange.isLowerBoundClosed() else \
                  (theRange.lowerEndpoint()+1)
                newUpperBd = theRange.upperEndpoint() if theRange.isUpperBoundClosed() else \
                  (theRange.upperEndpoint()-1)
                theRange = Range.closed(newLowerBd, newUpperBd)
            if checkRange.isConnected(theRange):
                # Add element to overlapping set if overlapping
                overlapping.add(element_id)
            elif theRange.lowerEndpoint() < checkRange.lowerEndpoint():
                # If the element is to the left of the search range
                range_dist = checkRange.getDistanceFromRange(theRange)
                if range_dist < less_than_closest_dist:
                    # Add the element to the less than set after clearing current
                    # contents 
                    less_than_closest_element.clear()
                    less_than_closest_element.add(element_id)
                    less_than_closest_dist = range_dist
                elif range_dist == less_than_closest_dist:
                    # Add the element to the current set
                    less_than_closest_element.add(element_id)
            elif theRange.lowerEndpoint() > checkRange.upperEndpoint():
                # If the element is to the right of the search range
                if len(overlapping) > 0:
                    # Return the overlapping set if there are any overlapping elements
                    return overlapping
                range_dist = checkRange.getDistanceFromRange(theRange)
                if range_dist > less_than_closest_dist:
                    # Return the closest elements if the current range distance is
                    # greater than the left distances
                    if less_than_closest_dist < greater_than_closest_dist:
                        return less_than_closest_element
                    elif greater_than_closest_dist < less_than_closest_dist:
                        return greater_than_closest_element
                    else:
                        return less_than_closest_element.union(greater_than_closest_element)
                elif range_dist < greater_than_closest_dist:
                    # Add the element to the greater than set after clearing
                    # current contents
                    greater_than_closest_element.clear()
                    greater_than_closest_element.add(element_id)
                    greater_than_closest_dist = range_dist
                elif range_dist == greater_than_closest_dist:
                    # Add the element to the current set
                    greater_than_closest_element.add(element_id)
                elif range_dist > greater_than_closest_dist:
                    # Return the closest elements if the current range distance
                    # is greater than the minimum right distance
                    if less_than_closest_dist < greater_than_closest_dist:
                        return less_than_closest_element
                    elif greater_than_closest_dist < less_than_closest_dist:
                        return greater_than_closest_element
                    else:
                        return less_than_closest_element.union(greater_than_closest_element)
        # Return the closest if finished going through all elements
        if len(overlapping) > 0:
            return overlapping
        elif less_than_closest_dist < greater_than_closest_dist:
            return less_than_closest_element
        elif greater_than_closest_dist < less_than_closest_dist:
            return greater_than_closest_element
        elif less_than_closest_dist == greater_than_closest_dist:
            return less_than_closest_element.union(greater_than_closest_element)                    
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
