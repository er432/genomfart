from bisect import bisect_right, bisect_left
from Ranger import Range
from genomfart.utils.bigDataFrame import BigDataFrame
import networkx as nx
import sys

if sys.version_info[0] > 2:
    xrange = range

## Instantiated class used to map between two versions of an assembly
class version_mapper:
    ## Graph, where nodes are named as (v1/v2, chrom, end_pt_of_region).
    #regions are connected by edges
    graph = None
    ## Dictionary of chrom -> [end points of regions in v1 assembly]
    # Same order as v2_end_dict for matching region
    v1_end_dict = {}
    ## Dictionary of chrom -> [end points of regions in v2 assembly]
    # Same order as v1_end_dict for matching region
    v2_end_dict = {}
    ## Instantiates the version mapper
    # @param map_file A file that maps between the two assemblies. It should have
    # columns v1_chrom, v1_start, v1_end, v2_chrom, v2_start, v2_end, orientation.
    # It is also assumed to be base-1
    def __init__(self, map_file):
        """
        Instantiates the version mapper

        Parameters
        ----------
        map_file : str
            A file that maps between the two assemblies. It should have columns
            v1_chrom, v1_start, v1_end, v2_chrom, v2_start, v2_end, orientation.
            It is also assumed to be base-1
        """
        map_frame = BigDataFrame(map_file,header=True, assume_uniform_types=False)
        ## Insantiate the dictionaries
        self.graph = nx.Graph()
        self.v1_end_dict = {}
        self.v2_end_dict = {}
        ## Fill the dictionaries
        for i,row in enumerate(map_frame):
            # Add chromosomes to endpoint dictionaries if necessary
            if row['v1_chrom'] not in self.v1_end_dict:
                self.v1_end_dict[row['v1_chrom']] = []
            if row['v2_chrom'] not in self.v2_end_dict:
                self.v2_end_dict[row['v2_chrom']] = []
            # Put in start and end points
            self.v1_end_dict[row['v1_chrom']].append(row['v1_end'])
            self.v2_end_dict[row['v2_chrom']].append(row['v2_end'])
            # Put in edges in graph
            node1 = ('v1',row['v1_chrom'],row['v1_end'])
            node2 = ('v2',row['v2_chrom'],row['v2_end'])
            self.graph.add_edge(node1,node2, strand=row['orientation'])
            self.graph.node[node1]['start'] = row['v1_start']
            self.graph.node[node1]['end'] = row['v1_end']
            self.graph.node[node1]['chrom'] = row['v1_chrom']
            self.graph.node[node2]['start'] = row['v2_start']
            self.graph.node[node2]['end'] = row['v2_end']
            self.graph.node[node2]['chrom'] = row['v2_chrom']
        ## Sort chromosome lists
        for dic in (self.v1_end_dict,self.v2_end_dict):
            for k in dic:
                dic[k] = sorted(dic[k])
    ## Gets the position in the version 2 genome if it exists (base 1 assumed)
    # @param v1_chr The version 1 chromosome
    # @param v1_pos The version 1 position
    # @return (chr,pos,orientation) if position exists in version 2 or None if not
    def v1_to_v2_map(self, v1_chr, v1_pos):
        """
        Gets the position of the version 2 genome if it exists (base 1 assumed)

        Parameters
        ----------
        v1_chr : hashable
            The version 1 chromosome
        v1_pos : int
            The version 1 position

        Returns
        -------
        (v2_chrom,v2_pos,orientation relative to v1) if the position exists in version 2
        or None if it doesn't exist in version 2
        """
        # Get the index of the end position within v1 in the end dict
        v1_end_ind = bisect_left(self.v1_end_dict[v1_chr],v1_pos)
        if v1_end_ind >= len(self.v1_end_dict[v1_chr]): return None
        # Define the region node for v1
        v1_node = ('v1',v1_chr,self.v1_end_dict[v1_chr][v1_end_ind])
        v2_node = self.graph.neighbors(v1_node)[0]
        # Check if the position is actually within the v1 interval
        # If not, return None
        if not self.graph.node[v1_node]['start'] <= v1_pos <= self.graph.node[v1_node]['end']:
            return None
        else:
            orientation = self.graph.edge[v1_node][v2_node]['strand']
            spacing = v1_pos - self.graph.node[v1_node]['start']
            if orientation == 1:
                return (v2_node[1], self.graph.node[v2_node]['start'] + spacing,orientation)
            else:
                return (v2_node[1], self.graph.node[v2_node]['end'] - spacing,orientation)
    ## Gets the position in the version 1 genome if it exists (base 1 assumed)
    # @param v2_chr The version 2 chromosome
    # @param v2_pos The version 2 position
    # @return (chr,pos,orientation) if position exists in version 2 or None if not
    def v2_to_v1_map(self, v2_chr, v2_pos):
        """
        Gets the position of the version 1 genome if it exists (base 1 assumed)

        Parameters
        ----------
        v2_chr : hashable
            The version 2 chromosome
        v2_pos : int
            The version 2 position

        Returns
        -------
        (v1_chrom,v1_pos,orientation relative to v2) if the position exists in version 1
        or None if it doesn't exist in version 1
        """
        # Get the index of the end position within v2 in the end dict
        v2_end_ind = bisect_left(self.v2_end_dict[v2_chr],v2_pos)
        if v2_end_ind >= len(self.v2_end_dict[v2_chr]): return None
        # Define the region node for v2
        v2_node = ('v2',v2_chr,self.v2_end_dict[v2_chr][v2_end_ind])
        v1_node = self.graph.neighbors(v2_node)[0]
        # Check if the position is actually within the v2 interval
        # If not, return None
        if not self.graph.node[v2_node]['start'] <= v2_pos <= self.graph.node[v2_node]['end']:
            return None
        else:
            orientation = self.graph.edge[v2_node][v1_node]['strand']
            spacing = v2_pos - self.graph.node[v2_node]['start']
            if orientation == 1:
                return (v1_node[1], self.graph.node[v1_node]['start'] + spacing,orientation)
            else:
                return (v1_node[1], self.graph.node[v1_node]['end'] - spacing,orientation)
    def v1_to_v2_seg_map(self, v1_chr, v1_start, v1_end):
        """
        Gets the ranges of positions in the version 2 genome if they exist (base 1 assumed)

        Parameters
        ----------
        v1_chr : hashable
            The version 1 chromosome
        v1_start : int
            The version 1 start position (inclusive)
        v1_end : int
            The version 1 end position (inclusive)

        Returns
        -------
        Dictionary of (v1_chrom,v1_start,v1_end)->(v2_chrom,v2_start,v2_end,orientatoin relative to v1) for
        ranges of positions existing in version 2 
        """
        # Get index of version 1 segment corresponding to end position in the end dict
        v1_end_ind = min(bisect_left(self.v1_end_dict[v1_chr],v1_end),len(self.v1_end_dict[v1_chr])-1)
        # Get index of version 1 segment corresponding to start position in the end dict
        v1_start_ind = min(bisect_left(self.v1_end_dict[v1_chr],v1_start),len(self.v1_end_dict[v1_chr])-1)
        return_dict = {}
        # Go through the range of v1 nodes
        for i in xrange(min(v1_start_ind,v1_end_ind),max(v1_start_ind,v1_end_ind)+1):
            # Define the nodes
            v1_node = ('v1',v1_chr,self.v1_end_dict[v1_chr][i])
            v2_node = self.graph.neighbors(v1_node)[0]
            # Check if there is actual overlap between v1 interval and the start+end points
            ovlap = Range.closed(v1_start,v1_end).isConnected(Range.closed(self.graph.node[v1_node]['start'], \
                                                        self.graph.node[v1_node]['end']))
            
            if ovlap:
                orientation = self.graph.edge[v1_node][v2_node]['strand']
                spacing = max(v1_start,self.graph.node[v1_node]['start']) - \
                  self.graph.node[v1_node]['start']
                if orientation == 1:
                    key = (v1_chr,self.graph.node[v1_node]['start']+spacing,min(v1_end,\
                                                self.graph.node[v1_node]['end']))
                    seg_length = key[2]-key[1]
                    val = (v2_node[1],self.graph.node[v2_node]['start']+spacing,\
                           self.graph.node[v2_node]['start']+spacing+seg_length,orientation)
                    return_dict[key] = val
                else:
                    key = (v1_chr,self.graph.node[v1_node]['start']+spacing,min(v1_end,\
                                                self.graph.node[v1_node]['end']))
                    seg_length = key[2]-key[1]
                    val = (v2_node[1],self.graph.node[v2_node]['end']-spacing, \
                           self.graph.node[v2_node]['end']-spacing-seg_length,orientation)
                    return_dict[key] = val
        return return_dict
    def v2_to_v1_seg_map(self, v2_chr, v2_start, v2_end):
        """
        Gets the ranges of positions in the version 1 genome if they exist (base 1 assumed)

        Parameters
        ----------
        v2_chr : hashable
            The version 2 chromosome
        v2_start : int
            The version 2 start position (inclusive)
        v2_end : int
            The version 2 end position (inclusive)

        Returns
        -------
        Dictionary of (v2_chrom,v2_start,v2_end)->(v1_chrom,v1_start,v1_end,orientation relative to v2) for
        ranges of positions existing in version 1 
        """
        # Get index of version 2 segment corresponding to end position in the end dict
        v2_end_ind = min(bisect_left(self.v2_end_dict[v2_chr],v2_end),len(self.v2_end_dict[v2_chr])-1)
        # Get index of version 2 segment corresponding to start position in the end dict
        v2_start_ind = min(bisect_left(self.v2_end_dict[v2_chr],v2_start),len(self.v2_end_dict[v2_chr])-1)
        return_dict = {}
        # Go through the range of v2 nodes
        for i in xrange(min(v2_start_ind,v2_end_ind),max(v2_start_ind,v2_end_ind)+1):
            # Define the nodes
            v2_node = ('v2',v2_chr,self.v2_end_dict[v2_chr][i])
            v1_node = self.graph.neighbors(v2_node)[0]
            # Check if there is actual overlap between v2 interval and the start+end points
            ovlap = Range.closed(v2_start,v2_end).isConnected(Range.closed(self.graph.node[v2_node]['start'], \
                                                        self.graph.node[v2_node]['end']))
            if ovlap:
                orientation = self.graph.edge[v2_node][v1_node]['strand']
                spacing = max(v2_start,self.graph.node[v2_node]['start']) - \
                  self.graph.node[v2_node]['start']
                if orientation == 1:
                    key = (v2_chr,self.graph.node[v2_node]['start']+spacing,min(v2_end,\
                                                self.graph.node[v2_node]['end']))
                    seg_length = key[2]-key[1]
                    val = (v1_node[1],self.graph.node[v1_node]['start']+spacing,\
                           self.graph.node[v1_node]['start']+spacing+seg_length,orientation)
                    return_dict[key] = val
                else:
                    key = (v2_chr,self.graph.node[v2_node]['start']+spacing,min(v2_end,\
                                                self.graph.node[v2_node]['end']))
                    seg_length = key[2]-key[1]
                    val = (v1_node[1],self.graph.node[v1_node]['end']-spacing, \
                           self.graph.node[v1_node]['end']-spacing-seg_length,orientation)
                    return_dict[key] = val
        return return_dict
        
