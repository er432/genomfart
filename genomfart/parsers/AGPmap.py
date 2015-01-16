import numpy as np
import re
import numba as nb
from numba import jit
from bisect import bisect_left

@jit(
    argtypes=(nb.int64, nb.int64, nb.int64, nb.int64, nb.int64[:]),
    restype=nb.int64[:]
    )
def _getInterval(chromosome, position, start, end, markerPosition):
    """ Gets the markers, position of marker in bp flanking a position
    on the chromosome

    Parameters
    ----------
    chromosome : int
       The chromosome
    position : int
       The position
    start : int
       Position to start search
    end : int
       Position to stop search
    markerPosition : np.ndarray
       Array of marker positions

    Returns
    -------
    [index of marker1, index of marker2]
    Note that if position past leftmost or rightmost marker, -1 will be returned
    for that position
    """
    if position < markerPosition[start]:
        return np.array([-1,start])
    if position > markerPosition[end]:
        return np.array([end,-1])
    if position == markerPosition[start]:
        return np.array([start,start])
    if position == markerPosition[end]:
        return np.array([end,end])
    left = start
    right = end
    while ((right - left) > 1):
        mid = left + (right-left)/2
        if (position == markerPosition[mid]):
            return np.array([mid,mid])
        if (position < markerPosition[mid]):
            right = mid
        else:
            left = mid
    return np.array([left, right])

@jit(
    argtypes=(nb.int64, nb.int64, nb.int64, nb.int64, nb.int64[:],
              nb.float64[:]),
    restype=nb.float64
    )
def _getCmFromPosition(chromosome, position, start, end, markerPosition,
                       markercm):
    """ Gets cM position from chromosome bp

    Parameters
    ----------
    chromosome : int
       The chromosome
    position : int
       The position
    start : int
       Position to start search
    end : int
       Position to stop search
    markerPosition : np.ndarray
       Array of marker positions
    markercm : np.ndarray
       Array of marker cM positions

    Returns
    -------
    cM position
    """
    if position < markerPosition[start]:
        g2pRatio = (markercm[start + 10] - markercm[start]) \
          / (markerPosition[start + 10] - markerPosition[start])
        return markercm[start] - (markerPosition[start] - position) * g2pRatio
    if position > markerPosition[end]:
        g2pRatio = (markercm[end] - markercm[end - 10]) / (markerPosition[end] - \
            markerPosition[end - 10])
        return markercm[end] + (position - markerPosition[end]) * g2pRatio
    if position == markerPosition[start]: return markercm[start]
    if position == markerPosition[end]: return markercm[end]
    left = start
    right = end
    while ((right-left) > 1):
        mid = left+(right-left)/2
        if (position == markerPosition[mid]): return markercm[mid]
        if (position < markerPosition[mid]): right = mid
        else:
            left = mid
    g2pRatio = (markercm[right] - markercm[left]) / (markerPosition[right]\
        - markerPosition[left]);
    return markercm[left] + (position - markerPosition[left]) * g2pRatio

@jit(
    argtypes=(nb.int64, nb.float64, nb.int64, nb.int64, nb.int64[:],
              nb.float64[:]),
    restype=nb.int64
    )
def _getPositionFromCm(chromosome, cM, start, end, markerPosition,
                       markercm):
    """ Gets chromosome bp from cM position

    Parameters
    ----------
    chromosome : int
       The chromosome
    cM : float
       The position in cM
    start : int
       Position to start search
    end : int
       Position to stop search
    markerPosition : np.ndarray
       Array of marker positions
    markercm : np.ndarray
       Array of marker cM positions

    Returns
    -------
    bp position
    """
    if (cM < markercm[start]):
        p2gRatio = (markerPosition[start + 10] - \
            markerPosition[start]) / (markercm[start + 10] - markercm[start])
        return int(markerPosition[start] - np.round((markercm[start] - cM) * p2gRatio))
    if (cM > markercm[end]):
        p2gRatio = (markerPosition[end] - markerPosition[end - 10]) / \
          (markercm[end] - markercm[end - 10])
        return int(markerPosition[end] + np.round((cM - markercm[end]) * p2gRatio))
    if (cM == markercm[start]): return markerPosition[start]
    if (cM == markercm[end]): return markerPosition[end]
    left = start
    right = end
    while ((right-left)>1):
        mid = left + (right-left)/2
        if (cM == markercm[mid]): return markerPosition[mid]
        if (cM < markercm[mid]): right = mid
        else:
            left = mid
    p2gRatio = (markerPosition[right] - markerPosition[left]) / (markercm[right] \
            - markercm[left])
    return markerPosition[left] + int(np.round((cM - markercm[left]) * p2gRatio))

class AGPMap:
    """ Class used to parse and manipulate data from an AGPmap with cM positions.

    The map can be in 1 of 2 different formats, each of which is tab-delimited

    The first has the chromosome first, the marker third, the cM position 4th,
    the AGP position 5th

    The second has the marker first, the chromosome second, the cM position third,
    and the AGP position 6th
    """
    def __init__(self, mapFile, useAgpV2=False):
        """ Instantiates the AGPMap parser

        Parameters
        ----------
        mapFile : str
            The filename for the map
        useAgpV2 : boolean
            Whether this is version 2 of the AGPmap format
        """
        colDict = {'chrom':1 if useAgpV2 else 0,
                   'marker':0 if useAgpV2 else 2,
                   'markercm':2 if useAgpV2 else 3,
                   'markerpos':5 if useAgpV2 else 4}
        # Keep track of where chromosomes end
        self.chrend = {}
        # Keep track of marker physical
        self.markerPosition = []
        # Keep track of marker chromosomes
        self.markerChromosome = []
        # Keep track of marker cM positions
        self.markercm = []
        # Keep track of marker names
        self.marker = []
        # Keep track of alternative marker names if necessary
        self.marker_alt_names = []
        # Load data
        mapFile = open(mapFile)
        header = mapFile.readline()
        for i,line in enumerate(mapFile):
            line = line.strip().split('\t')
            chrom = int(line[colDict['chrom']])
            self.chrend[chrom-1] = i
            self.marker.append(line[colDict['marker']])
            self.markercm.append(float(line[colDict['markercm']]))
            self.markerPosition.append(int(line[colDict['markerpos']]))
            if not useAgpV2:
                self.marker_alt_names.append(line[1])
            self.markerChromosome.append(chrom)
        # Convert arrays to numpy arrays
        self.markerPosition = np.array(self.markerPosition)
        self.markerChromosome = np.array(self.markerChromosome)
        self.markercm = np.array(self.markercm)
        self.marker = np.array(self.marker)
    def getInterval(self, chromosome, position):
        """ Gets the markers, position of marker in bp flanking a position
        on the chromosome

        Parameters
        ----------
        chromosome : int
           The chromosome
        position : int
           The position

        Returns
        -------
        (left marker, left marker position, right marker, right marker position)
        """
        end = self.chrend[chromosome-1]
        start = 0
        if chromosome > 1:
            start = self.chrend[chromosome-2]+1
        left_ind,right_ind = _getInterval(chromosome, position, start, end,
                                          self.markerPosition)
        return (self.marker[left_ind] if left_ind != -1 else None,
                self.markerPosition[left_ind] if left_ind != -1 else None,
                self.marker[right_ind] if right_ind != -1 else None,
                self.markerPosition[right_ind] if right_ind != -1 else None)
    def getFlankingMarkerIndices(self, chromosome, geneticPosition):
        """ Gets the indices of the markers flanking a given genetic position

        Parameters
        ----------
        chromosome : int
           The chromosome
        geneticPosition : float
           genetic position in cM

        Returns
        -------
        left flank index, right flank index
        """
        frm = 0
        if chromosome > 1:
            frm = self.chrend[chromosome-2]+1
        to = self.chrend[chromosome-1]+1
        ndx = bisect_left(self.markercm, geneticPosition, lo=frm, hi=to)
        if self.markercm[ndx] == geneticPosition:
            return ndx,ndx
        else:
            return (ndx-1,ndx)
    def getCmFromPosition(self, chromosome, position):
        """ Gets cM position from chromosome bp

        Parameters
        ----------
        chromosome : int
           The chromosome
        position : int
           The position

        Returns
        -------
        cM position
        """
        end = self.chrend[chromosome-1]
        start = 0
        if chromosome > 1:
            start = self.chrend[chromosome-2]+1
        return _getCmFromPosition(chromosome, position,
                                  start, end,
                                  self.markerPosition, self.markercm)
    def getPositionFromCm(self, chromosome, cM):
        """ Gets chromosome bp from cM position

        Parameters
        ----------
        chromosome : int
           The chromosome
        cM : float
           The position in cM

        Returns
        -------
        The chromosome bp position
        """
        end = self.chrend[chromosome-1]
        start = 0
        if chromosome > 1:
            start = self.chrend[chromosome-2]+1
        return _getPositionFromCm(chromosome, cM,
                                  start, end,
                                  self.markerPosition, self.markercm)
    def getFirstGeneticPosition(self, chromosome):
        """ Gets the first genetic position on a chromosome

        Parameters
        ----------
        chromosome : int
           The chromosome

        Returns
        -------
        The first genetic position (in cM)
        """
        if chromosome == 1: return self.markercm[0]
        else:
            return self.markercm[self.chrend[chromosome-2]+1]
    def getFirstMarkerName(self, chromosome):
        """ Gets the name of the first marker

        Parameters
        ----------
        chromosome : int
            The chromosome to get the first marker from
        
        Returns
        -------
        The name of the first marker
        """
        if chromosome == 1: return self.marker[0]
        else:
            return self.marker[self.chrend[chromosome-2]+1]
    def getLastMarkerName(self, chromosome):
        """ Gets the name of the last marker

        Parameters
        ----------
        chromosome : int
            The chromosome to get the last marker from

        Returns
        -------
        The name of the last marker
        """
        return self.marker[self.chrend[chromosome-1]]
    def getLastGeneticPosition(self, chromosome):
        """ Gets the last genetic position on a chromosome

        Parameters
        ----------
        chromosome : int
           The chromosome

        Returns
        -------
        The last genetic position (in cM)
        """
        return self.markercm[self.chrend[chromosome-1]]
    def getPhysPos(self, markerIndex):
        """ Gets the physical position of a marker

        Parameters
        ----------
        markerIndex : int
           Index of the marker

        Returns
        -------
        Physical position of marker
        """
        return self.markerPosition[markerIndex]
    def getGeneticPos(self, markerIndex):
        """ Gets the genetic position of a marker

        Parameters
        ----------
        markerIndex : int
           Index of the marker

        Returns
        -------
        Genetic position of the marker
        """
        return self.markercm[markerIndex]
    def getMarkerNumber(self, marker_name):
        """ Gets the number of the marker

        Parameters
        ----------
        marker_name : str
            The name of the marker

        Returns
        -------
        Number of the marker
        """
        return int(re.search('\d+',marker_name).group())
