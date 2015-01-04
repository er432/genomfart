from genomfart.parsers.AGPmap import AGPMap
from numba import jit
import numba as nb
import numpy as np

@jit(argtypes=(
    nb.int64[:], nb.int64[:], nb.double[:], nb.double[:,:], nb.int64,
    nb.int64, nb.int64, nb.double), restype=nb.void)
def _projectSnpBoolean(parents, popIndex, snpvalues, genotypes, nSamples,
                       leftmarker, rightmarker, pd):
    """ Projects SNPs if parents have boolean genotypes
    """
    for i in xrange(nSamples):
        if parents[popIndex[i]]:
            snpvalues[i] == 0
        else:
            leftval = genotypes[i,leftmarker]
            rightval = genotypes[i,rightmarker]
            if (leftval == rightval):
                snpvalues[i] = leftval
            else:
                snpvalues[i] = leftval*(1-pd) + rightval*pd

class snp_projector:
    """ Class used to project SNPs from a set of founders onto
    descendants
    """
    def __init__(self, chromosome, chrom_length, mapFile, rilFile, useAgpV2=False):
        """ Instantiates the projector

        Parameters
        ----------
        chromosome : int
            Chromosome to project onto
        chrom_length : int
            Length of the chromosome
        mapFile : str
            The filename for the map
        rilFile : str
            The filename for the RIL file
        useAgpV2 : boolean
            Whether this is version 2 of the AGPmap format        
        """
        self.chromosome = chromosome
        self.chrom_length = chrom_length
        self.theAGPMap = AGPMap(mapFile, useAgpV2)
        self.firstMarker = self.theAGPMap.getMarkerNumber(self.theAGPMap.getFirstMarkerName(chromosome))
        self.maxMarker = self.theAGPMap.getMarkerNumber(self.theAGPMap.getLastMarkerName(chromosome))
        self.sampleNameMap = {}
        self.genotypes = []
        self.samp_names = []        
        self.importMarkersForMap(rilFile)
    def importMarkersForMap(self, rilFile):
        """ Reads the data file for a chromosome with the sample allele states

        Parameters
        ----------
        rilFile : str
            The filename for the RIL file
        """
        self.genotypes = []        
        rilFile = open(rilFile)
        header = rilFile.readline().strip().split('\t')
        nMarkers = len(header)-1
        count = 0
        while 1:
            data = rilFile.readline().strip()
            if data == '': break
            data = data.split('\t')
            self.sampleNameMap[data[0]] = count
            self.genotypes.append([])
            self.samp_names.append(data[0])
            for j in xrange(nMarkers):
                self.genotypes[-1].append(float(data[j+1]))
        self.genotypes = np.array(self.genotypes)
        rilFile.close()
    def projectSnpBoolean(self, parents, pos, chrom_length, popIndex):
        """ Projects a SNP onto descendants if parent values are boolean

        Parameters
        ----------
        parents : np.ndarray, boolean
            Array of parent genotypes
        pos : int
            Position of the SNP
        popIndex : np.ndarray, int
            Indices of the population for each sample
        """
        left_mark, left_mark_pos, right_mark, right_mark_pos = self.theAGPMap.getInterval(self.chromosome,
                                                                                          pos)
        left = left_mark_pos
        if left is None:
            left = 0
        right = right_mark_pos
        if right is None:
            right = self.chrom_length
        # Proportion of distance of SNP between left and right markers
        pd = (float(pos-left))/(float(right-left))
        leftmarker = 0
        if left_mark:
            leftmarker = self.theAGPMap.getMarkerNumber(left_mark)
            leftmarker = leftmarker-self.firstMarker+1
        rightmarker = self.maxMarker
        if right_mark:
            rightmarker = self.theAGPMap.getMarkerNumber(right_mark)
            rightmarker = rightmarker - self.firstMarker + 1
        nSamples = len(popIndex)
        snpvalues = np.zeros(nSamples)
        _projectSnpBoolean(parents, popIndex, snpvalues, self.genotypes,
                           nSamples, leftmarker, rightmarker, pd)
        return snpvalues
