from genomfart.parsers.AGPmap import AGPMap
from numba import jit
import numba as nb
import numpy as np

@jit(argtypes=(
    nb.int64[:], nb.int64[:], np.float64[:], nb.float64[:,:], nb.int64,
    nb.int64, nb.int64, np.float64), restype=nb.void)
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
        self.firstMarker = self.getMarkerNumber(self.theAGPMap.getFirstMarkerName(chromosome))
        self.maxMarker = self.getMarkerNumber(self.theAGPMap.getLastMarkerName(chromosome))
        self.sampleNameMap = {}
        self.genotypes = []
    def importMarkersForMap(self, rilFile):
        """ Reads the data file for a chromosome with the sample allele states

        Parameters
        ----------
        rilFile : str
            The filename for the RIL file
        """
        self.genotypes = []        
        rilFile = open(rilFile)
        header = rilFile.readline.strip().split('\t')
        nMarkers = len(header)-1
        count = 0
        while 1:
            data = rilFile.readline().strip()
            if data == '': break
            data = data.split('\t')
            self.sampleNameMap[data[0]] = count
            self.genotypes.append([])
            for j in xrange(nMarkers):
                self.genotypes[-1].append(float(data[j+1]))
        self.genotypes = np.array(self.genotypes)
        rilFile.close()
    def projectSnpBoolean(parents, pos, chrom_length, popIndex):
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
        if left == -1:
            left = 0
        right = right_mark_pos
        if right == -1:
            right = self.chrom_length
        # Proportion of distance of SNP between left and right markers
        pd = (float(pos-left))/(float(right-left))
        leftmarker = self.theAGPMap.getMarkerNumber(left_mark)
        if (leftmarker == -1): leftmarker=0
        else:
            leftmarker = leftmarker-firstMarker+1
        rightmarker = self.theAGPMap.getMarkerNumber(right_mark)
        if (rightmarker == -1): rightmarker = self.maxMarker
        else:
            rightmarker = rightmarker - self.firstMarker + 1
        nSamples = len(popIndex)
        snpvalues = np.zeros(nSamples)
        _projectSnpBoolean(parents, self.popIndex, snpvalues, self.genotypes,
                           nSamples, leftmarker, rightmarker, pd)
        return snpvalues
