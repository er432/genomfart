from genomfart.parsers.AGPmap import AGPMap
from genomfart.utils.bigDataFrame import BigDataFrame
from numba import jit
import numba as nb
import numpy as np
import re

@jit(argtypes=(
    nb.int64[:], nb.int64[:], nb.double[:], nb.double[:,:], nb.int64,
    nb.int64, nb.int64, nb.double), restype=nb.void)
def _projectSnpBoolean(parents, popIndex, snpvalues, genotypes, nSamples,
                       leftmarker, rightmarker, pd):
    """ Projects SNPs if parents have boolean genotypes
    """
    for i in xrange(nSamples):
        if parents[popIndex[i]] == 0:
            snpvalues[i] = 0
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
        self.firstMarker = None
        self.sampleNameMap = {}
        self.genotypes = []
        self.samp_names = []        
        self.importMarkersForMap(rilFile)
        self.maxMarker = self.genotypes.shape[1]-1
    def importMarkersForMap(self, rilFile):
        """ Reads the data file for a chromosome with the sample allele states

        Parameters
        ----------
        rilFile : str
            The filename for the RIL file
        """
        rilFile = BigDataFrame(rilFile, header=True, assume_uniform_types=False)
        nGenos = len(rilFile.colLabels)-5
        self.genotypes = [[] for x in xrange(nGenos)]
        self.sampleNameMap = dict((k,v-5) for k,v in rilFile.colLabels.items() if v >= 5)
        self.samp_names = sorted(self.sampleNameMap.keys(),
                                 key=lambda x: self.sampleNameMap[x])
        count = 0
        nMarkers = 0
        self.firstMarker = self.theAGPMap.getMarkerNumber(self.theAGPMap.getFirstMarkerName(self.chromosome))
        for row in rilFile.iterrows(cache=False):
            for samp in self.samp_names:
                self.genotypes[self.sampleNameMap[samp]].append(row[samp])
        self.genotypes = np.array(self.genotypes)
        # header = rilFile.readline().strip().split('\t')
        # self.firstMarker = int(re.search('\d+',header[2]).group())
        # nMarkers = len(header)-1
        # count = 0
        # while 1:
        #     data = rilFile.readline().strip()
        #     if data == '': break
        #     data = data.split('\t')
        #     self.sampleNameMap[data[0]] = count
        #     self.genotypes.append([])
        #     self.samp_names.append(data[0])
        #     for j in xrange(nMarkers):
        #         self.genotypes[-1].append(float(data[j+1]))
        # self.genotypes = np.array(self.genotypes)
        # rilFile.close()
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
        pd = 0.
        if right != left:
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
    def projectAllSnps(self, chrom_length, popIndex, founder_data, boolean = True,
                       positions = None, enumerate_snps=False):
        """ Projects all SNPs on a chromsome onto descendants

        Parameters
        ----------
        chrom_length : int
            The length of the chromosome
        popIndex : np.ndarray, int
            Indices of the population for each sample
        founder_data : genomfart.parsers.SNPdata
            Object containing the founder SNP data
        boolean : boolean
            Whether the parents returned have boolean values for SNPs
        positions : set of ints
            A set of specific positions to project. If None, all SNPs will
            be projected
        enumerate_snps : boolean
            Whether to return both the absolute value of the SNP index along with the positions

        Returns
        -------
        Generator non-reference allele count for each RIL, as an array of
        doubles

        Examples
        --------
        >>> from genomfart.utils.snp_projector import snp_projector
        >>> from genomfart.parsers.SNPdata import SNPdata
        >>> import numpy as np
        >>> import re
        >>> import os
        >>> chrom = 10
        >>> chrom_length = 149686046
        >>> home_dir = os.path.expanduser('~')
        >>> map_dir = os.path.join(home_dir,'Documents','Zea','Data','Sequence',
        ... 'NAM_imputed_maps','imputedMarkers.0.1cm')
        >>> founder_dir = os.path.join(home_dir,'Documents','Zea','Data','Sequence',
        ... 'NAM_founder_genos')
        >>> founder_file = os.path.join(founder_dir,
        ... 'gwas_snps_union_chr%d.h1.h2_minors_indels.cnv_genic_2kwin_500bpbin.20130605.txt' % chrom)
        >>> map_file = os.path.join(map_dir,'Imputed_0.1cm_master.txt')
        >>> ril_file = os.path.join(map_dir,'Imputed_0.1cm_chr%d_clean.txt' % chrom)
        >>> projector = snp_projector(chrom,chrom_length,map_file,ril_file)
        >>> popIndex = np.array(map(lambda x: (17 if x[0] == 'M' else int(re.search('(?<=Z)\d+(?=E)',x).group()))-1,)
        ... projector.samp_names), dtype=np.int64)
        >>> founder_data = SNPdata(chrom, founder_file, 'B73')
        >>> projection = projector.projectAllSnps(chrom_length, popIndex, founder_data)
        """
        if not boolean:
            raise NotImplementedError("Non-boolean projection not yet implemented")
        founder_data.reset()
        ind = 0
        while (founder_data.next()):
            if positions is not None:
                if founder_data.getPosition() not in positions:
                    ind += 1
                    continue
            if boolean:
                parents = np.array(founder_data.getGenotype(), dtype=np.int64)
                if enumerate_snps:
                    yield ind,self.projectSnpBoolean(parents, founder_data.getPosition(),
                        chrom_length, popIndex)
                else:
                    yield self.projectSnpBoolean(parents, founder_data.getPosition(),
                        chrom_length, popIndex)
            else:
                pass
            ind += 1

        
