import numpy as np
import sys

if sys.version_info.major > 2:
    xrange = range

class SNPdata:
    def __init__(self, chromosome, snp_file, ref_samp, samp_start_col=11,
                 totalSNPnumber=None):
        """ Instantiates a reader of SNP data on founders

        Parameters
        ----------
        chromosome : int
            The chromosome
        snp_file : str
            Name of file containing SNP data for the founders
        ref_samp : str
            The name of the reference sample (one of the columns, the one containing
            the allele state for 0)
        """
        self.chromosome = chromosome
        self.snp_file = snp_file
        self.parsedLine = None
        self.samp_start_col = samp_start_col
        self.numberOfSnps = 0
        # Open up a handle on the SNP file
        self.br = open(snp_file)
        header = self.br.readline().strip().split('\t')
        # Get the names of the samples
        self.samps = header[samp_start_col:]
        self.ref_samp = ref_samp
        self.ref_ind = header.index(self.ref_samp)
        self.col_length = len(header)
        if not totalSNPnumber:
            self.findTotalSnpNumber()
        self.reset()
    def next(self):
        """ Reads the next line of the SNPs file
        """
        if self.br is None: return False
        inLine = self.br.readline().strip()
        if inLine == '':
            self.br.close()
            return False
        else:
            self.parsedLine = inLine.strip().split('\t')
        return True
    def getGenotype(self, skip_ref = True):
        """ Gets the genotype of the founders at the current SNP

        Parameters
        ----------
        skip_ref : boolean
            Whether to skip the reference sample
        
        Returns
        -------
        Array of genotypes for each of the founders
        """
        geno = np.zeros(len(self.samps) -(1 if skip_ref else 0))
        subtract_num = 0
        for count,i in enumerate(xrange(self.samp_start_col,self.col_length)):
            if skip_ref and i == self.ref_ind:
                subtract_num = 1
                continue
            geno[count-subtract_num] = float(self.parsedLine[i])
        return geno
    def getPosition(self):
        """ Gets the current position on the chromosome

        Returns
        -------
        Position on the chromosome
        """
        return int(self.parsedLine[3])
    def getNumberOfSnps(self):
        """ Gets the number of SNPs in the file

        Returns
        -------
        The number of SNPs in the file
        """
        return self.numberOfSnps
    def findTotalSnpNumber(self):
        """ Finds the total number of SNPs in the file
        """
        self.reset()
        for i,line in enumerate(self.br):
            pass
        self.numberOfSnps = i+1
    def getAllele(self):
        """ Gets the allele configuration

        Returns
        -------
        The allele configuration
        """
        return self.parsedLine[1]
    def reset(self):
        """ Resets the reader to the first data line
        """
        self.br.close()
        self.br = open(self.snp_file)
        header = self.br.readline()
