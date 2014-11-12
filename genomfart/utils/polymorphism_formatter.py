from os import path
from collections import deque

class plink:
    """ Static class used to work with PLINK-formatted data
    """
    @staticmethod
    def write_geno_mat_to_ped(write_base, geno_mat, samples, positions,
                              chrom, families=None):
        """ Writes MAP and PED files for a given set of genotypes

        Parameters
        ----------
        write_base : str
            The base name of the resulting plink files
        geno_mat : np.ndarray
            A numpy matrix of integer genotypes (0,1,2). All markers are assumed
            to be biallelic. Rows are samples,  and columns are positions
        samples : iterable
            The sample names, in the same order as matrix rows
        positions : iterable
            The positions of the markers, in the same order as matrix columns
        chrom : str
            The chromosome name
        families : iterable, optional
            Family name for each individual, if they are grouped by family
        """
        # First write the MAP file
        map_file = open(write_base+'.map','w')
        for i,position in enumerate(positions):
            line = [chrom,"marker%d" % i,"0",str(position)]
            map_file.write(' '.join(line)+'\n')
        map_file.close()
        # Write the ped file
        ped_file = open(write_base+'.ped','w')
        for i,samp in enumerate(samples):
            if families:
                write_line = deque([families[i],samp,'0','0','0','0'])
            else:
                write_line = deque([str(i+1),samp,'0','0','0','0'])
            for j,pos in enumerate(positions):
                if geno_mat[i,j] == 0:
                    write_line.append('1')
                    write_line.append('1')
                elif geno_mat[i,j] == 1:
                    write_line.append('1')
                    write_line.append('2')
                elif geno_mat[i,j] == 2:
                    write_line.append('2')
                    write_line.append('2')
            ped_file.write(' '.join(write_line)+'\n')
        ped_file.close()
        
