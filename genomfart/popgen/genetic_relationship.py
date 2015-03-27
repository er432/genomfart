from genomfart.utils.snp_projector import snp_projector
from genomfart.parsers.SNPdata import SNPdata
import numpy as np
import numba as nb
from numba import jit

@jit(argtypes=(
    nb.double[:], nb.double[:,:], nb.double, nb.int64
), restype=nb.boolean)
def _append_genetic_relationship(genotypes, relationship_mat, min_MAF, taxa_count):
    """ Calculates portion of genetic relationship accounted for by
    a single locus. Adds onto current matrix of genetic relationships

    Parameters
    ----------
    genotypes : np.ndarray, float
        Genotypes for the samples (between 0 and 2)
    relationship_mat : np.ndarray, float, square
        Matrix of current genetic relationships
    min_MAF : float
        Minimum minor allele frequence allowed for a marker
    taxa_count : int
        Number of taxa in matrix

    Returns
    -------
    Whether SNP was included
    """
    # Calculate p for this locus
    p = 0.
    for i in xrange(taxa_count):
        p += genotypes[i]
    p /= (2*taxa_count)
    if (p < min_MAF or (1.-p) < min_MAF):
        return False
    E_geno = 2*p
    V_geno = 2*p*(1-p)
    # Get genetic relationship values
    for j in xrange(taxa_count):
        for k in xrange(j+1):
            relationship_mat[j,k] += ((genotypes[j]-E_geno)*(genotypes[k]-E_geno))/V_geno
    return True

@jit(argtypes=(
    nb.double[:], nb.double[:,:], nb.int64[:,:], nb.double, nb.int64, nb.int64
), restype=nb.boolean)
def _append_genetic_relationship_with_missing(genotypes, relationship_mat, snp_counts,
                                              min_MAF, max_missing, taxa_count):
    """ Calculates portion of genetic relationship accounted for by a single locus.
    Adds onto current matrix of genetic relationships, but can account for missing data,
    by not adding to a pair if at least one of the genotypes for that pair is missing. SNP
    counts used for each pair are updated at the same time

    Parameters
    ----------
    genotypes : np.ndarray, float
        Genotypes for the samples (between 0 and 2, or -1 if missing)
    relationship_mat : np.ndarray, float, square
        Matrix of current genetic relationships
    snp_counts : np.ndarray, int, square
        Matrix of SNP count used for each pair. Should be same shape as
        relationship_mat
    min_MAF : float
        Minimum minor allele frequence allowed for a marker
    max_missing : int
        The maximum number of taxa allowed to be missing for a marker
    taxa_count : int
        Number of taxa in matrix

    Returns
    -------
    Whether SNP was included
    """
    # Calculate p for this locus
    p = 0.
    missing_count = 0
    for i in xrange(taxa_count):
        if genotypes[i] >= 0:
            p += genotypes[i]
        else:
            missing_count += 1
    p /= (2*(taxa_count-missing_count))
    if (p < min_MAF or (1.-p) < min_MAF):
        return False
    elif missing_count > max_missing:
        return False
    E_geno = 2*p
    V_geno = 2*p*(1-p)
    # Get genetic relationship values
    for j in xrange(taxa_count):
        if genotypes[j] < 0: continue
        for k in xrange(j+1):
            if genotypes[k] < 0: continue
            relationship_mat[j,k] += ((genotypes[j]-E_geno)*(genotypes[k]-E_geno))/V_geno
            snp_counts[j,k] += 1
    return True

class genetic_relationship:
    """ Class used to calculate pairwise genetic relationships between taxa
    """
    @staticmethod
    def get_genetic_relationships(snp_generator, min_MAF=0.025, verbosity = None):
        """ Gets the genetic relationship matrix, as defined in the GCTA paper
        (e.g. eq. 3 from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/)

        Note here that all pairs are assumed to have the same number of SNPs
        used for the calculation of their relationship (i.e. not accounting for
        missing genotypes)

        Parameters
        ----------
        snp_generator : generator
            Generator of numpy arrays of double, giving the genotype of each
            sample as a number between 0 and 2.
        min_MAF : float
            Minimum minor allele frequency for a marker
        verbosity : int, optional
            If not None, how often to print to screen how many snps have
            been processed
        Returns
        -------
        (Numpy square array, in which the lower diagonal contains the
        genetic relationship values. Order of rows and columns is the same
        as in the generator), number of markers used
        """
        relationship_array = None
        taxa_count = 0
        marker_count = 0
        iter_count = 0
        while 1:
            try:
                genotypes = next(snp_generator)
                if relationship_array is None:
                    taxa_count = len(genotypes)
                    relationship_array = np.zeros((taxa_count, taxa_count))
                included = _append_genetic_relationship(genotypes, relationship_array,
                            min_MAF, taxa_count)
                if included:
                    marker_count += 1
                iter_count += 1
                if verbosity:
                    if iter_count % verbosity == 0:
                        print "%d markers processed..." % iter_count
                        print (np.max(relationship_array), np.min(relationship_array),
                               np.any(np.isnan(relationship_array)))
            except StopIteration:
                break
            except ValueError:
                continue
        return (1./marker_count)*relationship_array, marker_count
    @staticmethod
    def get_genetic_relationships_with_missing(snp_generator, min_MAF=0.025,
            max_missing = 500, verbosity = None):
        """ Gets the genetic relationship matrix, as defined in the GCTA paper
        (e.g. eq. 3 from http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3014363/)

        Here, pairs can have different numbers of non-missing SNPs. Essentially,
        missing data gets imputed to the mean genotype, which zeros that entry out
        in the relationship equation.

        Parameters
        ----------
        snp_generator : generator
            Generator of numpy arrays of double, giving the genotype of each
            sample as a number between 0 and 2.
        min_MAF : float
            Minimum minor allele frequency for a marker
        max_missing : int
            The maximum number of taxa allowed to be missing for a marker            
        verbosity : int, optional
            If not None, how often to print to screen how many snps have
            been processed
        Returns
        -------
        (Numpy square array, in which the lower diagonal contains the
        genetic relationship values. Order of rows and columns is the same
        as in the generator), (Numpy square array of same dimension as the
        relationship matrix, where lower diagonal contains number of SNPs
        used to calculate that relationship)
        """
        relationship_array = None
        taxa_count = 0
        marker_count = None
        iter_count = 0
        while 1:
            try:
                genotypes = next(snp_generator)
                if relationship_array is None:
                    taxa_count = len(genotypes)
                    relationship_array = np.zeros((taxa_count, taxa_count))
                    marker_count = np.zeros((taxa_count, taxa_count),
                                            dtype=np.int64)
                included = _append_genetic_relationship_with_missing(genotypes, relationship_array,
                            marker_count, min_MAF, max_missing, taxa_count)
                iter_count += 1
                if verbosity:
                    if iter_count % verbosity == 0:
                        print "%d markers processed..." % iter_count
                        print (np.max(relationship_array), np.min(relationship_array),
                               np.any(np.isnan(relationship_array)))
            except StopIteration:
                break
            except ValueError:
                continue
        return np.multiply((1./marker_count),relationship_array), marker_count
        
