#!/usr/bin/env python
from genomfart.utils.bigDataFrame import BigDataFrame
import numpy as np
import argparse
import os

def shuffle_pheno_within_family(in_file, out_file):
    """
    Shuffle a phenotype within each family, and write to file

    Parameters
    ----------
    in_file : str
        Name of the PLINK-formatted phenotype file
    out_file : str
        Name to give the output file
    """
    # Open up file to shuffle
    frame = BigDataFrame(in_file, header=False, assume_uniform_types=False)
    ## Create dictionary of family -> values
    pheno_vals = {}
    for row in frame:
        if row[0] not in pheno_vals:
            pheno_vals[row[0]] = []
        pheno_vals[row[0]].append(row[2])
    ## Make permuted file
    write_file = open(out_file,'w')
    count = 0
    for fam, vals in pheno_vals.iteritems():
        np.random.shuffle(vals)
        for val in vals:
            write_line = [str(fam),str(count),"%0.8g" % val if type(val) != str else val]
            write_file.write('\t'.join(write_line)+'\n')
            count += 1
    write_file.close()

def main():
    parser = argparse.ArgumentParser(description="Shuffle a phenotype within family in a PLINK-phenotype file")
    parser.add_argument('in_file',metavar='in_file',nargs='?',help='Input file (PLINK format)')
    parser.add_argument('out_file',metavar='out_file',nargs='?',help='Output file (PLINK format)')
    args = parser.parse_args()
    shuffle_pheno_within_family(args.in_file,args.out_file)

if __name__ == "__main__":
    main()
