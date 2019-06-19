#!/usr/bin/env python

import sys

def bed_to_region(bed_line):
    chrom, start, end = bed_line.strip().split()
    try:
        #regions are 1-based, so push the start in one
        new_start = int(start)+1
    except ValueError:
        sys.exit("Second coloumn of bed doesn't look like an number\n:\t", 
                 bed_line)        
                
    return("{}:{}-{}\n".format(chrom, int(start)+1, end))


if __name__ == "__main__":
    try:
        infile = sys.argv[1]
    except IndexError:
        sys.exit("Usage: bam_to_region.py [in.bed] > out.regions")
    for L in open(infile):
        sys.stdout.write(bed_to_region(L))
    sys.exit(0)


