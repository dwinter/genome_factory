#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.SeqUtils import GC

def gather_data(rec, n):
    return(
        [ r.id, 
          len(r), 
          round(GC(r.seq),1),
          str(r.seq[:n]), 
          str(r.seq[-n:])
        ]
    )


if __name__ == "__main__":
    try:
        infile = sys.argv[1]
    except ValueError:
        print("Usage: assembly_quick_look.py [ref.fasta]")
    recs = SeqIO.parse(infile, "fasta")
    sys.stdout.write("contig\tlen\tGC\tstart\tend\n")
    for r in recs:
        contig_data = gather_data(r,20)
        sys.stdout.write("{}\t{}\t{}\t{}\t{}\n".format(*contig_data))
    sys.exit(0)




