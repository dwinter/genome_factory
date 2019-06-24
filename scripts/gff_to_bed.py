#! /usr/bin/env python
import sys

def gff_to_bed(line):
    """ Take a gff line, return a bed line """

    toks = line.split("\t")
    start, end, score, strand = toks[3:7] 
    tmpl = '{}\t{}\t{}\t.\t{}\t{}\n'
    return(tmpl.format(toks[0], int(start)-1, end, score, strand))

if __name__ == "__main__":
    try:
        in_file = sys.argv[1]
    except ValueError:
        print("Usage gff_to_bed.py [in.gff] > out.bed")
    with open(in_file) as input:
        n = 0
        for L in input:
            if L.startswith("#"):
                pass
            else:
                sys.stdout.write(gff_to_bed(L))
                n += 1
    sys.stderr.write("Wrote {} intervals".format(n))
