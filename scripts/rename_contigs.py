#!/usr/bin/env python

import sys
import argparse
import os
from Bio import SeqIO


parser = argparse.ArgumentParser(
        description="Convert assembly to standard naming convention"
)

parser.add_argument("--assembly", dest="infile", help="Assembly file (fasta)",  required= True)
parser.add_argument("--species", dest="species", help="Species name", required=True)
parser.add_argument("--strain", dest="strain", help="strain ID", required=True)
parser.add_argument("--version", dest="version", help="Assembly version (no'v')", required=True)
parser.add_argument("-sort", dest="sort", action="store_true", help="Should contigs be sorted an named by size")
parser.add_argument("-final", dest="final", action="store_true", help="Is this a final assembly (in which case total assembly length is added to output name)")




def name_all(recs, species_name, strain, version, sort=False, keep_m=True):
    """
    THIS IS A GENERATOR FUNCTION, it doesn't return anything, rather you can
    iterate through the contents of recs, producing one updated record after
    anotehr without loading all records into memory

    recs: generator of Biopython SeqRecords from SeqIO.parse or similar
    species_name: species biomial to include in header info
    strain: strain ID to inlclude in header info
    version: assemble versoin ("v" will be preprended in description)
    sort: if True, sort contigs by size an name accordingly (default False)
    keep_m: if True, a single contig with ID ending "_m" will not be renamed to
    have a number (for keeping mtDNA seperate).
    """
    try:
        gen,sp = species_name.split(" ")
    except ValueError:
        sys.exit("species should be in form 'Genus species', can't parse", sp)
    tag = "{}{}_{}".format(gen[0], sp[:2], strain)
    if sort:
        recs = sorted(recs, key=lambda x: len(x),reverse=True)
    # what to add to idx below to get a chrom number. Note, if we get a mito
    # contig before the last one we don't want add anything to the counter
    # {0,1,2} -> {1,mt,2}
    offset = 1
    seen_mito = False
    for idx, rec in enumerate(recs):
        if rec.id.endswith("_m"): #it's the mito
            if seen_mito:
                sys.exit("Only mitochrondial loci should end '_m', but multiple contigs ends this way")
            offset = 0
            num = "m"
        else:
            num = idx+ offset
        
        rec.id = "{}_{}".format(tag, num)
        new_desc = "{} {} {} {} v{}".format(species_name, strain, num, len(rec), version)
        rec.description = new_desc
        yield(rec)



if __name__ == "__main__":
    args = parser.parse_args()
    recs = SeqIO.parse(args.infile, "fasta")
    out_recs = name_all(recs, args.species, args.strain, args.version, args.sort)
    if args.final:
        #need to know total genome size, so load all updated recs into memory
        out_recs = list(out_recs)
        total_len = sum([len(r) for r in out_recs])
    else:
        total_len = "XXX"
    try:
        gen,sp = args.species.split(" ")
    except ValueError:
        sys.exit("species should be in form 'Genus species', can't parse", sp)
    tag = "{}{}{}".format(gen[0], sp[:2], args.strain)
    out_name = "{}_{}_{}_{}_v{}.fna".format(tag, args.species.replace(" ", "_"), args.strain, total_len, args.version)
    if os.path.exists(out_name):
        sys.exit("Assembly file", out_name, "already exists")
    n = SeqIO.write(out_recs, out_name, "fasta")
    sys.stdout.write("Wrote {} records to {}\n".format(n, out_name))
    sys.exit(0)



    




