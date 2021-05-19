import sys
import argparse
import pyfastx
import csv
from uuid import uuid4
import os
import subprocess as sp
from tqdm import tqdm


def revcom(s):
        """Return reverse complement of a sequence"""
        def complement(s):
                        basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
                        letters = list(s)
                        letters = [basecomplement[base] for base in letters]
                        return ''.join(letters)
        return complement(s[::-1])

def main(args):
    blacklist = set([s.strip() for s in open(args.blacklist).readlines()]) if args.blacklist else set()
    fasta = pyfastx.Fasta(args.msa)
    fasta_dict = {}
    for seq in fasta:
        fasta_dict[seq.name] = str(fasta[seq.name]).upper()

    seqnames = set([seq.name for seq in fasta])
    inclade = set()
    for row in csv.DictReader(open(args.clade_csv)):
        if row["Clade"]==args.clade and row["Taxa"] in fasta:
            inclade.add(row["Taxa"])
    inclade = inclade-blacklist
    outclade = seqnames - inclade - blacklist

    print("GCA_001865635" in inclade)
    print("GCA_000262165" in outclade)

    sys.stderr.write("Num inclade sequences: %s\n" % len(inclade))
    sys.stderr.write("Num outclade sequences: %s\n" % len(outclade))
    args.uuid = str(uuid4())

    print(args.uuid)
    sp.call("dsk -kmer-size %(kmer_size)s -file %(msa)s -out %(uuid)s > /dev/null 2> /dev/null" % vars(args),shell=True)
    sp.call("dsk2ascii -file %(uuid)s.h5 -out %(uuid)s.txt > /dev/null 2> /dev/null" % vars(args),shell=True)

    kmers = set()
    for l in tqdm(open(args.uuid+".txt")):
        kmers.add(l.strip().split()[0])

    with open(args.out+".kmers.txt", "w") as O:
        for kmer in tqdm(kmers):
            rkmer = revcom(kmer)
            inclade_alleles = []
            outclade_alleles = []
            for seq in inclade:
                inclade_alleles.append(1 if (kmer in fasta_dict[seq] or rkmer in fasta_dict[seq]) else 0)
            for seq in outclade:
                outclade_alleles.append(1 if (kmer in fasta_dict[seq] or rkmer in fasta_dict[seq]) else 0)

            if sum(inclade_alleles)==len(inclade_alleles) and sum(outclade_alleles)==0:
                print(kmer)
                O.write("%s\t%s\n" % (kmer,args.clade))


            # if sum(inclade_alleles)/len(inclade_alleles)>0.9 and sum(outclade_alleles)/len(outclade_alleles)<0.1:
            #     print([list(inclade)[i] for i in range(len(inclade_alleles)) if inclade_alleles[i]==0])
            #     print([list(outclade)[i] for i in range(len(outclade_alleles)) if outclade_alleles[i]==1])
            #     import pdb; pdb.set_trace()
            # if "GTACCGTCACTTTCGCTTCGTCCCTACT" in kmer:
            #     import pdb; pdb.set_trace()


    os.remove(args.uuid+".h5")
    os.remove(args.uuid+".txt")




parser = argparse.ArgumentParser(description='XXX pipeline',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--msa',help='VCF file',required=True)
parser.add_argument('--out',help='VCF file',required=True)
parser.add_argument('--clade-csv',help='VCF file',required=True)
parser.add_argument('--clade',help='VCF file',required=True)
parser.add_argument('--kmer-size',default=31,type=int,help='VCF file')
parser.add_argument('--blacklist',help='VCF file')
parser.set_defaults(func=main)
args = parser.parse_args()
args.func(args)
