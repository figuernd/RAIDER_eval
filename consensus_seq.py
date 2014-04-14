#!/software/python/3.3.3/bin/python3.3

############################################################################
# consensus_seq.py
# program for generating consensus sequences
# Written by: Rachael Morton

import random
import argparse
import unittest
import sys
import re
import Bio
from Bio import Seq
from Bio import SeqIO
from collections import defaultdict

def generate_dictionary():
    return {x:0 for x in ['A', 'C', 'G', 'T', 'N']}

def get_consensus(family_num, R, genome_string, wp):
    family_length = R[0][1] - R[0][0]
    C = [generate_dictionary() for i in range(family_length)]
    for start, finish in R:
        assert finish - start == family_length
        for i in range(family_length):
            try:
                C[i][genome_string[start+i].upper()] += 1
            except:

                C[i]['N'] += 1

    consensus = [max([(c,b) for b,c in D.items()])[1] for D in C]
    wp.write("%d\t%s\n" % (family_num, "".join(consensus)))

def main(seq,elements,output):  
    genomeFile = seq
    genome_string = next(SeqIO.parse(genomeFile, "fasta"))
    
    wp = open(output, "w")

    # Opening elements file
    elementFile = elements
    f = open(elementFile, 'r')
    f.readline()
    current_family = 0
    R = []
    while True:
        line = f.readline()
        if not line:
            get_consensus(current_family, R, genome_string, wp)
            break
        A = re.split("\s+", line)
        if int(A[0]) != current_family:
            get_consensus(current_family, R, genome_string, wp)
            R = []
            current_family = int(A[0])
        R.append((int(A[-2]), int(A[-1])))

        
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = "Generate Simulated Sequence")
    parser.add_argument(
        '-s', '--seq', action='store', required=False, type=str, default="../hg19.chr22.fa",
        help="sequence fasta file"
    )
    parser.add_argument(
        '-e', '--elements',action='store',required=False,type=str,
        default="../out.d/elements",
        help = "elements file"
    )
    parser.add_argument("output", type=str, help="output file")
    args = parser.parse_args()
    main(args.seq,args.elements,args.output)

