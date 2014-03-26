#!/usr/bin/python2.7
############################################################################
# simulation.py
# Code for generating a simulated chromsome based on a k-degree markov chain and a model sequences
# Written by: Jiajun Wang, John Karro
# Date: March, 2014
# Email: karroje@miamiOH.edu

import re
import sys
import argparse
import random
import unittest
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

fn={
    'A':0, 'a':0,
    'C':1, 'c':1,
    'G':2, 'g':2,
    'T':3, 't':3,
    'N':-1, 'n':-1
}

bases = "acgtACGT"
baseSet = set(bases)
base2index = {x:y for x,y in zip("acgtnACGTN", [0,1,2,3,-1,0,1,2,3,-1])}

def ProbTuple(Mark,kmo):
    total = 0
    tlist=[0]*4
    for t in range(4):
        total += Mark[(kmo<<2)|t]
    if total == 0 :
        return [0,0,0,0]     
    probsum = 0
    for t in range(4):
        tlist[t] = probsum + Mark[(kmo<<2)|t]/float(total)
        probsum = tlist[t]
    return tlist
    
def pickFromProbTuple(probTuple):
    """Pick a random index based on a probability tuple.  Assumes last element of tuple is 1."""
    r = random.random()
    i = 0
    while r > probTuple[i]:
        i += 1
    return i

#return the index of the first base without 'N' in following it for k length 
def NextIndex(i,k,seq): # Should be having inSeq as a parameter
    """Return the smallest j >= i such that inSeq[j:j+k] contains only bases"""
    j = i
    while j + k < len(seq) and seq[j] not in baseSet:
        j += 1
    if k == 0:
        return j
    while j + k < len(seq):
        p = 1
        while p < k and seq[j+p] in baseSet:
            p += 1
        if p == k:
            return j
        j += p + 1
    return len(seq)

   
def Markov(k,inSeq): # inSeq as parameter?
    """Create Markov chain
    Input: k is the degree of the Markov chain
    """
    Mark= [0]*(4**(k+1))   
    mask=(1<<((k+1)*2))-1
    t=0
    #get first k+1 subsequence
    i=NextIndex(0,k,inSeq)
    for j in range(k+1):
        t=(t<<2) | base2index[inSeq[i+j]] 

    Mark[t] = 1
    #go through the whole seq
    i += k+1     
    while i<len(inSeq) : 
        if inSeq[i] not in bases:
             i=NextIndex(i,k,inSeq)
             if (i+k)>=len(inSeq):
                 break
             t=0
             for j in range(k+1):
                 t = (t<<2) | base2index[inSeq[i+j]]
             i += k+1 
        else:
             t = ((t<<2) | base2index[inSeq[i]]) & mask     
             i += 1             
        Mark[t] += 1   
        
    
    Klist=[[0,0,0,0]]*(4**k)  
    for KMOne in range(4**k):         
        Klist[KMOne] = ProbTuple(Mark,KMOne)
    return Klist       




def generate_sequence(seq, k, l = None, rng_seed = None):
    """Generate a random sequence of length l using a k-order Markov chain modeled on sequence seq"""
    random.seed(rng_seed)
    seq_len = l if l else len(seq)

    ## Generate order Markov chains
    markov_list = [Markov(i, seq) for i in range(k+1)]

    mask = (1 << (k*2)) - 1

    j = 0
    seq = ['']*seq_len
    key = 0
    for i in range(0, seq_len):
        next_tuple = markov_list[j][key]
        if next_tuple[-1] < 0.9999999:
            j = 0
            key = 0
            next_tuple = markov_list[j][key]
        base = pickFromProbTuple(next_tuple)
        key = ((key << 2) | base) & mask
        seq[i] = bases[base]
        j = min(j+1,k)

    return "".join(seq)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Generate simulated sequence using a k-order Markov chain")
    parser.add_argument('-k', type = int, help = "Order of Markov chain", default = 5)
    parser.add_argument('-l', type = int, help = "Sequence length", default = None)
    parser.add_argument('-r', '--rng_seed', type = int, help = "rng seed", default = None)
    parser.add_argument('-o', '--output', help = "Output file (stdout by default -- fasta format if file name ends in .fa or .fasta)", default = "-")
    parser.add_argument('seq_file', help = "Model sequence file (fasta format)")
    args = parser.parse_args()

    s = "".join([str(r) for r in SeqIO.read(args.seq_file, 'fasta')])
    s2 = generate_sequence(s, args.k, args.l, args.rng_seed)
    if args.output == '-':
        sys.stdout.write(s2)
    elif re.search(".fa(asta)?$", args.output):
        SeqIO.write(SeqRecord(seq = Seq(s2), id = args.seq_file, description = "simulation"), args.output, "fasta")
    else:
        open(args.output, "w").write(s2)

                
