import argparse
from Bio import SeqIO
from Bio import Seq
import random
import sys
import markov_gen
import re

def parse_params(args):
    parser = argparse.ArgumentParser(description = "Generate simulated chromsome")
    parser.add_argument('-k', type = int, help = "Order of Markov chain", default = 5)
    parser.add_argument('-s', '--seed', type = int, help = "RNG seed", default = None)
    parser.add_argument('-l', '--length', type = int, help = "Simulated sequence length", default = None)
    parser.add_argument("chromsome", help = "Chromosome file")
    parser.add_argument("repeat_file", help = "Repeat file")
    return parser.parse_args(args)


def nextRepeat(rpt_file):
    fp = open(rpt_file)
    fp.readline()
    fp.readline()
    fp.readline()
    for line in fp.readline():
        A = re.split(line.strip(), "\s+")
        yield(A[4], int(A[5]), int(A[6]), A[8], A[9])


def (seq_file, rpt_file, k, length = None):

    s = "".join([str(r) for r in SeqIO.read("seq_file", 'fasta')])
    simulate_sequence = Markov.generate_sequence(s, k, length)

    open(rpt_file)
