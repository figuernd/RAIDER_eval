import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
import sys
import markov_gen
import re

def parse_params(args):
    parser = argparse.ArgumentParser(description = "Generate simulated chromsome")
    parser.add_argument('-c', '--cutoff', type = int, help = "Limit model to first c non-N bases")
    parser.add_argument('-k', type = int, help = "Order of Markov chain", default = 5)
    parser.add_argument('-s', '--seed', type = int, help = "RNG seed", default = None)
    parser.add_argument('-l', '--length', type = int, help = "Simulated sequence length", default = None)
    parser.add_argument('-n', '--negative_strand', action = "store_true", help = "Use repeats on negative string", default = True)
    parser.add_argument('-f', '--family_file', help = "List of repeat families to use", default = None)
    parser.add_argument('-o', '--output', help = "Output file (Default: replace chomosome file \".fa\" with \".sim.fa\")")
    parser.add_argument("chromsome", help = "Chromosome file")
    parser.add_argument("repeat_file", help = "Repeat file")
    return parser.parse_args(args)


def nextRepeat(rpt_file, use_negative, S):
    fp = open(rpt_file)
    fp.readline()
    fp.readline()
    for line in fp:
        if line.rstrip():
            A = re.split("\s+", line.strip())
            chr, start, finish, strand, family = A[4], int(A[5])-1, int(A[6]), A[8], A[9]
            if (strand == '+' or use_negative) and (family in S or not S):
                yield chr, start, finish, strand, family


def generate_chromosome(seq_file, rpt_file, k, output, length = None, negative_strand = False, family_file = None, seed = None, cutoff = None):

    original_sequence = "".join([str(r) for r in SeqIO.read(seq_file, 'fasta')])
    if cutoff:
        i = 0
        while i < len(original_sequence) and original_sequence[i] not in "ACGTacgt":
            i += 1
        j = i + 1;
        count = 1;
        while j < len(original_sequence) and count < cutoff:
            if original_sequence[j] in "ACGTacgt":
                count += 1
            j += 1
        original_sequence = original_sequence[i:j]
            
    
    if not length:
        length = len(original_sequence)

    simulated_sequence = markov_gen.generate_sequence(original_sequence, k, length, seed)

    S = {x.rstrip() for x in open(family_file)} if family_file else {}
    for chr, start, finish, strand, family in nextRepeat(rpt_file, negative_strand, S):
        if finish < length:
            simulated_sequence = simulated_sequence[:start] + original_sequence[start:finish] + simulated_sequence[finish:]

    output_file = ouput if output else re.sub(".fa$", ".sim.fa", seq_file)
    SeqIO.write([SeqRecord(seq = Seq(simulated_sequence), id = re.sub(".fa", "", seq_file), description = "%d-order Markov chain simulation" % (k))], output_file, 'fasta')


if __name__ == "__main__":
    args = parse_params(sys.argv[1:])
    generate_chromosome(args.chromsome, args.repeat_file, args.k, args.output, args.length, args.negative_strand, args.family_file, args.seed, args.cutoff)
