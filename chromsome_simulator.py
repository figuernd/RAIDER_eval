import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
import sys
import re
import os.path

import markov_gen


def parse_params(args):
    parser = argparse.ArgumentParser(description = "Generate simulated chromsome")
    # parser.add_argument('-c', '--cutoff', type = int, help = "Limit model to first c non-N bases")
    parser.add_argument('-l', '--length', type = int, help = "Simulated sequence length", default = None)
    parser.add_argument('-k', type = int, help = "Order of Markov chain", default = 5)
    parser.add_argument('-s', '--seed', type = int, help = "RNG seed", default = None)
    parser.add_argument('-n', '--negative_strand', action = "store_true", help = "Use repeats on negative string", default = False)
    parser.add_argument('-f', '--family_file', help = "List of repeat families to use", default = None)
    parser.add_argument('-m', '--mask', action = "store_true", help = "Turn masking on (all repeats printed as lower case).", default = False)
    parser.add_argument('-S', '--suppress_pmck', action = "store_true", help = "Suppress the generation of a .pmc<k> file to store the markov chain for re-use")
    parser.add_argument('--mi', '--max_interval', dest = "max_interval", type = int, help = "Maximum allowed length of interval between repeats; -1 value (default) means no maximum", default = -1)
    #parser.add_argument('-o', '--output', help = "Output file (Default: replace chomosome file \".fa\" with \".sim.fa\")")
    parser.add_argument("seq_file", help = "Sequence file (must be .fa)")
    parser.add_argument("repeat_file", help = "RepeatMasker file (.fa.out)")
    parser.add_argument("output", help = "Output file")
    return parser.parse_args(args)


def nextRepeat(rpt_file, use_negative = True, S = {}):
    """Generator: each invokation returns the chromosome, start, finish, starand, 
    and family for the next repeat of the repeatmasker .fa.out files.  S, if not empty,
    is a filter for which repeats to use."""
    fp = open(rpt_file)
    fp.readline()
    fp.readline()
    for line in fp:
        if line.rstrip():
            A = re.split("\s+", line.strip())
            chr, start, finish, strand, family, rpt_class = A[4], int(A[5])-1, int(A[6]), A[8], A[9], A[10]
            if (strand == '+' or use_negative) and (family in S or not S):
                yield chr, start, finish, strand, family, rpt_class

# fa_out_header: The fixed header lines for the .fa.out file
fa_out_header = "\tSW\tperc\tperc\tperc\tquery\tposition in query\tmatching\trepeat\tposition in  repeat\nscore\tdiv.\tdel.\tins.\tsequence\tbegin\tend\t(left)\trepeat\tclass/family\tbegin\tend (left)\tID\n"
# fa_out_template: A template for creating lines for the .fa.out file.
fa_out_template = "\t0\t0\t0\t0\t{chr}\t{start}\t{finish}\t(0)\t{strand}\t{family}\t{rpt_class}\t0\t0\t(0)\t1"
def generate_chromosome(seq, markov_list, rpt_gen, mask = False, length = None, max_interval = None):
    """
    Generate a syntehtic sequence with real repeats:
    * seq: A sequence (as a string).
    * markov_list: List of the k+1 i-th order markov chains (from the markov_gen module).
    * rpt_gen: A generating function returning the repeat information (created by nextRepeat)
    * mask: If true, all repeats will be lower-case.  Otherwise, upper case.)
    """
    coord_adjust = 0
    current_coord = 0
    if max_interval == -1:
        max_interval = len(seq)
    s = []                # Hold the sequence (in chunks)
    fa_out = [fa_out_header]           # Hold the new .fa.out file contents (by line)
    if not length:
        length = len(seq)

    for chr, start, finish, strand, family, rpt_class in rpt_gen:
        if start > length:
            break
        if start >= current_coord:
            s.append(markov_gen.generate_sequence(markov_list, min(start - current_coord, max_interval)))
            coord_adjust += max(0, start-current_coord-max_interval)

            s.append(seq[start:finish].lower() if mask else seq[start:finish].upper())
            fa_out.append(fa_out_template.format(chr=chr, start=start+1-coord_adjust, finish=finish-coord_adjust, strand=strand, family=family, rpt_class=rpt_class))
        current_coord = finish
    if length > current_coord:
        s.append(markov_gen.generate_sequence(markov_list, length - current_coord))
    return "".join(s), "\n".join(fa_out)

def loadSeqAndChain(seq_file, k, suppress_save = False):
    """Load the sequence and the Markov Chain List.
    Load the MC list from a file if it exists.  If not, create the chain
    and save it to the file for the next use (skip the save if suppressed)."""
    template_seq = str(SeqIO.read(seq_file, 'fasta').seq)
    
    mc_file = re.sub("\.(fa|fasta)$", ".pmc%d" % (k), seq_file)
    if os.path.exists(mc_file):
        return template_seq, markov_gen.read_pmck(mc_file)
    else:
        markov_list = markov_gen.MarkovArray(k, template_seq)
        if not suppress_save:
            markov_gen.pickle_markov_list(markov_list, mc_file)
        return template_seq, markov_list

def create_chromosome_file(seq_file, repeat_file, output_file, k = 5, use_3prime = True, filter_file = "rpt_list.txt", mask = False, seed = None, length = None, suppress = False, max_interval = -1):
    """
    Create a simualted chrosome with real repeat sequences from a chromsoe file.
    Parameters:
    * seq_file: fasta <seq>.fa, file containing the template sequence.
      -- Assumed to exist a file <seq>.fa.out containing the repeatmasker annotations.
    * k: Use a k-order markov chain.  There must exists a markov chain file <seq>.pmc<k>.
    * output_file: Fasta file to print sequence to.
    * use_3prime: If false, only sequence on the 5' strand will be used.  Default: True
    * filter_file: A list of the repeats that should be used.  If empty: all repeats.  Default: "rpt_list.txt"
    * mask: If true: copied repeats will be in lower case.  Default: False
    * seed: RNG seed
    """
    random.seed(args.seed)
    template_seq, markov_list = loadSeqAndChain(args.seq_file, args.k, suppress)
    filter_set = {y.strip() for line in open(filter_file) for y in re.split("\s+", line.rstrip())} if filter_file else {}
    rpt_gen = nextRepeat(repeat_file, use_3prime, filter_set)
    simulated_sequence, fa_out = generate_chromosome(template_seq, markov_list, rpt_gen, mask, length, max_interval)
            
    SeqIO.write([SeqRecord(seq = Seq(simulated_sequence), id = "seq_file", description = "Simulated sequence from %s using order %d markov chain" % (seq_file, len(markov_list)-1))], output_file, 'fasta')
    open(output_file + ".out", "w").write(fa_out)

if __name__ == "__main__":
    args = parse_params(sys.argv[1:])
    create_chromosome_file(seq_file = args.seq_file, k = args.k, output_file = args.output, 
                           repeat_file = args.repeat_file, use_3prime = args.negative_strand, 
                           filter_file = args.family_file, mask = args.mask, seed = args.seed,
                           length = args.length, max_interval = args.max_interval)
