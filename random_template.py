import sys
import re
import random
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

bases = "CGT"
def random_sequence(j):
    return "".join([random.choice(bases) for i in range(j)])

l = int(sys.argv[1])    # Total length of simulatd genome
n = int(sys.argv[2])    # A repeat will start at every coordinate c*n
m = int(sys.argv[3])    # The repeat will be m bases long

output_file = sys.argv[4]

header = """   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

"""
sample_line="  0  0  0.0  0.0  chr_sim     %d %d (%d) +  X           SINE/X                0  0  (0)      0\n"


rpt_seq = random_sequence(m)

fp = open(output_file + ".fa.out", "w")
fp.write(header)

s = random_sequence(n)
for i in range(n, l-m, n):
    s = s + rpt_seq + random_sequence(n-m)
    fp.write(sample_line % (int(i+1), int(i+m), int(l-(i+m))))

if len(s) < l:
    s += random_sequence(l - len(s))

SeqIO.write([SeqRecord(seq = Seq(s), id = "Simulated template", description = "l = %d, n=%d, m=%d" % (l,n,m))], output_file + ".fa", 'fasta')
