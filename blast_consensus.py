from Bio.Blast.Applications import NcbiblastxCommandline
import sys
import re
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import subprocess
import os.path
from coverageBLAST import coverageBLAST

"""
BLAST the consensus sequences against the known cloations of repeats to calculate coverage.
Input:
* consensus sequence file
* sequence file
* .fa.out file
"""
isRerun = False    # Only set True for debugging purposes

def create_query(consensus_file): 
    L = []
    for line in open(consensus_file):
        id, s = re.split("\s+", line.rstrip())
        L.append(SeqRecord(seq = Seq(s), id = 'er.'+id, description = "elementry repeat"))
    SeqIO.write(L, consensus_file + ".fa", 'fasta')
    return consensus_file + ".fa"

def create_target(seq_file, rm_file):
    assert rm_file.endswith(".fa.out")
    target_file = re.sub(".fa.out", ".out.fa", rm_file)

    D = {}
    if not os.path.isfile(target_file):

        chr_seq = SeqIO.read(seq_file, 'fasta').seq

        L = []
        fp = open(rm_file)
        fp.readline()
        fp.readline()
        for line in fp:
            if line.rstrip():
                A = re.split("\s+", line.strip())
                start = int(A[5])
                stop = int(A[6]) - 1
                location = "{chr}:{start}-{stop}".format(chr=A[4], start=start, stop=stop)
                id = "{class_name}|{family}|{location}".format(class_name=A[9], family=A[10], location=location)
                D[id] = 0 if id not in D else D[id] + 1
                id += " " + str(D[id])
                description = ""
                R = SeqRecord(seq = chr_seq[start:stop], id = id, description=description)
                L.append(R)
        SeqIO.write(L, target_file, "fasta")

    return target_file, SeqIO.parse(target_file, 'fasta')


    return NCBIXML.parse(open(query_file + ".xml"))

                                         
def main(consensus_file, seq_file, rm_file, e_val = 0.0001):
    #consensus_file = create_query(consensus_file)

    # Consensus_records is a generator of the consensus sequence records
    consensus_records = SeqIO.parse(consensus_file, 'fasta')
    CSD = {r.id:str(r.seq).upper() for r in consensus_records}


    # target_file is assignd a file name; Rpts is a sequene record generator of the known repeats
    target_file = "seq.out.fa"
    target_file, Rpts = create_target(seq_file, rm_file);
    RPD = {r.id:str(r.seq).upper() for r in Rpts}
    

    # Run blast
    qc, tc = coverageBLAST(consensus_file, target_file, e_val, isRerun = isRerun);

    print("qc:\n")
    sum, total = 0,0
    for k,v in qc.items():
        f = v[0]/v[1] if v[1] > 0 else -1
        sum += f
        total += 1
        print("{er}: {covered}, {length}, {frac}".format(er = k, covered = v[0], length = v[1], frac = round(f,4)))
    print("Average coverage: {avg}".format(avg = round(sum / total, 4)))

    print("-"*50)
    print("tc:\n")
    D = {}
    for k,v in sorted(tc.items()):
        fam = k[:k.find("|")]
        f = v[0]/v[1] if v[1] > 0 else -1
        if fam not in D:
            D[fam] = (0,0)
        print("{rpt}: {covered}, {length}, {frac}".format(rpt = k, covered = v[0], length = v[1], frac = round(f,4)))
        D[fam] = D[fam][0] + f, D[fam][1] + 1

    print("-"*50)
    for fam,T in sorted(D.items()):
        print("{fam}: {avg}".format(fam=fam, avg=round(T[0]/T[1] if T[1] > 0 else -1,4)))
    


if __name__ == "__main__":
    consensus_file = sys.argv[1]
    seq_file = sys.argv[2]
    rm_file = sys.argv[3]
                
    main(consensus_file, seq_file, rm_file);

