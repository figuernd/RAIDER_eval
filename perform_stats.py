#!/software/python/3.3.3/bin/python3.3

############################################################################

import argparse
import sys
import re
import Bio
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_repeat_file(rep_file):
    rf = open(rep_file, 'r')
    rf.readline()
    rf.readline()
    rf.readline()
    indices = []
    while True:
        line = rf.readline()
        if not line:
            break
        A = re.split("\s+", line)
        indices.append((int(A[6]) - 1, int(A[7]) - 1, A[10]))
    return indices

def get_real_repeats(real_repeats, seq_file):
    info = read_repeat_file(real_repeats)
    reps = []
    record = SeqIO.read(seq_file, "fasta")
    for i in info:
        reps.append(record.seq[i[0]:i[1]])
    return reps

def get_generated_repeats(masker_output, consensus_output):
    info = read_repeat_file(masker_output)
    reps = []
    record_dict = SeqIO.to_dict(SeqIO.parse(consensus_output, "fasta"))
    for i in info:
        reps.append(record_dict[i[2]].seq)
    return reps
    
def get_stats(real_reps, gen_reps):
    tps = []
    fns = []
    for i in range(len(real_reps)):
        matched = False
        for j in range(len(gen_reps)):
            if str(gen_reps[j]) == str(real_reps[i]):
                matched = True
                gen_reps.pop(j)
                tps.append(str(real_reps[i]))   
                break
        if matched == False:
            fns.append(str(real_reps[i]))
    fps = []
    for i in gen_reps:
        fps.append(str(i))

    
    tp = len(tps)
    fp = len(fps)
    fn = len(fns)
    tpr = tp/float(tp + fn)
    ppv = tp/float(tp + fp)
    fdr = 1 - ppv
    return tp, fp, fn, tpr, ppv, fdr, tps, fps, fns


def main(seq_file, real_repeats, consensus_output, masker_output, output_file, print_reps):
    # assuming masker_repeats references repeats named in consensus_output file
    real_reps = get_real_repeats(real_repeats, seq_file)
    gen_reps = get_generated_repeats(masker_output, consensus_output)
    pos = len(gen_reps)
    tp, fp, fn, tpr, ppv, fdr, tps, fps, fns= get_stats(real_reps, gen_reps)
    f = open(output_file, 'w')
    f.write("TP: %d \n" % (tp))
    f.write("FP: %d \n" % (fp))
    f.write("FN: %d \n" % (fn))
    f.write("TPR: %f \n" % (tpr))
    f.write("PPV: %f \n" % (ppv))
    f.write("FDR: %f \n" % (fdr))
    f.write("\n")
    if print_reps:
        f.write("\nThese repeats were correctly identified (true positives):\n")
        for rep in tps:
            f.write("\t %s \n" % (rep))
        f.write("\nThese repeats were incorrectly identified (false positives):\n")
        for rep in fps:
            f.write("\t %s \n" % (rep))
        f.write("\nThese repeats were missed (false negatives):\n")
        for rep in fns:
            f.write("\t %s \n" % (rep))
    f.close()    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Generate Sensitivity and Specificity Stats")
    parser.add_argument("--print_reps", action = "store_true", help = "Print out the repeats", default= False)
    parser.add_argument("seq_file", help = "Sequence file")
    parser.add_argument("repeat_file", help = "Repeats file")
    parser.add_argument("consensus", help = "Consensus Sequence for the Sequence")
    parser.add_argument("masker_output", help = "Masker output using consensus sequence and sequence file")
    parser.add_argument("output_file", help = "Statistics output file")
    args = parser.parse_args()
    main(args.seq_file, args.repeat_file, args.consensus, args.masker_output, args.output_file, args.print_reps)

    


