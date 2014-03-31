#!/usr/bin/python3.3
import sys
import subprocess
import os.path
import os.path
import argparse
from redhawk import *
import tempfile
import re

def parse_params(args):
    """Parse command line arguments using the argparse library"""
    parser = argparse.ArgumentParser(description = "Evaluate RAIDER against RepeatScout")
    parser.add_argument('-f', type = int, help = "e.r. occurance threshold", default = 2)
    parser.add_argument('-s', '--seed', help = "spaced seed string", default = "1111011110111101111")
    parser.add_argument('-r', '--raider', action = "store_true", help = "Run raider", default = True)
    parser.add_argument('-d', '--output_dir', help = "Raider output directory", default = None)
    parser.add_argument('-e', '--output_ext', help = "Output Extension", default = None)
    parser.add_argument('-C', '--cleanup_off', dest = "cleanup", action = "store_false", help = "Turn off file cleanup", default = True)
    parser.add_argument('seq_files', nargs = "+", help = "Files containing genomic sequence")
    return parser.parse_args(args)
    
def run_raider(seed, f, input_file, output_dir = None):
    """Given raider parameters and an input file, run RAIDER and put the output into
    the directory specified in output_dir (creating a random name is none is
    specified."""
    output_dir = output_dir if output_dir else tempfile.mkdtemp(dir = ".")

    if not os.path.exists('./%s' % (args.output_dir)):
        os.makedirs('./%s' % (args.output_dir))

    cmd = "./raider -q -c %d %s %s %s" % (f, seed, input_file, output_dir)
    print(cmd)
    p = pbsJobHandler(batch_file = "raider", executable = cmd)
    p.submit()
    p.file = input_file
    p.raider_output = output_dir
    return p

def create_raider_consensus(p, output):
    """Given the pbs object (from redhawk.py) used to start a RAIDER job, this
    waits until the job is done, then invokes consensus_seq.py on the ouput and
    writes the results to the output directory."""
    p.wait();
    cmd = "./consensus_seq.py -s %s -e %s %s" % (p.file, p.raider_output + "/elements", output)
    print(cmd)
    p2 = pbsJobHandler(batch_file = "%s.batch" % (output), executable = cmd)
    p2.submit()
    p2.output = output
    return p2


if __name__:
    args = parse_params(sys.argv[1:])

    ### For each listed file (in args.seq_files): invoke RAIDER on the file
    ### and put the resuling pbs object ito the J list.
    J = []
    for file in args.seq_files:
        assert file.endswith(".fa") or file.endswith('.fasta')
        J.append(run_raider(seed = args.seed, f = args.f, input_file = file, output_dir = args.output_dir))
        
    ### Sequentially work the the J list and, when the next job has finished,
    ### start the conensus sequence job running and stick the new pbs job nto
    ### the J2 job. 
    J2 = []
    for p in J:
              # The velow line is creating the output name in such a way as to avoid conflict.
              # It will change X.fa to X.conensus.fa (if args.output_ext == None), stick that output_ext
              # string in there if it does exist.
        output = re.sub("((\.fa)|(\.fasta))$", ".consensus%sfa" % ("." if not args.output_ext else "." + output_ext + "."), p.file)
        J2.append(create_raider_consensus(p, output))

    for p in J2:
        p.wait()






        
