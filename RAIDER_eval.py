#!/software/python/3.3.3/bin/python3.3
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
    parser.add_argument('-r', '--raider', action = "store_true", help = "Run raider", default = True)
    parser.add_argument('-f', type = int, help = "E.R. occurrence threshold", default = 2)
    parser.add_argument('-s', '--seed', help = "spaced seed string", default = "1111011110111101111")
    parser.add_argument('-d', '--output_dir', help = "Raider output directory", default = None)
    parser.add_argument('-e', '--output_ext', help = "Output Extension", default = None)
    parser.add_argument('-C', '--cleanup_off', dest = "cleanup", action = "store_false", help = "Turn off file cleanup", default = True)
    # REPEAT MASKER ARGUMENTS
    parser.add_argument('--masker_dir', help = "Repeat masker output directory", default = None)
    parser.add_argument('-p', '--pa', type = int, help = "Number of processors will be using")
   
    subparsers = parser.add_subparsers(dest="subparser_name")
    parser_seqs = subparsers.add_parser("seq_files")
    parser_seqs.add_argument('seq_files', nargs = "+", help = "Files containing genomic sequence")
    parser_chrom = subparsers.add_parser("chrom_sim")
    parser_chrom.add_argument('chromosome', help = "Template chromosome file")
    parser_chrom.add_argument('repeat', help = "Repeat file")
    parser_chrom.add_argument('num_sims', type = int, help ="Number of simulations")
    parser_chrom.add_argument('--rng_seed', type = int, help = "RNG seed", default = None) 
    parser_chrom.add_argument('-l', '--length', type = int, help = "Simulated sequence length", default = None)
    parser_chrom.add_argument('-n', '--negative_strand', action = "store_true", help = "Use repeats on negative string", default = False)
    parser_chrom.add_argument('-f', '--family_file', help = "List of repeat families to use", default = None)
    parser_chrom.add_argument('-o', '--output', help = "Output file (Default: replace chromosome file \".fa\" with \".sim.fa\")")
    parser_chrom.add_argument('--chrom_dir', help = "Simulated chromosome output directory", default = None) 
    return parser.parse_args(args)

def simulate_chromosome(chromosome, repeat, rng_seed, length, neg_strand, fam_file, output_file, file_index, chrom_dir):
    """Given chromosome file and repeat file and rng_seed, runs chromosome 
    simulator and then passes raider params (including new simulated chromosome 
    file) into run_raider"""
    #output_file = output_file if output_file else tempfile.mkstemp(dir = ".")[1]
    output_file = output_file if output_file else re.sub(".fa$", "."+str(file_index)+".sim.fa", chromosome) 
    output_part = "-o %s " % (output_file)
    seed_part = "-s %d " % (rng_seed) if rng_seed else ""
    length_part = "-l %d " % (length) if length else ""
    neg_strand_part = "-n " if neg_strand else ""
    fam_file_part = "-f %s " % (fam_file) if fam_file else ""
    cmd = "./chromsome_simulator.py " + output_part + seed_part + length_part + neg_strand_part + fam_file_part + "%s %s" % (chromosome, repeat)
    
    print(cmd)
    p = pbsJobHandler(batch_file = "%s.batch" % (output_file), executable = cmd)
    p.submit()
    p.chrom_output = output_file
    p.index = file_index
    return p

def run_raider_chrom(p, seed, f, output_dir):
    output_dir = output_dir if output_dir else "raider_output_%d" % (p.index)
    
    p.wait();
    return run_raider(seed, f, p.chrom_output, output_dir)


def run_raider(seed, f, input_file, output_dir):
    """Given raider parameters and an input file, run RAIDER and put the output into
    the directory specified in output_dir (creating a random name is none is
    specified."""
    output_dir = output_dir if output_dir else tempfile.mkdtemp(prefix = "raider_output", dir = ".")

    if not os.path.exists('./%s' % (output_dir)):
        os.makedirs('./%s' % (output_dir))

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
    p2.seq_file = p.file
    p2.output = output
    return p2

def run_repeat_masker(p, num_processors, masker_dir):
    p.wait()
    if masker_dir and not os.path.exists('./%s' % (masker_dir)):
        os.makedirs('./%s' % (masker_dir))
    lib_output = re.sub("((\.fa)|(\.fasta))$", ".lib.fa" , p.output)
    library_part = "-lib %s " % (lib_output)
    processor_part = "-pa %d " % (num_processors) if num_processors else ""
    output_part = "-dir %s " % (masker_dir) if masker_dir else ""
    cmd = "RepeatMasker " + library_part + processor_part + output_part + "%s" % (p.seq_file)
    print(cmd)
    p2 = pbsJobHandler(batch_file = "repeatmasker", executable = cmd, RHmodules = ["RepeatMasker", "python-3.3.3"]) 
    p2.submit()
    p2.masker_output = masker_dir
    return p2

if __name__:
    args = parse_params(sys.argv[1:])
    
    if args.subparser_name == "seq_files" and args.seq_files:	
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
        J3 = []
        for p in J2:
            J3.append(run_repeat_masker(p, args.pa, args.masker_dir))
        for p in J3:
            p.wait()

    elif args.subparser_name == "chrom_sim" and args.chromosome and args.repeat and args.num_sims:
        ### Create n simulated chromosomes, invoke RAIDER on each chromosome
        ### and put the resulting pbs object into the J list.
        J =[]

        if args.chrom_dir and not os.path.exists('./%s' % (args.chrom_dir)):
            os.makedirs('./%s' % (args.chrom_dir))


        for i in range(int(args.num_sims)):
            J.append(simulate_chromosome(chromosome = args.chromosome, repeat = args.repeat, \
                     rng_seed = args.rng_seed, length = args.length, neg_strand = args.negative_strand, \
                     fam_file = args.family_file, output_file= args.output, file_index = i, chrom_dir = args.chrom_dir))
        J2 = []
        for p in J:
            J2.append(run_raider_chrom(p, seed = args.seed, f = args.f, output_dir = args.output_dir))
        J3 = []
        for p in J2:
            output = re.sub("((\.fa)|(\.fasta))$", ".consensus%sfa" % ("." if not args.output_ext else "." + output_ext + "."), p.file)
            J3.append(create_raider_consensus(p, output))
        J4 = []
        for p in J3:
            J4.append(run_repeat_masker(p, args.pa, args.masker_dir))
        for p in J4:
            p.wait()





        
