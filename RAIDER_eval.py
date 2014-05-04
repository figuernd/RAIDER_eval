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
    
    # GENERAL ARGUMENTS
    parser.add_argument('--organize', action = "store_true", help = "Create directory for all Raider Eval output", default = False)
    parser.add_argument('-r', '--raider', action = "store_true", help = "Run raider", default = True)
    
    # RAIDER ARGUMENTS
    parser.add_argument('-f', type = int, help = "E.R. occurrence threshold", default = 2)
    parser.add_argument('-m', '--min', type = int, help = "Minimum repeat length. Defaults to pattern length.", default = None)
    parser.add_argument('-s', '--seed', help = "Spaced seed string", default = "1111011110111101111")
    parser.add_argument('-d', '--output_dir', help = "Raider output directory", default = None)
    parser.add_argument('-e', '--output_ext', help = "Output Extension", default = None)
    parser.add_argument('-C', '--cleanup_off', dest = "cleanup", action = "store_false", help = "Turn off file cleanup", default = True)
    
    # REPEAT MASKER ARGUMENTS
    parser.add_argument('--masker_dir', help = "Repeat masker output directory", default = None)
    parser.add_argument('-p', '--pa', type = int, help = "Number of processors will be using")
    
    # STATISTICS ARGUMENTS
    parser.add_argument('--stats_dir', help = "Statistics output directory", default = None)
    parser.add_argument('--print_reps', action = "store_true", help = "Print out repeats in statistics file", default = False)

    subparsers = parser.add_subparsers(dest="subparser_name")
    
    # SEQUENCE FILE OPTION ARGUMENTS
    parser_seqs = subparsers.add_parser("seq_files")
    parser_seqs.add_argument('seq_files', nargs = "+", help = "Files containing genomic sequence")
    
    # CHROMOSOME SIMULATION OPTION ARGUMENTS
    parser_chrom = subparsers.add_parser("chrom_sim")
    parser_chrom.add_argument('chromosome', help = "Template chromosome file")
    parser_chrom.add_argument('repeat', help = "Repeat file")
    parser_chrom.add_argument('num_sims', type = int, help ="Number of simulations")
    parser_chrom.add_argument('--chrom_dir', help = "Simulated chromosome directory", default = None)
    parser_chrom.add_argument('--rng_seed', type = int, help = "RNG seed", default = None) 
    parser_chrom.add_argument('-l', '--length', type = int, help = "Simulated sequence length", default = None)
    parser_chrom.add_argument('-k', type = int, help = "Order of markov chain", default = 5)  # KARRO: Added this
    parser_chrom.add_argument('-n', '--negative_strand', action = "store_true", help = "Use repeats on negative string", default = False)
    parser_chrom.add_argument('-f', '--family_file', help = "List of repeat families to use", default = None)
    parser_chrom.add_argument('-o', '--output', help = "Output file (Default: replace chromosome file \".fa\" with \".sim.fa\")")
    return parser.parse_args(args)

def simulate_chromosome(chromosome, repeat, rng_seed, length, neg_strand, fam_file, chrom_dir, output_file, file_index, curr_dir, k):  # KARRO: Added k value
    """Given chromosome file and repeat file and rng_seed, runs chromosome 
    simulator and then passes raider params (including path to new simulated chromosome 
    file) into run_raider"""
    if chrom_dir:
        chrom_dir = "%s/%s" % (curr_dir, chrom_dir) if curr_dir else chrom_dir
        if not os.path.exists('./%s' % (chrom_dir)):
            os.makedirs('./%s' % (chrom_dir))


    # KARRO: I did a few things in the following:
    #        1) Modified the construction of cmd to reflect the changed in chromosome_simulator.py.
    #        2) Added the k_arg
    #        4) Rearranged the order of the *_arg assignments for readability
    #        5) Used the str.format() method to construct command -- improves readability, making it east to change later if necessary.

    # Output file is either specified or replace .fa with .sim.#.fa
    length_arg = "-l %d" % (length) if length else ""
    k_arg = "-k %d" % (k)
    seed_arg = "-s %d" % (rng_seed) if rng_seed else ""
    neg_arg = "-n" if neg_strand else ""
    fam_arg = "-f %s" % (fam_file) if fam_file else ""
    seq_arg = chromosome
    repeat_arg = repeat

    output_file = (output_file if output_file else re.sub(".fa$", ".sim.%d.fa"%(file_index), chromosome)) 
    output_path = "%s/%s" % (chrom_dir, output_file) if chrom_dir else "%s/%s" % (curr_dir, output_file) if curr_dir else output_file
    output_arg = "%s" % (output_path)

    cmd = "python3.3 chromosome_simulator.py {length} {k} {seed} {neg} {fam} {seq} {repeat} {output}".format(length=length_arg, k=k_arg, seed=seed_arg, neg=neg_arg, fam=fam_arg, seq=seq_arg, repeat=repeat_arg, output=output_arg)
    print(cmd)
    p = pbsJobHandler(batch_file = "%s.batch" % (output_file), executable = cmd)
    p.submit()
    p.chrom_output = output_path
    p.curr_dir = curr_dir
    p.index = file_index
    return p

def run_raider_chrom(p, seed, f, m, output_dir):
    """Given the pbs object used to start simulated chromosome job as well as 
    raider parameters, wait until the job is done and then call run_raider 
    on the newly generated chromosome output file"""
    p.wait();
    return run_raider(seed, f, m, p.chrom_output, output_dir, p.curr_dir)


def run_raider(seed, f, m, input_file, output_dir, curr_dir):
    """Given raider parameters and an input file, run RAIDER and put the output into
    the directory specified in output_dir (creating a random name is none is
    specified."""
    if output_dir:
        output_dir = "%s/%s" % (curr_dir, output_dir) if curr_dir else output_dir
    else:
        dir_part = curr_dir if curr_dir else "."
        output_dir = output_dir if output_dir else tempfile.mkdtemp(prefix = "raider_output", dir = dir_part)
    if not os.path.exists('./%s' % (output_dir)):
        os.makedirs('./%s' % (output_dir))
    min_arg = "-m %d" % (m) if m else ""
    cmd = "./raider -q -c %d %s %s %s %s" % (f, min_arg, seed, input_file, output_dir)
    print(cmd)
    p = pbsJobHandler(batch_file = "raider", executable = cmd)
    p.submit()
    p.file = input_file
    p.curr_dir = curr_dir
    p.raider_output = output_dir
    return p

def create_raider_consensus(p, output):
    """Given the pbs object (from redhawk.py) used to start a RAIDER job, this
    waits until the job is done, then invokes consensus_seq.py on the ouput and
    writes the results to the output directory."""
    p.wait();
    cmd = "python3.3 consensus_seq.py -s %s -e %s/elements %s" % (p.file, p.raider_output, output)
    print(cmd)
    p2 = pbsJobHandler(batch_file = "%s.batch" % (os.path.basename(output)), executable = cmd)
    p2.submit()
    p2.curr_dir = p.curr_dir
    p2.seq_file = p.file
    p2.output = output
    return p2

def run_repeat_masker(p, num_processors, masker_dir):
    """Given the pbs object used to start a consensus sequence job as well as
    repeatmasker arguments, wait until the job is done and then call repeatmasker 
    on the output and put results in masker_dir (current dir if unspecified)"""
    p.wait()
    # Determine actual path the masker directory (if specified)
    if masker_dir:
        masker_dir = "%s/%s" % (p.curr_dir, masker_dir) if p.curr_dir else masker_dir
        if not os.path.exists('./%s' % (masker_dir)):
            os.makedirs('./%s' % (masker_dir))
    masker_file = os.path.basename(p.seq_file) + ".out"
    masker_output = "%s/%s" % (masker_dir, masker_file) if masker_dir else "%s/%s" % (p.curr_dir, masker_file) if p.curr_dir else masker_file
    lib_output = re.sub("((\.fa)|(\.fasta))$", ".lib.fa" , p.output)
    library_part = "-lib %s " % (lib_output)
    processor_part = "-pa %d " % (num_processors) if num_processors else ""
    output_part = "-dir %s " % (masker_dir) if masker_dir else "-dir %s " % (p.curr_dir) if p.curr_dir else ""
    cmd = "RepeatMasker %s %s %s %s" % (library_part, processor_part, output_part , p.seq_file)
    print(cmd)
    p2 = pbsJobHandler(batch_file = "repeatmasker", executable = cmd, RHmodules = ["RepeatMasker", "python-3.3.3"]) 
    p2.submit()
    p2.seq_file = p.seq_file
    p2.consensus = lib_output
    p2.curr_dir = curr_dir
    p2.masker_dir = masker_dir
    p2.masker_output = masker_output
    return p2

def performance_stats(p, true_repeats, stats_dir, print_rpts):
    """Given the pbs object used to start a repeatmasker job as well as the original
    repeat file, wait until the job is done and then invoke perform_stats.py on
    the original chromosome file, the original repeat file, the simulated sequence
    file, and the masker output file. Put results in stats_dir (Current dir if unspecified)"""
    p.wait()
    stats_file = re.sub("((\.fa)|(\.fasta))$", ".stats", os.path.basename(p.seq_file))
    stats_out = "%s/%s"%(stats_dir, stats_file) if stats_dir else "%s/%s" % (p.curr_dir, stats_file) if p.curr_dir else stats_file
    print_part = "--print " if print_rpts else "" 
    cmd = "python3.3 perform_stats.py %s %s %s %s %s %s" % (print_part, args.chromosome, true_repeats, p.seq_file, p.masker_output, stats_out)
    print(cmd)
    p2 = pbsJobHandler(batch_file = "stats", executable = cmd)
    p2.submit()
    p2.curr_dir = p.curr_dir
    p2.stats_dir = stats_dir
    p2.stats_output = stats_out
    return p2

def performance_sum(stats_jobs, stats_dir, curr_dir):
    """Given a list of all of the statistics jobs, uses the statistics output files to
    generate a summary file indicative of overall performance. Put results in stats_dir
    (Current dir if unspecified)"""
    tps = 0
    tns = 0
    fps = 0
    fns = 0
    for p in stats_jobs:
        sf = open(p.stats_output, "r")
        tps += int(re.split("\s+", sf.readline())[1])
        fps += int(re.split("\s+", sf.readline())[1])
        tns += int(re.split("\s+", sf.readline())[1])
        fns += int(re.split("\s+", sf.readline())[1])
        sf.close()
    smry_path = "%s/summary.stats" % (stats_dir) if stats_dir else "%s/summary.stats" %(curr_dir) if curr_dir else "summary.stats"
    
    smry = open(smry_path, 'w')
    smry.write("Evaluation completed.\n")
    smry.write("True Positives (TP): \t %d \n" % (tps))
    smry.write("False Positives (FP): \t %d \n" % (fps))
    smry.write("True Negatives (TN): \t %d \n" % (tns))
    smry.write("False Negatives (FN): \t %d \n" % (fns))
    smry.write("\nPerformance of Repeat Classification Tool\n")
    smry.write("Sensitivity (TPR): \t\t %f %%\n" % (tps/(tps + fns)))
    smry.write("Specificity (TNR): \t\t %f %%\n" % (tns/(fps + tns)))
    smry.write("Precision (PPV): \t\t %f %%\n" % (tps/(tps + fps)))
    smry.write("Neg. Pred. Val. (NPV): \t %f %%\n" % (tns/(tns + fns)))
    smry.write("Fall-Out (FPR): \t\t %f %%\n" % (fps/(fps + tns)))
    smry.write("False Disc. Rate (FDR): \t  %f %%\n" % (fps/(tps + fps)))
    smry.close()


if __name__:
    args = parse_params(sys.argv[1:])
    if args.organize:
        prefix_part = "REV_%d_%d_" % (args.f, len(args.seed))
        curr_dir = tempfile.mkdtemp(prefix = prefix_part, dir = ".")
        if not os.path.exists('./%s' % (curr_dir)):
            os.makedirs('./%s' % (curr_dir))
    else:
        curr_dir = None
    if args.subparser_name == "seq_files" and args.seq_files:	
        ### For each listed file (in args.seq_files): invoke RAIDER on the file
        ### and put the resuling pbs object ito the J list.
        J = []
        for file in args.seq_files:
            assert file.endswith(".fa") or file.endswith('.fasta')
            J.append(run_raider(seed = args.seed, f = args.f, m = args.min, input_file = file, output_dir = args.output_dir, curr_dir = curr_dir))
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
        ### Sequentially work through the J2 list and, when next job has finished,
        ### start repeatmasker job running and stick new pbs job onto J3 list
        J3 = []
        for p in J2:
            J3.append(run_repeat_masker(p, args.pa, args.masker_dir))
        for p in J3:
            p.wait()

    elif args.subparser_name == "chrom_sim" and args.chromosome and args.repeat and args.num_sims:
        ### Create n simulated chromosomes, invoke RAIDER on each chromosome
        ### and put the resulting pbs object into the J list.
        J =[]
        for i in range(int(args.num_sims)):
            J.append(simulate_chromosome(chromosome = args.chromosome, repeat = args.repeat, 
                                         rng_seed = args.rng_seed, length = args.length, 
                                         neg_strand = args.negative_strand, fam_file = args.family_file, 
                                         chrom_dir = args.chrom_dir, output_file= args.output, file_index = i, 
                                         curr_dir = curr_dir, k = args.k))  # KARRO: Added k, reformated for readability.  (FYI: the EOL \ not needed if you 
                                                                            # are in the middle of a structure (e.g. a parameter list), where the parser can tell the line is not done.

        ### Sequentially work through the J list and, when next job is finished,
        ### start raider job running and put resulting pbs job object into J2 list
        J2 = []
        for p in J:
            J2.append(run_raider_chrom(p, seed = args.seed, f = args.f, m = args.min,  output_dir = args.output_dir))
        ### Sequentially work through the J2 list and, when next job is finished,
        ### start consensus sequence job running and stick new pbs job onto J3 list
        J3 = []
        for p in J2:
            output = re.sub("((\.fa)|(\.fasta))$", ".consensus%sfa" % ("." if not args.output_ext else "." + output_ext + "."), p.file)
            J3.append(create_raider_consensus(p, output))
        ### Sequentially work through the J3 list and, when next job is finished,
        ### start repeatmasker job running and stick new job onto J4 list
        J4 = []
        for p in J3:
            J4.append(run_repeat_masker(p, args.pa, args.masker_dir))
        ### Sequentially work through the J4 list and, when next job is finished,
        ### start statistics job running and stick new job onto J5 list
        J5 = []
        stats_dir = "%s/%s" % (curr_dir, args.stats_dir) if curr_dir and args.stats_dir else args.stats_dir
        if args.stats_dir and not os.path.exists('./%s' % (stats_dir)):
            os.makedirs('./%s' % (stats_dir))
        for p in J4:
            J5.append(performance_stats(p, args.repeat, stats_dir, args.print_reps))
        ### Wait for all statistics jobs to complete. Submit the list of files to the
        ### performance_sum method to generate the final summary of statistics
        for p in J5:
            p.wait()
        performance_sum(J5, stats_dir, curr_dir)


