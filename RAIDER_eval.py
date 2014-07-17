#!/software/python/3.3.3/bin/python3.3
import sys
import subprocess
import os
import os.path
import argparse
from redhawk import *
import tempfile
import re

#################################################################
# The following global variables are related to debugging only.
show_progress = False
simulate_only = False
#################################################################

def parse_params(args):
    """Parse command line arguments using the argparse library"""
    parser = argparse.ArgumentParser(description = "Evaluate RAIDER against RepeatScout")
    
    # GENERAL ARGUMENTS
    parser2 = parser.add_mutually_exclusive_group()
    parser2.add_argument('--organize', action = "store_true", help = "Create directory for all Raider Eval output", default = False)
    parser2.add_argument('--no', '--named_organize', dest = "named_organize", help = "Organize under a named directory", default = None)

    # TOOL SELECTION
    parser_tools = parser.add_argument_group("tool selection (all on by default)")
    parser_tools.add_argument('-R', '--raider_off', dest = 'run_raider', action = 'store_false', help = 'Turn RAINDER off', default = True)
    parser_tools.add_argument('--RS', '--repscout_off', dest = 'run_repscout', action = 'store_false', help = 'Turn RAINDER off', default = True)
    # Will later add: RepeatModeler, RECON, PILER (other?)


    # I/O ARGUMENTs
    parser_io = parser.add_argument_group("i/o arguments")
    parser_io.add_argument('-r', '--results_dir', dest = "results_dir", help = "Directory containing all results", default = "RAIDER_EVAL")
    parser_io.add_argument('--rd', '--radier_dir', dest = "raider_dir", help = "Subdirectory containing raider results", default = "RAIDER")
    parser_io.add_argument('--rsd', '--rpt_scout_dir', dest = 'rpt_scout_dir', help = "Subdirectory containing rpt scout results", default = "RPT_SCT")


    
    # RAIDER ARGUMENTS
    raider_argument = parser.add_argument_group("RAIDER parameters")
    raider_argument.add_argument('-f', type = int, help = "E.R. occurrence threshold", default = 2)
    raider_argument.add_argument('-m', '--min', type = int, help = "Minimum repeat length. Defaults to pattern length.", default = None)
    raider_argument.add_argument('-s', '--seed', help = "Spaced seed string", default = "1111011110111101111")
    raider_argument.add_argument('-d', '--output_dir', help = "Raider output directory", default = None)
    raider_argument.add_argument('-e', '--output_ext', help = "Output Extension", default = None)
    raider_argument.add_argument('-C', '--cleanup_off', dest = "cleanup", action = "store_false", help = "Turn off file cleanup", default = True)
    

    # REPEAT MASKER ARGUMENTS
    parser.add_argument('--masker_dir', help = "Repeat masker output directory", default = None)
    parser.add_argument('-p', '--pa', type = int, help = "Number of processors will be using")
    
    # STATISTICS ARGUMENT
    stats_group = parser.add_argument_group(title = "Statistics argument")
    stats_group.add_argument('--stats_dir', help = "Statistics output directory", default = None)
    stats_group.add_argument('--print_reps', action = "store_true", help = "Print out repeats in statistics file", default = False)

    # DEBUGGING ARGUMENTS
    debug_group = parser.add_argument_group(title = "debugging")
    debug_group.add_argument('--sp', '--show_progress', dest = 'show_progress', action = 'store_true', help = "Print reports on program progress to stdout", default = False)
    debug_group.add_argument('--so', '--simulate_only', dest = 'simulate_only', action = 'store_true', help = "Quit after creating simulated file", default = False)

    

    subparsers = parser.add_subparsers(dest="subparser_name")
    
    # SEQUENCE FILE OPTION ARGUMENTS
    parser_seqs = subparsers.add_parser("seq_files")
    parser_seqs.add_argument('seq_files', nargs = "+", help = "Files containing genomic sequence (for running on real data)")
    
    # CHROMOSOME SIMULATION OPTION ARGUMENTS
    parser_chrom = subparsers.add_parser("chrom_sim")
    parser_chrom.add_argument('chromosome', help = "Template chromosome file")
    parser_chrom.add_argument('repeat', help = "Repeat file")
    parser_chrom.add_argument('num_sims', type = int, help ="Number of simulations")
    parser_chrom.add_argument('--sim_dir', help = "Directory containing the resulting simulated chromosome", default = "SIM")
    parser_chrom.add_argument('--rng_seed', type = int, help = "RNG seed", default = None) 
    parser_chrom.add_argument('-l', '--length', type = int, help = "Simulated sequence length", default = None)
    parser_chrom.add_argument('-k', type = int, help = "Order of markov chain", default = 5)  # KARRO: Added this
    parser_chrom.add_argument('-n', '--negative_strand', action = "store_true", help = "Use repeats on negative string", default = False)
    parser_chrom.add_argument('-f', '--family_file', help = "List of repeat families to use", default = None)
    parser_chrom.add_argument('-o', '--output', help = "Output file (Default: replace chromosome file \".fa\" with \".sim.fa\")")
    
    
    arg_return =  parser.parse_args(args)
    
    #### The following is to set the global debugging variables 
    global show_progress
    show_progress = arg_return.show_progress

    global simulate_only
    simulate_only = arg_return.simulate_only
    ###

    return arg_return

def simulate_chromosome(chromosome, repeat, rng_seed, length, neg_strand, fam_file, sim_dir, output_file, file_index, k):
    """Given chromosome file and repeat file and rng_seed, runs chromosome 
    simulator and then passes raider params (including path to new simulated chromosome 
    file) into run_raider"""

    print("sim_dir: ", sim_dir)
    if not os.path.exists(sim_dir):
        os.makedirs(sim_dir)

    # Output file is either specified or replace .fa with .sim.#.fa
    length_arg = "-l %d" % (length) if length else ""
    k_arg = "-k %d" % (k)
    seed_arg = "-s %d" % (rng_seed) if rng_seed else ""
    neg_arg = "-n" if neg_strand else ""
    fam_arg = "-f %s" % (fam_file) if fam_file else ""
    seq_arg = chromosome
    repeat_arg = repeat

    output_file = (output_file if output_file else re.sub(".fa$", ".sim.%d.fa" % (file_index), os.path.basename(chromosome)))
    output_path = "%s/%s" % (sim_dir, output_file)

    cmd = "python3.3 chromsome_simulator.py -m {length} {k} {seed} {neg} {fam} {seq} {repeat} {output}".format(length=length_arg, k=k_arg, seed=seed_arg, neg=neg_arg, fam=fam_arg, seq=seq_arg, repeat=repeat_arg, output=output_path)
    if show_progress:
        print("Creating simulation:\n", cmd)

    #p = pbsJobHandler(batch_file = "%s.batch" % (output_file), executable = cmd)
    p = pbsJobHandler(batch_file = "%s.batch" % (output_file), executable = cmd)
    p.submit()
    p.sim_output = output_path
    p.index = file_index
    return p

# def run_raider_chrom(p, seed, f, m, output_dir):
#    """Given the pbs object used to start simulated chromosome job as well as 
#    raider parameters, wait until the job is done and then call run_raider 
#    on the newly generated chromosome output file"""
#    p.wait()
#    return run_raider(seed, f, m, p.chrom_output, output_dir, p.curr_dir)

def run_raider(seed, f, m, input_file, output_dir, curr_dir):
    """Given raider parameters and an input file, run RAIDER and put the output into
    the directory specified in output_dir (creating a random name is none is
    specified."""
    if output_dir:
        output_dir = "%s/%s" % (curr_dir, output_dir) if curr_dir else output_dir
    else:
        dir_part = curr_dir if curr_dir else "."
        output_dir = output_dir if output_dir else tempfile.mkdtemp(prefix = "raider_%s_"%(os.path.basename(input_file)), dir = dir_part)
    if not os.path.exists('./%s' % (output_dir)):
        os.makedirs('./%s' % (output_dir))
    min_arg = "-m %d" % (m) if m else ""
    cmd = "./raider -q -c %d %s %s %s %s" % (f, min_arg, seed, input_file, output_dir)
    if show_progress:
        print("Launching raider: ", cmd)
    p = pbsJobHandler(batch_file = "raider", executable = cmd)
    p.submit()
    p.seq_file = input_file
    p.curr_dir = curr_dir
    p.raider_output = output_dir
    return p

def create_raider_consensus(p, output):
    """Given the pbs object (from redhawk.py) used to start a RAIDER job, this
    waits until the job is done, then invokes consensus_seq.py on the ouput and
    writes the results to the output directory."""
    p.wait();
    cmd = "python3.3 consensus_seq.py -s %s -e %s/elements %s" % (p.seq_file, p.raider_output, output)
    #cmd = "python3.3 consensus_seq.py -s %s -e %s/elements %s" % (p.file, p.raider_output, output)
    if show_progress:
        print("Creating consensus sequence: ", cmd)
    p2 = pbsJobHandler(batch_file = "%s.batch" % (os.path.basename(output)), executable = cmd)
    p2.submit()
    p2.curr_dir = p.curr_dir
    p2.seq_file = p.seq_file
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
    if show_progress:
        print("Launch repeatmasker: ", cmd)
    p2 = pbsJobHandler(batch_file = "repeatmasker", executable = cmd, RHmodules = ["RepeatMasker", "python-3.3.3"]) 
    p2.submit()
    p2.seq_file = p.seq_file
    p2.consensus = lib_output
    p2.curr_dir = curr_dir
    p2.masker_dir = masker_dir
    p2.masker_output = masker_output
    return p2

def run_scout_chrom(p, f, m):
    p.wait()
    return run_scout(p.chrom_output, f, m, p.curr_dir)

def run_scout(input_file, f, m, curr_dir):
    output_file = re.sub(".fa$", ".scout.fa",os.path.basename(input_file))
    
    output_path = "%s/%s" %(curr_dir, output_file)
    cmd = "./RunRepeatScout.sh {sequence} {length} {output} {min}".format(sequence=input_file, length=m, output=output_path, min=f)
    print(cmd)
    p = pbsJobHandler(batch_file = "%s.batch" % output_file, executable = cmd)
    p.submit()
    p.seq_file = input_file
    p.scout_output = output_path
    p.curr_dir = curr_dir
    return p

def run_rm_scout(p, num_processors, masker_dir):
    p.wait()
    if masker_dir:
        masker_dir = "%s/%s" % (p.curr_dir, masker_dir) if p.curr_dir else masker_dir
        if not os.path.exists('./%s' % (masker_dir)):
            os.makedirs('./%s' % (masker_dir))
    masker_file = os.path.basename(p.seq_file) + ".out"
    pa_part = num_processors if num_processors else 4
    masker_output="%s/%s" % (masker_dir, masker_file) if masker_dir else "%s/%s" % (p.curr_dir, masker_file) if p.curr_dir else masker_file
    cmd = "./MaskRepeatScout.sh {sequence} {output} {pas} {dir}".format(sequence=p.seq_file, output=p.scout_output, pas=pa_part, dir=masker_dir)
    print(cmd)
    p2 = pbsJobHandler(batch_file = "repeatmaskerscout", executable = cmd, RHmodules = ["RepeatMasker", "python-3.3.3"]) 
    p2.submit()
    p2.seq_file = p.seq_file
    p2.consensus = p.scout_output
    p2.curr_dir = p.curr_dir
    p2.masker_dir = masker_dir
    p2.masker_output = masker_output
    return p2


def performance_stats(p, true_repeats, stats_dir, print_rpts, test):
    """Given the pbs object used to start a repeatmasker job as well as the original
    repeat file, wait until the job is done and then invoke perform_stats.py on
    the original chromosome file, the original repeat file, the simulated sequence
    file, and the masker output file. Put results in stats_dir (Current dir if unspecified)"""
    p.wait()
    stats_file = re.sub("((\.fa)|(\.fasta))$", ".%s.stats"%(test) , os.path.basename(p.seq_file))
    stats_out = "%s/%s"%(stats_dir, stats_file) if stats_dir else "%s/%s" % (p.curr_dir, stats_file) if p.curr_dir else stats_file
    print_part = "--print " if print_rpts else "" 
    cmd = "python3.3 perform_stats.py %s %s %s %s %s %s" % (print_part, args.chromosome, true_repeats, p.seq_file, p.masker_output, stats_out)
    if show_progress:
        print("Launching analysis: ", cmd)
    p2 = pbsJobHandler(batch_file = "stats", executable = cmd)
    p2.submit()
    p2.curr_dir = p.curr_dir
    p2.stats_dir = stats_dir
    p2.stats_output = stats_out
    return p2

def performance_sum(stats_jobs, stats_dir, curr_dir, test):
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
    stats_file = "summary.%s.stats" % (test)
    smry_path = "%s/%s" % (stats_dir, stats_file) if stats_dir else "%s/%s" %(curr_dir, stats_file) if curr_dir else stats_file
    
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


if __name__ == "__main__":
    args = parse_params(sys.argv[1:])

    # ### Setup directory organization
    # if args.organize or args.named_organize:
    #     if args.organize:
    #         prefix_part = "REV_%d_%d_" % (args.f, len(args.seed))
    #         curr_dir = tempfile.mkdtemp(prefix = prefix_part, dir = ".")
    #     else:
    #         curr_dir = args.named_organize
    #         if not os.path.exists(curr_dir):
    #             os.makedirs(curr_dir)

    #     if not os.path.exists('./%s' % (curr_dir)):
    #         os.makedirs('./%s' % (curr_dir))
    # else:
    #     curr_dir = None
    if not os.path.exists(args.results_dir):
        os.makedirs(args.results_dir)



    ### Generate simulated file(s) and run to completion
    if args.subparser_name == "chrom_sim":
        sim_output = args.results_dir + "/" + args.sim_dir

        # Launch the jobs
        f = lambda i: simulate_chromosome(chromosome = args.chromosome, repeat = args.repeat, 
                                          rng_seed = args.rng_seed, length = args.length, 
                                          neg_strand = args.negative_strand, fam_file = args.family_file, 
                                          sim_dir = args.results_dir + "/" + args.sim_dir, output_file = args.output, file_index = i, 
                                          k = args.k)
        J = [f(i) for i in range(args.num_sims)]

        # Run jobs to completion
        [j.wait() for j in J]    # Let all the simulations finish

        # Quit (if done)
        if simulate_only:
            [j.erase_files() for j in J]
            sys.exit(0)

        # Get the list of simulated file names
        file_list = [j.sim_output for j in J]

    else:
        # Get the list of file names
        file_list = seq_files

    ### Run RAIDER
    #RAIDER_JOBS = [run_raider(seed = args.seed, f = args.f, m = args.m, input_file = file, output_dir = args.output_dir, curr_dir = curr_dir) for files in file_list]
    #CONSENSUS_JOBS = [create_raider_consensus(p, re.sub("((\.fa)|(\.fasta))$", ".consensus%sfa" % ("." if not args.output_ext else "." + output_ext + "."), p.seq_file)) for p in RAIDER_JOBS]
    #REPEATMASKER_JOBS = 


    # ####################################################################
    # if args.subparser_name == "seq_files" and args.seq_files:

    #     ### For each listed file (in args.seq_files): invoke RAIDER on the file
    #     ### and put the resuling pbs object ito the J list.
    #     J = []
    #     for file in args.seq_files:
    #         assert file.endswith(".fa") or file.endswith('.fasta')
    #         J.append(run_raider(seed = args.seed, f = args.f, m = args.min, input_file = file, output_dir = args.output_dir, curr_dir = curr_dir))
    #     ### Sequentially work the the J list and, when the next job has finished,
    #     ### start the conensus sequence job running and stick the new pbs job nto
    #     ### the J2 job. 
    #     J2 = []
    #     for p in J:
    #         # The below line is creating the output name in such a way as to avoid conflict.
    #         # It will change X.fa to X.conensus.fa (if args.output_ext == None), stick that output_ext
    #         # string in there if it does exist.
    #         output = re.sub("((\.fa)|(\.fasta))$", ".consensus%sfa" % ("." if not args.output_ext else "." + output_ext + "."), p.file)
    #         J2.append(create_raider_consensus(p, output))

    #     ### Sequentially work through the J2 list and, when next job has finished,
    #     ### start repeatmasker job running and stick new pbs job onto J3 list
    #     J3 = []
    #     for p in J2:
    #         J3.append(run_repeat_masker(p, args.pa, args.masker_dir))
    #     for p in J3:
    #         p.wait()

    # elif args.subparser_name == "chrom_sim" and args.chromosome and args.repeat and args.num_sims:
    #     ### Create n simulated chromosomes, invoke RAIDER on each chromosome
    #     ### and put the resulting pbs object into the J list.
    #     J =[]
    #     for i in range(int(args.num_sims)):
    #         J.append(simulate_chromosome(chromosome = args.chromosome, repeat = args.repeat, 
    #                                      rng_seed = args.rng_seed, length = args.length, 
    #                                      neg_strand = args.negative_strand, fam_file = args.family_file, 
    #                                      chrom_dir = args.chrom_dir, output_file= args.output, file_index = i, 
    #                                      curr_dir = curr_dir, k = args.k)) 

    #     if simulate_only:
    #         [r.wait() for r in J]
    #         [r.erase_files() for r in J]
    #         sys.exit(0)


    #     ### Sequentially work through the J list and, when next job is finished,
    #     ### start raider job running and put resulting pbs job object into J2 list
    #     J2 = []
    #     for p in J:
    #         J2.append(run_raider_chrom(p, seed = args.seed, f = args.f, m = args.min,  output_dir = args.output_dir))

    #     ### Sequentially work through the J2 list and, when next job is finished,
    #     ### start consensus sequence job running and stick new pbs job onto J3 list

    #     J3 = []
    #     for p in J2:
    #         output = re.sub("((\.fa)|(\.fasta))$", ".consensus%sfa" % ("." if not args.output_ext else "." + output_ext + "."), p.seq_file)
    #         J3.append(create_raider_consensus(p, output))

    #     rm_raider_dir = "RM_raider" #, dir = (curr_dir if curr_dir else "."))
    #     ### Sequentially work through the J3 list and, when next job is finished,
    #     ### start repeatmasker job running and stick new job onto J4 list
    #     J4 = []
    #     for p in J3:
    #         J4.append(run_repeat_masker(p, args.pa, rm_raider_dir))
    #     ### Sequentially work through the J4 list and, when next job is finished,
    #     ### start statistics job running and stick new job onto J5 list
    #     J5 = []
    #     stats_dir = "%s/%s" % (curr_dir, args.stats_dir) if curr_dir and args.stats_dir else args.stats_dir
    #     if args.stats_dir and not os.path.exists('./%s' % (stats_dir)):
    #         os.makedirs('./%s' % (stats_dir))
    #     for p in J4:
    #         J5.append(performance_stats(p, args.repeat, stats_dir, args.print_reps, "raider"))
    #     ### Wait for all statistics jobs to complete. Submit the list of files to the
    #     ### performance_sum method to generate the final summary of statistics
    #     J6 = []
    #     for p in J5:
    #         p.wait()
    #     performance_sum(J5, stats_dir, curr_dir, "raider")   # J6.append(p.stats_output)
        
        
    #     J7 = []
    #     for p in J:
    #         J7.append(run_scout_chrom(p, f = args.f, m = len(args.seed)))

    #     rm_scout_dir = "RM_scout" #, dir = (curr_dir if curr_dir else "."))
    #     J8 = []
    #     for p in J7:
    #         J8.append(run_rm_scout(p, args.pa, rm_scout_dir))
    #     J9 = []
    #     for p in J8:
    #         J9.append(performance_stats(p, args.repeat, stats_dir, args.print_reps, "scout"))
    #     J10 = []
    #     for p in J9:
    #         p.wait()
    #     performance_sum(J9, stats_dir, curr_dir, "scout")
    #     #summary_stats(stats, stats_dir, curr_dir)
    #     #summary_job.wait()
    #     #performance_sum(J5, stats_dir, curr_dir)


