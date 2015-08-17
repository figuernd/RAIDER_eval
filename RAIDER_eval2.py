import sys
import os
import argparse
import locations
import subprocess
import pickle
from redhawk import *

# Second version of the script for testing RAIDER.
# Generates runs only; will use a seperate script for testing.



#######################
# Globals: So we don't have to pass everything through parameters.
args = None     # Command line arguments object
Locations = locations.Locations
progress_fp = None
job_log_dir = None
seed_map = {}         # Maps seed to (index,short rep.) pairs
seed_list = None      # Sorted list of seeds

#######################
# Defaults 
walltime_default = "4:00:00"
delay_default = 120           # Number of seconds to sleep when cycling on a redhawk wait

#######################
# Useful utilitiy functions
def file_base(file):
    """Extract the name of a file froma directory"""
    return os.path.basename(file)

def file_dir(file):
    """Extract the directory a file is contained in"""
    return file.rstrip(file_base(file)).rstrip("/")

def convert_seed(seed):
    """Convert an abriviated seed to a full seed (e.g. "1{2}0{3}1{2}" => "1100011" """
    i = 0
    while (i < len(seed)-1):
        if seed[i+1] == '^':
            j = i+2
            assert seed[j] == "{"
            k = j+1
            while seed[k] != '}':
                k += 1
            n = int(seed[j+1:k])
            seed = seed[:i] + seed[i]*n + seed[k+1:]
        i += 1
    return seed

def seed_weight(seed):
    return sum([x == 1 for x in seed])

def cmp_seed(s1, s2):
    d = len(s1) - len(s2)
    if (d != 0):
        return d < 0
    d = weight(s1) - weight(s2)
    if f != 0:
        return d < 0
    return s1 < s2

#######################
# Command line parsing
def parse_params():
    """Parse command line arguments, set them equal to the global args, and set certain global variables"""
    parser = argparse.ArgumentParser(description = "Evaluate RAIDER against RepeatScout")

    # TOOL SELECTION
    parser_tools = parser.add_argument_group("tool selection (all on by default)")
    parser_tools.add_argument('-R', '--raider_on', dest = 'run_raider', action = 'store_true', help = 'Turn RAIDER on', default = False)
    parser_tools.add_argument('--R2', '--raider2_on', dest = 'run_raider2', action = 'store_true', help = 'Turn RAIDERV2 on', default = False)
    parser_tools.add_argument('--AR', '--araider_on', dest = 'run_araider', action = 'store_true', help = 'Turn ARAIDER on', default = False)
    parser_tools.add_argument('--RS', '--repscout_on', dest = 'run_repscout', action = 'store_true', help = 'Turn RAIDER on', default = False)
    parser_tools.add_argument('-B', '--bigfoot_on', dest = 'run_bigfoot', action = 'store_true', help = 'Turn BIGFOOT on', default = False)
    parser_tools.add_argument('-P', '--piler_on', dest = 'run_piler', action = 'store_true', help = 'Turn PILER on', default = False)
    parser_tools.add_argument('-A', '--all_tools', dest = 'all_tools', action = 'store_true', help = 'Turn all tools on (overide all other tool arguments)', default = False)
    parser_tools.add_argument('--A2', '--all_tools2', dest = 'all_tools2', action = 'store_true', help = 'Turn all tools on except araider (overide all other tool arguments)', default = False)
    parser_tools.add_argument('--tl', '--time_limit', dest = 'time_limit', help = 'Redhawk time limit (max: 400:00:00 default: 4:00:00)', default = walltime_default)

    # I/O ARGUMENTs
    parser_io = parser.add_argument_group("i/o arguments")
    parser_io.add_argument('-r', '--results_dir', dest = "results_dir", help = "Directory containing all results", default = "EVAL")
    parser_io.add_argument('--nuke', dest ='nuke', action = "store_true", help = "Nuke the results directory", default = False)

    # RAIDER ARGUMENTS
    raider_argument = parser.add_argument_group("RAIDER parameters")
    raider_argument.add_argument('-f', type = int, help = "E.R. occurrence threshold", default = 5)
    raider_argument.add_argument('--pre', '--pre_scan', action = 'store_true', help = "Use pre-scan version of raider", default = False)
    raider_argument.add_argument('--mem', action = 'store_true', help = "Use large memory-nodes", default = False);
    raider_argument.add_argument("--mn", '--max_nodes', dest = "max_nodes", action="store_true", help="Reserve all nodes of a processor for each tool (disabled by default).", default=False)
    seed_group = raider_argument.add_mutually_exclusive_group(required = False)     
    seed_group.add_argument('-s', '--seed', dest = "seed", help = "Spaced seed string", default = "111111111111111111111111111111")    
    seed_group.add_argument('--sf', '--seed_file', dest = 'seed_file', help = 'File containing raider seeds', default = None)

    # RAIDER2 ARGUMENTS
    raider2_argument = parser.add_argument_group("RAIDER2 parameters")
    raider2_argument.add_argument('--age', type = int, help="Use older version of raider2", default=1) 
    raider2_argument.add_argument('--aa', '--all_ages', dest="all_ages", action="store_true", help="Run all ages of raider2", default=False) # type = int, help="Use older version of raider", default=0)
    #raider2_argument.add_argument('--multi', '--multi_seed', dest="multi_seed", action="store_true", help="Run all seeds in seed file concurrently",default=False)
    raider2_argument.add_argument('--na', '--no_family_array', dest="family_array", action="store_false", help="Disable family array in Raider2", default=True)
    raider2_argument.add_argument('--ex', '--excise', dest="excising", action="store_true", help="Enable excising in RAIDER2", default=False)
    raider2_argument.add_argument('--no', '--no_overlaps', dest="overlaps", action="store_false", help="Do not require overlaps in RAIDER2", default=True)
    raider2_argument.add_argument('--tu', '--tie_up', dest="tieup", action="store_true", help="Enable alternative tie ups", default=False)
    raider2_argument.add_argument('--ps', '--prosplit', dest="prosplit", action="store_true", help="Enable proactive splitting(disabled by default).", default=False)
    raider2_argument.add_argument("--pf", '--prevfam', dest="prevfam", action="store_true", help="Enable pointers to prev family (disabled by default).", default=False)

    # REPSCOUT ARGUMENTS
    repscout_argument = parser.add_argument_group("REPSCOUT parameters")
    repscout_argument.add_argument('--repscout_min', type = int, help = "Minimum repeat length for repscout.", default = 10)
    repscout_argument.add_argument('--rs_min_freq', type = int, help = "Minimum repeat length for repscout.", default = 3)

    # REPEAT MASKER ARGUMENTS
    repeatmasker_arguments = parser.add_argument_group("RepeatMasker parameters")
    repeatmasker_arguments.add_argument('--RM', '--suppress_rm', dest = 'run_rm', action="store_false", help = "suppress RepeatMasker run", default = True)
    repeatmasker_arguments.add_argument('--masker_dir', help = "Repeat masker output directory", default = None)
    repeatmasker_arguments.add_argument('-p', '--pa', type = int, help = "Number of processors will be using", default = 1)    
    repeatmasker_arguments.add_argument('--rwt', '--rm_walltime', dest = "rm_walltime", help = "Wall time limit for repeat masker", default = "10:00:00")

    # DEBUGGING ARGUMENTS
    debug_group = parser.add_argument_group(title = "debugging")
    debug_group.add_argument('--sp', '--show_progress', dest = 'show_progress', action = 'store_true', help = "Print reports on program progress to stderr", default = False)
    debug_group.add_argument('--dry', '--dry_run', dest = 'dry_run', action = 'store_true', help = "Dry run -- don't actually launch jobs", default = False)

    # Positional arguments
    parser.add_argument('data_files', nargs = '+', help = "Data files to process")

    global args;
    args = parser.parse_args(args)


def print_progress(s):
    """Print s to progress_fp.  Also print it to stdout if args.sp is true."""
    progress_fp.write(s + "\n");
    progress_fp.flush()
    if args.show_progress:
        sys.stdout.write(s + "\n")
        sys.stdout.flush()

def launch_job(cmd, title, base_dir, walltime = walltime_default, ppn = 1, bigmem = False, depend = None, modules = None):
    """Launch a redhawk jobs.
    * cmd: Command to be launched
    * title: Used as the bases for creating a job name and redhawk-related files
    * base_dir: Directory to place all redhawk-related files
    * walltime: walltime limit
    * ppn: Number of processors requested
    * bigmem: If true, run only on a large-memory node
    * depend: A list of job objects which must terminate before first
    """
    print_progress(cmd + "\n");
    
    batch_file = base_dir + "/" + title + ".job"
    job_name = title;
    stdout_file = base_dir + "/" + title + ".stdout"
    stderr_file = base_dir + "/" + title + ".stderr"

    p = pbsJobHandler(batch_file = batch_file, executable = cmd, job_name = job_name,
                      stdout_file = stdout_file, stderr_file = stderr_file,
                      walltime = walltime, depends = depend,
                      mem = Locations['high_mem_arch'] if args.mem else False,
                      RHmodules = modules)

    if (not args.dry_run):
        p.submit(preserve = True, delay = default_delay)

    with open(job_log_dir + "/" + title, "w") as fp:
        print("\n".join([str(x) + "\t" + str(getattr(p,x)) + "\t" + str(type(getattr(p,x))) for x in dir(p)]))
        pickle.dump([p], fp)
        fp.flush()

    return p




def setup():
    global args

    global job_log_dir;
    job_log_dir = args.results_dir + "/job_log/"
    
    global debug_file
    debug_file = args.results_dir + "/debug.txt"

    if args.nuke and os.path.exists(args.results_dir):
        subprocess.call("rm -r %s" % args.results_dir, shell = True)
 
    if not os.path.exists(args.results_dir):
        os.makedirs(args.results_dir)
    if not os.path.exists(job_log_dir):
        os.makedirs(job_log_dir)

    global progress_fp
    progress_fp = open(debug_file, "w")

    global seed_map
    if args.seed_file:
        seed_map = {convert_seed(seed):(i,seed) for i,line in enumerate(open(args.seed_file)) for seed in [line.strip()]}
    else:
        seed_map = {convert_seed(args.seed):(0,args.seed)}

    global seed_list
    seed_list = sorted(seed_map.values(), key = lambda(x): x[0])
    with open(args.results_dir + "/seed_file.txt", "w") as fp:
        fp.write("\n".join([str(x[0]) + "\t" + x[1] for x in seed_list]))


        
raider_cmd = "{raider} -q -c {f} -s {seed} {input_file} {output_dir}"
consensus_cmd = "{python} consensus_seq.py -s {data_file} -e {elements_file} {consensus_txt} {consensus_fa}"
repeat_masker_cmd = "{RepeatMasker} -nolow -lib {library_file} -pa {pa} -dir {output_dir} {seq_file}"

def create_raider2_pipeline(input_file, seed, f):
    ##########################
    # SETUP
    input_base = file_base(input_file).strip(".fa")
    raider2_dir = args.results_dir + "/RAIDERV2";    # Directory of RAIDER results
    if not os.path.exists(raider2_dir):
        os.makedirs(raider2_dir)

    title = "{file}.{seed_index}.{f}".format(file=input_base, seed_index=seed_map[seed][0], f=f)

    seed_index = seed_map[seed][0];                     
    consensus_name = input_base + ".s" + str(seed_index) + ".f" + str(f)

    elements_dir = raider2_dir + "/" + consensus_name.upper() + ".RM"
    if not os.path.exists(elements_dir):
        os.makedirs(elements_dir)

    consensus_txt = raider2_dir + "/" + consensus_name + ".consensus.txt"
    consensus_fa = raider2_dir + "/" + consensus_name + ".consensus.fa"

    rm_dir = raider2_dir + "/" + consensus_name
    if not os.path.exists(rm_dir):
        os.makedirs(rm_dir)


    ##########################
    # Step 1: Run phRAIDER
    cmd1 = raider_cmd.format(raider=Locations['raider2'], f=f, seed=seed,
                             input_file=input_file, output_dir=elements_dir)
    title1 = "phRA." + title;
    p1 = launch_job(cmd=cmd1, title=title1, base_dir=elements_dir, ppn = Locations['high_mem_arch'] if args.max_nodes else 1, bigmem = args.mem)

    

    # Step 2: Generate consensus sequence
    cmd2 = consensus_cmd.format(python=Locations['python'], data_file=input_file,
                                elements_file=elements_dir + "/elements",
                                consensus_txt=consensus_txt,
                                consensus_fa=consensus_fa)
    title2 = "con." + title
    p2 = launch_job(cmd=cmd2, title=title2, base_dir=raider2_dir, depends = [p1])

    # Step 3: Apply repeat masker consensus output
    if (args.run_rm):
        cmd3 = repeat_masker_cmd.format(RepeatMasker = Locations['RepeatMasker'],
                                        library_file = consensus_fa, pa = args.pa,
                                        output_dir = rm_dir,
                                        seq_file = input_file)
        title3 = "phRA-RM" + title
        p3 = launch_job(cmd=cmd3, title=title3, base_dir=rm_dir, walltime = args.rm_walltime, ppn = args.pa, bibmem = False, modules = Locations['rm_modules'])
    




if __name__ == "__main__":
    parse_params()
    setup()
    A, B = create_raider2_pipeline(args.data_files[0], args.seed, args.f)
    print("\n".join(A))
