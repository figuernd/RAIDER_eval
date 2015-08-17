import sys
import os
import redhawk
import argparse
import locations

# Second version of the script for testing RAIDER.
# Generates runs only; will use a seperate script for testing.



#######################
# Globals: So we don't have to pass everything through parameters.
args = None     # Command line arguments object
default_time_limit = "4:00:00"
Locations = locations.Locations
progress_fp = None
seed_map = {}         # Maps seed to (index,short rep.) pairs
seed_list = None      # Sorted list of seeds

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
    parser_tools.add_argument('--tl', '--time_limit', dest = 'time_limit', help = 'Redhawk time limit (max: 400:00:00 default: 4:00:00)', default = default_time_limit)
    parser_tools.add_argument("--mn", '--max_nodes', dest = "max_nodes", action="store_true", help="Reserve all nodes of a processor for each tool (disabled by default).", default=False)

    # I/O ARGUMENTs
    parser_io = parser.add_argument_group("i/o arguments")
    parser_io.add_argument('-r', '--results_dir', dest = "results_dir", help = "Directory containing all results", default = "EVAL")
    parser_io.add_argument('--nuke', dest ='nuke', action = "store_true", help = "Nuke the results directory", default = False)

    # RAIDER ARGUMENTS
    raider_argument = parser.add_argument_group("RAIDER parameters")
    raider_argument.add_argument('-f', type = int, help = "E.R. occurrence threshold", default = 5)
    raider_argument.add_argument('-d', '--output_dir', help = "Raider output directory", default = None)
    raider_argument.add_argument('-e', '--output_ext', help = "Output Extension", default = None)
    raider_argument.add_argument('-C', '--cleanup_off', dest = "cleanup", action = "store_false", help = "Turn off file cleanup", default = True)
    raider_argument.add_argument('--raider_min', '--raider_min', type = int, help = "Minimum repeat length. Defaults to pattern length.", default = None)
    raider_argument.add_argument('--pre', '--pre_scan', action = 'store_true', help = "Use pre-scan version of raider", default = False)
    raider_argument.add_argument('--mem', action = 'store_true', help = "Use large memory-nodes", default = False);
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

    # DEBUGGING ARGUMENTS
    debug_group = parser.add_argument_group(title = "debugging")
    debug_group.add_argument('--sp', '--show_progress', dest = 'show_progress', action = 'store_true', help = "Print reports on program progress to stderr", default = False)

    # Positional arguments
    parser.add_argument('data_files', nargs = '+', help = "Data files to process")

    global args;
    args = parser.parse_args(args)

def setup():
    global args
    if args.nuke and os.path.exists(args.results_dir):
        subprocess.call("rm -r %s" % args.results_dir, shell = True)
 
    os.makedirs(args.results_dir)

    global debug_fp;
    progress_fp = open(args.results_dir + "/debug.txt", "w")

    global seed_map
    if args.seed_file:
        seed_map = {convert_seed(seed):(i,seed) for i,line in enumerate(open(args.seed_file)) for seed in [line.strip()]}
    else:
        seed_map = {consert_seed(args.seed):(0,args.seed)}


        
raider_cmd = "{raider} -q -c {f} -s {seed} {input_file} {output_dir}"
consesus_cmd = "{python} consensus_seq.py -s {data_file} -e {elements_file} {consensus_txt} {consensus_fa}"
repeat_masker_cmd = "{RepeatMasker} -nolow -lib {library_file} -pa {pa} -dir {output_dir} {seq_file}"

def create_raider2_pipeline(input_file, seed, f):
    raider2_dir = args.results_dir + "/RAIDERV2/";    # Directory of RAIDER results
    if not os.path.exists(results_dir):
        os.makedirs(outputdir)
        
    seed_index = map_seed[seed];                     
    input_base = file_base(input_file).rstirp(".fa")  # Name of simulation file
    consensus_name = input_base + ".s" + str(seed_num) + ".f" + str(f)
    elements_dir = raider2_dir + "/" + consensus_name.upper() + ".RM"
    rm_dir = raider2_dir + "/" + consensus_name
    consensus_txt = current_dir + "/" + consensus_name + ".txt"
    consensus_fa = current_dir + "/" + consensus_name + ".fa"

    # Step 1: Run phRAIDER
    cmd1 = raider_cmd.format(raider=Locations['raider2'], f=f, seed=seed,
                             input_file=input_file, output_dir=elements_dir)
    title1 = "phRA.{file}.{seed_num}.{f}".format(file=input_base, seed_num=seed_num)

    # Step 2: Generate consensus sequence
    cmd2 = consensus_cmd.format(python=Locations['python'], data_file=input_file,
                                elements_file=elements_dir + "/" + elements,
                                consensus_txt=consensus_txt,
                                consensys_fa=consensus_fa)
    title2 = "phRA.{file}.{seed_num}.{f}"

    # Step 3: Apply repeat masker consensus output
    if (args.run_rm):
        cmd3 = repeat_masker_cmd.format(RepeatMasker = Locations['RepeatMasker'],
                                        lib = consensus_fa, pa = args.ps,
                                        dir = rm_dir,
                                        seq_file = input_file)
        title3 = "phRA-RM.{file}.{feed_num}.{f}"
    else:
        cmd3 = ""
        title3 = ""

    return [cmd1, cmd2, cmd3], [title1, title2, title3]


if __name__ == "__main__":
    parse_params()
    io_setup()
    A, B = create_raider2_pipeline(input, seed, f)
    print("\n".join(a))
