import sys
import perform_stats
import glob
import re
import redhawk
import parse_pra_output
import os
import subprocess
import argparse
tool_prefix = {'phRAIDER':"phRA", "RepeatScout":"RS"}

args = None
seed_map = None
data_map = None
f_list = set()
stats_map = {}


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

stats = ["tool", "seed_num", "l", "w", "w/l", "org", "f", "tp", "fp", "fn", "tn", "tpr", "tnr", "ppv", "npv", "fpr", "fdr", "ToolCpuTime", "ToolWallTime", "ToolMem", "ToolVMem", "RMCpuTime", "RMWallTime", "RMMem", "RMVMem", "ConCoverage", "QuCoverage"]   
print_stats = ["org", "tool", "seed_num", "l", "w", "w/l", "f", "tpr", "tnr", "ToolCpuTime", "ToolMem", "ToolVMem", "ConCoverage", "QuCoverage"]
def create_stats_hash(tool, org, seed_num, f):
    H = {x:None for x in stats}
    H["tool"] = tool
    H["seed_num"] = -1 if seed_num is None else int(seed_num)
    H["org"] = org
    H["f"] = int(f)
    return H


def print_stats_hash(H, all_stats = False, fp = sys.stdout):
    for x in stats:
        if all_stats or not (H[x] is None):
            fp.write(x + ": " + str(H[x]) + "\n")


def load_seeds(DIR):
    global seed_map
    seed_map = {}
    for line in open(DIR + "/seed_file.txt"):
        index, seed = re.split("\s+", line.strip())
        if args.seeds is None or int(index) in args.seeds:
            seed_map[int(index)] = seed


def load_map(DIR):
    global data_map
    data_map = {}
    for line in open(DIR + "/data_files.txt"):
        short_name, long_name = re.split("\s+", line.strip())
        if args.data is None or any([re.search(x, short_name) for x in args.data]):
            data_map[short_name] = long_name

def load_f_list(DIR):
    global f_list
    for line in open(DIR + "/f.txt"):
        for f in re.split("\s+", line.strip()):
            if args.f is None or int(f) in args.f:
                f_list.add(int(f))

        
def collectRaider(DIR, tool):
    # FIRST: Get any applicable RM output stats
    for org in data_map.keys():
        for seed_num in seed_map.keys():
            for f in f_list:
                RAIDER_job_file = DIR + "/../job_log/{prefix}.{org}.s{seed_num}.f{f}".format(prefix = tool_prefix[tool], org=org, seed_num=seed_num, f=f)        
                RM_job_file = DIR + "/../job_log/rm.{prefix}.{org}.s{seed_num}.f{f}".format(prefix = tool_prefix[tool], org=org, seed_num=seed_num, f=f)
                RM_dir = DIR + "/" + ("{org}.s{seed}.f{f}".format(org=org, seed=seed_num, f=f)).upper()
                RM_file = RM_dir + "/" + "{org}.fa.out".format(org=org, seed=seed_num, f=f)
                blast_file = RM_dir + "/" + "{org}.s{seed}.f{f}.blast.6.txt".format(org=org, seed=seed_num, f=f)
                pra_output = RM_dir + "/" + "{org}.s{seed}.f{f}.pra.txt".format(org=org, seed=seed_num, f=f)

                tool_output = RM_file
                real_repeats = data_map[org] + ".out"

                H = create_stats_hash(tool, org, int(seed_num), int(f))

                seed = convert_seed(seed_map[seed_num])
                seed_len = len(seed)
                seed_weight = seed.count("1")
                seed_ratio = seed_weight / (float(seed_len))


                H['l'] = seed_len
                H['w'] = seed_weight
                H['w/l'] = seed_ratio


                # Get stats from RM run
                try:
                    Counts, Stats, Sets = perform_stats.perform_stats(real_repeats, tool_output, None)
                    H['tp'], H['fp'], H['fn'], H['tn'] = Counts
                    H['tpr'], H['tnr'], H['ppv'], H['npv'], H['fpr'], H['fdr']  = Stats
                except Exception as E:
                    pass
                    #raise E;


                # Get resource usage from RAIDER run
                if os.path.exists(RAIDER_job_file):
                    p = redhawk.loadPBS(open(RAIDER_job_file, "rb"))[0]
                    try:
                        if p.efile_exists():
                            H['ToolCpuTime'], H['ToolWallTime'], H['ToolMem'], H['ToolVMem'] = p.getResources()
                    except:
                        pass
                    redhawk.storePBS([p], open(RAIDER_job_file, "wb"))

                # Get resource usage from RM run
                if os.path.exists(RM_job_file):
                    p = redhawk.loadPBS(open(RM_job_file, "rb"))[0]
                    try:
                        if p.efile_exists():
                            H['RMCpuTime'], H['RMWallTime'], H['RMMem'], H['RMVMem'] = p.getResources()
                    except:
                        pass
                    redhawk.storePBS([p], open(RM_job_file, "wb"))

                if os.path.exists(blast_file):
                    cmd = "./pra_analysis2 {blast_output} {output}".format(blast_output=blast_file, output=pra_output)
                    subprocess.call(re.split("\s+", cmd))
                    query_cover, target_cover, Used = parse_pra_output.parse_pra_output(pra_output, "exclude.txt")
                    H['ConCoverage'], H['QuCoverage'] = query_cover, target_cover


                stats_map[(tool,org,seed_num,f)] = H
    return None

def collectRptScout(DIR, tool):
    # FIRST: Get any applicable RM output stats
    for org in data_map.keys():
        for f in f_list:
            RS_job_file = DIR + "/../job_log/{prefix}.{org}.s0.f{f}".format(prefix = tool_prefix[tool], org=org, f=f)        
            RM_job_file = DIR + "/../job_log/rm.{prefix}.{org}.s0.f{f}".format(prefix = tool_prefix[tool], org=org, f=f)
            RS_dir = DIR + "/" + ("{org}.s0.f{f}".format(org=org, f=f)).upper()
            RM_file = RS_dir + "/" + "{org}.fa.out".format(org=org, f=f)
            blast_file = RS_dir + "/" + "{org}.s0.f{f}.RS.blast.6.txt".format(org=org, f=f)
            pra_output = "{DIR}/{org}.s0.f{f}.pra.txt".format(DIR=RS_dir, org=org, f=f)

            tool_output = RM_file
            real_repeats = data_map[org] + ".out"

            H = create_stats_hash(tool, org, None, int(f))


            # Get stats from RM run
            try:
                Counts, Stats, Sets = perform_stats.perform_stats(real_repeats, tool_output, None)
                H['tp'], H['fp'], H['fn'], H['tn'] = Counts
                H['tpr'], H['tnr'], H['ppv'], H['npv'], H['fpr'], H['fdr']  = Stats
            except Exception as E:
                pass
                #raise E;


            # Get resource usage from RPT_SCOUT run
            if os.path.exists(RS_job_file):
                p = redhawk.loadPBS(open(RS_job_file, "rb"))[0]
                try:
                    if p.efile_exists():
                        H['ToolCpuTime'], H['ToolWallTime'], H['ToolMem'], H['ToolVMem'] = p.getResources()
                except:
                    pass
                redhawk.storePBS([p], open(RS_job_file, "wb"))

            # Get resource usage from RM run
            if os.path.exists(RM_job_file):
                p = redhawk.loadPBS(open(RM_job_file, "rb"))[0]
                try:
                    if p.efile_exists():
                        H['RMCpuTime'], H['RMWallTime'], H['RMMem'], H['RMVMem'] = p.getResources()
                except:
                    pass
                redhawk.storePBS([p], open(RM_job_file, "wb"))

            if os.path.exists(blast_file):
                cmd = "./pra_analysis2 {blast_output} {output}".format(blast_output=blast_file, output=pra_output)
                subprocess.call(re.split("\s+", cmd))
                query_cover, target_cover, Used = parse_pra_output.parse_pra_output(pra_output, "exclude.txt")
                H['ConCoverage'], H['QuCoverage'] = query_cover, target_cover


            stats_map[(tool,org,None,f)] = H
    return None






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Collect stats. from testing_pipeline.run")
    parser.add_argument('-f', dest = 'f', type = int, nargs = '+', help = 'help values', default = None)
    parser.add_argument('-o', '--orgs', dest = 'data', type = str, nargs = '+', help = "orgs", default = None)
    parser.add_argument('-s', '--seed_list', dest = 'seeds', type = int, nargs = '+', default = None)
    parser.add_argument('-t', '--tool_list', dest = 'tools', nargs = '+', default = None)
    parser.add_argument('-b', '--bogus', dest = 'bogus', action = 'store_false', 
                        help = 'Nonsense argument for argument termination', default = False)
    parser.add_argument('--out', '--output', dest = 'output', type = str, default = 'stats.txt')
    parser.add_argument('DIR', help = 'base directory', default = False)

    args = parser.parse_args()
    DIR = args.DIR



    load_seeds(DIR + "/")
    load_map(DIR + "/")
    load_f_list(DIR + "/")


    if not args.tools or 'phRAIDER' in args.tools:
        collectRaider(DIR + "/" + "PHRAIDER", 'phRAIDER')
    if not args.tools or 'RepeatScout' in args.tools:
        collectRptScout(DIR + "/" + "RPT_SCT", "RepeatScout")


    with open(DIR + "/" + args.output, "w") if args.output != '-' else sys.stdout as fp:
        fp.write("\t".join(["{:<15}".format(x) for x in print_stats]) + "\n")
        for H in sorted(stats_map.values(), key = lambda x: (x["org"], x['f'], x["tool"], x["l"], x["w/l"])):
            line = "";
            for col in print_stats:
                v = H[col]
                if v is None:
                    v = "NA"
                elif type(v) == float:
                    v = round(v,3)
                line += "{:<15}\t".format(str(v))
            line += "\t" + seed_map[x['seed_num']]
            fp.write(line.rstrip() + "\n");
