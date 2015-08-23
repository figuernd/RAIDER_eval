import sys
import perform_stats
import glob
import re
import redhawk
import parse_pra_output
import os
import subprocess

tool_prefix = {'phRAIDER':"phRA", "RepeatScout":"RS"}

seed_map = None
data_map = None
f_list = None
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
    H["seed_num"] = int(seed_num)
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
        seed_map[int(index)] = seed


def load_map(DIR):
    global data_map
    data_map = {}
    for line in open(DIR + "/data_files.txt"):
        short_name, long_name = re.split("\s+", line.strip())
        data_map[short_name] = long_name

def load_f_list(DIR):
    global f_list
    f_list = [f for line in open("DIR + /f.txt") for f in re.split("\s+", line.strip())]

        
regex = re.compile("S(\d+)\.F(\d+)\.RM\/(.*).fa.out$")
def collectRaider(DIR, tool):
    # FIRST: Get any applicable RM output stats
    for org in data_map.keys():
        for seed_num in seed_map.keys():
            for f in f_list:
                RAIDER_job_file = DIR + "/../job_log/{prefix}.{org}.{seed_num}.{f}".format(prefix = tool_prefix[tool], org=org, seed_num=seed_num, f=f)        
                RM_job_file = DIR + "/../job_log/rm.{prefix}.{org}.{seed_num}.{f}".format(prefix = tool_prefix[tool], org=org, seed_num=seed_num, f=f)        
                RM_dir = ("{org}.s{seed}.f{f}.rm".format(org=org, seed=seed_num, f=f)).upper
                RM_file = RM_dir + "/" + "{org}.{seed}.{f}.fa.out".format(org=org, seed=seed, f=f)
                blast_file = "{org}.s{seed}.f{f}.blast.6.txt".format(org=org, seed=seed, f=f)
                pra_output = "{DIR}/{org}.s{seed}.f{f}.pra.txt".format(DIR=DIR, org=org, seed=seed_num, f=f)


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
                    output = 
                    cmd = "./pra_analysis2 {blast_output} {output}".format(blast_output=blast_file, output=pra_output)
                    query_cover, target_cover, Used = parse_pra_output.parse_pra_output(pra_output, "exclude.txt")
                    H['ConCoverage'], H['QuCoverage'] = query_cover, target_cover


                stats_map[(tool,org,seed_num,f)] = H
        return None

regex2 = re.compile("\./(.*).fa.out$")
def collectRptScout(DIR, tool):
    # FIRST: Get any applicable RM output stats
    for file in glob.glob(DIR + "/*.fa.out"):
        print("HERE: ", file)
        r = regex2.search(file)
        assert(r)
        seed_num, f, org = r.group(1,2,3)
        f = int(f)
        tool_output = file
        real_repeats = data_map[org] + ".out"

        H = create_stats_hash(tool, org, int(seed_num), int(f))

        # Get stats from RM run
        try:
            Counts, Stats, Sets = perform_stats.perform_stats(real_repeats, tool_output, None)
            H['tp'], H['fp'], H['fn'], H['tn'] = Counts
            H['tpr'], H['tnr'], H['ppv'], H['npv'], H['fpr'], H['fdr']  = Stats
        except Exception as E:
            pass
            #raise E;


        # Get resource usage from RAIDER run
        job_file = DIR + "/../job_log/{prefix}.{org}.{f}".format(prefix = tool_prefix[tool], org=org, seed_num=seed_num, f=f)        
        if os.path.exists(job_file):
            p = redhawk.loadPBS(open(job_file, "rb"))[0]
            try:
                if p.efile_exists():
                    H['ToolCpuTime'], H['ToolWallTime'], H['ToolMem'], H['ToolVMem'] = p.getResources()
            except:
                pass
            redhawk.storePBS([p], open(job_file, "wb"))

        # Get resource usage from RM run
        job_file = DIR + "/../job_log/rm.{prefix}.{org}.{f}".format(prefix = tool_prefix[tool], org=org, seed_num=seed_num, f=f)        
        if os.path.exists(job_file):
            p = redhawk.loadPBS(open(job_file, "rb"))[0]
            try:
                if p.efile_exists():
                    H['RMCpuTime'], H['RMWallTime'], H['RMMem'], H['RMVMem'] = p.getResources()
            except:
                pass
            redhawk.storePBS([p], open(job_file, "wb"))

        blast_output = "{DIR}/{org}.f{f}.blast.6.txt".format(DIR=DIR, org=org, seed=seed_num, f=f)
        if os.path.exists(blast_output):
            output = "{DIR}/{org}.f{f}.pra.txt".format(DIR=DIR, org=org, seed=seed_num, f=f)
            cmd = "./pra_analysis2 {blast_output} {output}".format(blast_output=blast_output, output=output)
            subprocess.call(re.split("\s+", cmd))
            query_cover, target_cover, Used = parse_pra_output.parse_pra_output(output, "exclude.txt")
            H['ConCoverage'], H['QuCoverage'] = query_cover, target_cover
            peint_stats_hash(H)

        stats_map[(tool,org,seed_num,f)] = H
        
        return None



if __name__ == "__main__":
    DIR = sys.argv[1]
    load_seeds(DIR + "/")
    load_map(DIR + "/")
    load_f_list(DIR + "/")

    for f in glob.glob(DIR+"/*"):
        if f.endswith('PHRAIDER'):
            collectRaider(f, 'phRAIDER')
        #elif f.endswith('RPT_SCT'):
        #    collectRptScout(f, "RepeatScout")


    with open(DIR + "/stats.txt", "w") as fp:
        fp.write("\t".join(["{:<15}".format(x) for x in print_stats]) + "\n")
        for H in stats_map.values():
            line = "";
            for col in print_stats:
                v = H[col]
                if v is None:
                    v = "NA"
                elif type(v) == float:
                    v = round(v,3)
                line += "{:<15}\t".format(str(v))
            fp.write(line.rstrip() + "\n");
