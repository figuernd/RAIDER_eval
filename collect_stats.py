import sys
import perform_stats
import glob
import re
import redhawk.py

seed_map = None
data_map = None
stats_map = []

def load_seeds(DIR):
    global seed_map
    seed_map = {}
    for line in open(DIR + "/seed.txt"):
        index, seed = re.split("\s+", line.strip())
        seed_map[index] = seed
        
    return seed_map

def load_map(DIR):
    global data_map
    data_map = {}
    for line in open(DIR + "/data_files.txt"):
        short_name, long_name = re.split("\s+", line.strip())
        data_map[short_name] = long_name
        
regex = re.compile("S(\d+)\.F(\d+).rm/(.*).fa.out")
def collectRaider(DIR):
    # FIRST: Get any applicable RM output stats
    for file in glob(DIR + "/*.RM/*.fa.out"):
        r = regex.search(file)
        seed_num, f, org = r.group(1,2,3)
        seed = seed_map[int(seed_num)]
        f = int(f)
        tool_output = file
        real_repeats = data_map[org]
        try:
            Counts, Stats, Sets = perform_stats(real_repeats, tool_output, None)
        except:
            Stats = []

        try:
            p = loadPBS(open("phRA." + org + "." + str(seed_num) + "." + str(f)))
            Resources = list(p.getResources(clean=False))
        except:
            Resources = []
                    
        
        stats_map[(org,s,f)] = [Stats, Resources)

        

if __name__ == "__main__":
    DIR = sys.argv[1]
    load_seeds(DIR)
    
    for f in glob.glob("*"):
        if f.endswith('PHRAIDER')
            collectRaider(f)
        #elif f.sendswith('RPT_SCT'):
        #    collectRptsct(f)
            

