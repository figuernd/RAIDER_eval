import subprocess

chr_list = [21, 19, 18, 15, 12, 11, 10]
cmd = "python RAIDER_eval.py --nuke -r RESULTS.seed34/chr{chr}.sim1.all --A2 --sf seed34.txt --mem chrom_sim data/chr{chr}.fa > nohup.{chr}.out"

P = []
for i in chr_list:
    P.append(subprocess.Popen(cmd.format(chr=i), shell = True))

for p in P:
    P.wait()

