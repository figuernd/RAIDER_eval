import subprocess

chr_list = [21, 19, 18, 15, 12, 11, 10]
cmd = "python RAIDER_eval.py --nuke -r RESULTS.seedlen/chr{chr}.sim1.all --R --sf seed.len.txt --mem chrom_sim data/chr{chr}.fa > nohup.{chr}.out"

P = []
for i in chr_list:
    c = cmd.format(chr=i)
    print(c)
    P.append(subprocess.Popen(c, shell = True))

for p in P:
    P.wait()

