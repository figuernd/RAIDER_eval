import subprocess

chr_list = [21, 20, 19, 18, 15, 12, 11]
cmd = "python3.3 RAIDER_eval.py --nuke -r RESULTS.len3/chr{chr}.sim1.all -R --sf seed.len3.txt --mem chrom_sim --st 0 data/chr{chr}.fa > nohup.{chr}.out"

P = []
for i in chr_list:
    c = cmd.format(chr=i)
    print(c)
    P.append(subprocess.Popen(c, shell = True))

for p in P:
    p.wait()

