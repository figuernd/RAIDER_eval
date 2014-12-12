import subprocess

chr_list = [21, 20, 19, 18, 17, 16]
cmd = "python3.3 RAIDER_eval.py --nuke -r RESULTS.len4/chr{chr}.sim1.all -R --sf seed.len4.txt --mem chrom_sim --st 0 data/chr{chr}.fa > nohup.{chr}.out"

P = []
for i in chr_list:
    c = cmd.format(chr=i)
    print(c)
    P.append(subprocess.Popen(c, shell = True))

for p in P:
    p.wait()

