import glob
import re

TD = {'RS':'RepeatScout', 'phRA':'phRAIDER', 'RA':'RAIDER', 'prephRA':'pre-phRAIDER'}

D = {x:y for line in open("sizes.txt") for x,y in [re.split("\s+", line.rstrip())]}

for file in glob.glob("*/*/*.stderr"):
   if "add_files" in file:
      continue
   r = re.search("(\w+)\.(\w+)\.(chr\w+)\.s(\d+)\.f(\d+)\.stderr$", file);
   if not r:
      continue
   tool, build, chr, seed, f = r.group(1,2,3,4,5)
   if build + "." + chr not in D:
      continue

   time = "NA"
   mem = "NA"
   for line in open(file):
      r = re.search("User time \(seconds\):\s+(\d+\.\d+)", line)
      if r:
         time=r.group(1)
         continue
      r = re.search("Maximum resident set size \(kbytes\):\s+(\d+)", line)
      if r:
         mem = int(r.group(1))/4

   print("\t".join([TD[tool], build, chr, build + "." + chr, seed, D[build + "." + chr], str(time), str(mem)]))
