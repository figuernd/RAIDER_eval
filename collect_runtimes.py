import glob
import re

ToolsDic = {'RS':'RepeatScout', 'phRA':'phRAIDER', 'RA':'RAIDER', 'prephRA':'pre-phRAIDER'}
SizeDic = {x:y for line in open("sizes.txt") for x,y in [re.split("\s+", line.rstrip())]}

def readFile(file):
   for line in open(file):
      r = re.search("User time \(seconds\):\s+(\d+\.\d+)", line)
      if r:
         time=int(r.group(1))
         continue
      r = re.search("Maximum resident set size \(kbytes\):\s+(\d+)", line)
      if r:
         mem = int(r.group(1))/4
         return (time,mem)
      return ('NA', 'NA')


R = {}
for file in glob.glob("*/*/*.stderr"):
   if "add_files" in file:
      continue
   r = re.search("(\w+)\.(\w+)\.(chr\w+)\.s(\d+)\.f(\d+)\.stderr$", file);
   if not r:
      continue
   tool, build, chr, seed, f = r.group(1,2,3,4,5)
   key = (tool,build,chr,seed,f)

   
   label = build + "." + chr + "." = 
   if build + "." + chr not in D:
      continue

   
   time,mem = readFile(file)
   if time == 'NA' or mem == 'NA':
      next

   if key not in R:
      R[key] = (time,mem)
   else:
      t,m = R[key]
      R[key][0] += t
      R[key][1] = min(R[key][1], m)

for key in sorted(R.key):
   print("\t".join([tool, build + "." + chr, SizeDic[build + "." + chr]
