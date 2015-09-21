import glob
import re

ToolsDic = {'RS':'RepeatScout', 'phRA':'phRAIDER', 'RA':'RAIDER', 'prephRA':'pre-phRAIDER'}
SizeDic = {x:y for line in open("sizes.txt") for x,y in [re.split("\s+", line.rstrip())]}

def readFile(file):
   for line in open(file):
      r = re.search("User time \(seconds\):\s+(\d+\.\d+)", line)
      if r:
         time= float(r.group(1))
         continue
      r = re.search("Maximum resident set size \(kbytes\):\s+(\d+)", line)
      if r:
         mem = float(r.group(1))/4
         return (time,mem)
   return ('NA', 'NA')


R = {}
for file in glob.glob("*/*/*.stderr"):
   if "add_files" in file:
      continue
   r = re.search("(\w+(?:\.\w+))\.(\w+)\.(chr\w+)\.s(\d+)\.f(\d+)\.stderr$", file);
   if not r:
      continue
   tool, build, chr, seed, f = r.group(1,2,3,4,5)
   if build + "." + chr not in SizeDic:
      continue
   key = (tool,build,chr,seed,f)

      
   time,mem = readFile(file)
   if time == 'NA' or mem == 'NA':
      continue

   if key not in R:
      R[key] = (time,mem)
   else:
      t,m = R[key]
      R[key] = (t+time, max(m,mem))

print("\t".join(["tool", "label", "seed", "f", "size", "time", "mem"]) + "\n")
for key in sorted(R.keys()):
   print("\t".join([ToolsDic[key[0]], key[1] + "." + key[2], str(key[3]), str(key[4]), SizeDic[key[1] + "." + key[2]], str(R[key][0]), str(R[key][1])]))
