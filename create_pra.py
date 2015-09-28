import sys
import redhawk
import glob
import os

reg = re.compile("^(.*)\//([^\/]()\.s(\d+)\.f(\d+)\.blast.6.txt(\.bz2)?$")

def main(DIR):
    for file in glob.glob("{DIR}/*/*/*.blast.6.txt*".format(DIR)):
        r = reg.search(file)
        if not r:
            sys.stderr.write("Bad file: " + file + "\n")
            exit(1);
        D = r.group(1)
        org = r.group(2)
        s = int(r.group(3))
        f = int(r.group(4))

        pra_output =  D + "/" + "{org}.s{seed}.f{f}.pra.txt".format(org=org, seed=seed_num, f=f)

        if os.file.exists(pra_output):
            sys.stdout.write("Skipping: " + pra_output + "\n")
            continue

        if file.endswith(".bz2"):
            cmd = "bzcat {blast_output} | ./pra_analysis2 {output}".format(blast_output=file, output=pra_output)
        else:
            cmd = "cat {blast_output} | ./pra_analysis2 {output}".format(blast_output=file, output=pra_output)

        sys.stdout.write("CMD: " + cmd)
        
if __name__ == "__main__":
    main(sys.arv[1])
