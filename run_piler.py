import sys
import subprocess
import re
import os
import glob

"""Helper script for running piler"""

if __name__ == "__main__":
    seq_file_dir = sys.args[1]
    seq_file = sys.argv[2]
    piler_dir = sys.argv[3]

    gff_file = re.sub("\.fa$|\.fasta$", ".gff", seq_file)
    trs_file = re.sub("\.fa$|\.fasta$", "trs.gff", seq_file)
    fam_dir = piler_dir + "/" + re.sub("\.fa$|\.fasta$", seq_file)
    output_lib = re.sub("\.fa$|\.fasta$", ".lib", seq_file)

    os.mkdir(fam_dir)

    cmd1 = "./pals -self {seq_file} {gff_file}"
    cmd2 = "./piler2 -trs {gff_file} -out {trs_file}"
    cmd3 = "./piler2 --trs2fasta {trs_file} -seq {seq_file} -path {fam_dir}"
    cmd = "; ".join([cmd1, cmd2, cmd3]).format(seq_file=seq_file, gff_file=gff_file, trs_file = trs_file, fam_dir=fam_dir)
    print(cmd);

    cmd = ""
    for fam in glob(fam_dir + "/*"):
        cmd += "./muscle -in {fam} -out {fam}.aligned -maxiters 1 -diags1; ".format(fam=fam)
        cmd += "./piler2 -cons {fam}.aligned -out {fam}.cons; ".format(fam=fam)
    cmd += "cat {fam_dir}/*.cons > {piler_dir}/{output_file}".format(fam_dir = fam_dir, piler_dir = piler_dir, output_file = output_lib)
    print(cmd)
    #subprocess.call(cmd, shell = True)

    

    
