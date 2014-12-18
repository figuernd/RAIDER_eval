import glob
import re

wp = open("collected_stats.txt", "w")
headings = ["chr", "tool", "seed", "tpr", "tnr", "ppv", "npv" ,"fpr", "fdr", "Time", "Memory", "RM_Time", "RM_Mem"]
wp.write(("{:<14}"*len(headings)).format(*headings) + "\n")
def main():
    for file in glob.glob("chr*/stats.txt"):
        chr  = re.match("(chr\d+)", file).group(1)
        with open(file) as fp:
            fp.readline()
            for line in fp:
                A = re.split("\s+", line.rstrip())
                if len(A) != 3:
                    O = A[0], A[1], A[6], A[7], A[8], A[9], A[10], A[11], A[12], A[14], A[16], A[18]
                    wp.write(("{:<14}"*(len(O)+1)).format(chr, *O) + "\n")

if __name__ == "__main__":
    main()
