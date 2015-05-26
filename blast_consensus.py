from Bio.Blast.Applications import NcbiblastxCommandline
import sys
import re
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
import subprocess
import os.path
from coverageBLAST import coverageBLAST

e_val = 0.0001   # Needed constant for blast invocation.  Specifies the expected number of false positive matches.


def invokeBlast(query_file, db_file, out_file = None):
    """Blast each ancestor against the elements database; return name of output file"""
    output = out_file if out_file else query_file + ".xml"
    blastn_cline = NcbiblastnCommandline(query=query_file, db=db_file, outfmt=5, out=out_file if out_file else query_file + ".xml", evalue = e_val)
    stdout, stderr = blastn_cline()

    #assert not stderr, "BLAST: " + stderr

    return output



def main(consensus_file, rm_fa_file, database_file, output_file):
    """Blast the consensus sequences against known instance of repeats to calculate coverage.
    * consensus_file: Contains the consensus sequences created by the tool.
    * rm_fa_file: A fasta file of the know repeat instances.
    * database_file: A blast database created from rm_fa_file (already created)
    * Output_file: Where to put the output.
    """
    #### First: run blast
    blast_output_file = invokeBlast(consensus_file, database_file, consensus_file.rstrip(".fa") + ".blast.xml")

    #### Second: Check coverage of each seqeunce
    fp = open(blast_output_file)
    blast_records = NCBIXML.parse(fp)

    consensus_map = {}
    rpt_map = {}
    for blast_obj in blast_records:
        consensus_id = re.match("(\S+)", blast_obj.query).group(1)
        if not consensus_id in consensus_map:
            consensus_map[consensus_id] = set()

        for align_obj in blast_obj.alignments:
            rpt_id = re.match("\S+\s+(\S+)", align_obj.title).group(1)

            if not rpt_id in rpt_map:
                rpt_map[rpt_id] = set()

            for hsp_obj in align_obj.hsps:
                consensus_start = hsp_obj.query_start
                consensus_finish = hsp_obj.query_end
                if consensus_finish < consensus_start:
                    consensus_start, consensus_finish = consensus_finish, consensus_start
                consensus_map[consensus_id] |= set(range(consensus_start-1, consensus_finish))

                subject_start = hsp_obj.sbjct_start
                subject_finish = hsp_obj.sbjct_end
                if subject_finish < subject_start:
                    subject_finish, subject_start = subject_start, subject_finish
                rpt_map[rpt_id] |= set(range(subject_start-1,subject_finish))

    #### Third: Calculate coverage of consensus sequences.
    #### At the end, for each consensus_sequence id, we will have a tuple (c,l), where
    #### c is the number of bases covered, and l is the number of bases.
    consensus_coverage = {}
    for seq_record in SeqIO.parse(consensus_file, 'fasta'):
        if seq_record.id in consensus_map:
            consensus_coverage[seq_record.id] = len(consensus_map[seq_record.id]), float(len(seq_record))
        else:
            consensus_coverage[seq_record.id] = 0, float(len(seq_record))

    #### Fourth: Calculate coverage each rpt sequences.
    #### At the end, for each rpt_sequence instance, we will have a tuple (c,l), where
    #### c is the number of bases covered, and l is the number of bases.
    rpt_coverage = {}
    for seq_record in SeqIO.parse(rm_fa_file, 'fasta'):
        if seq_record.id in rpt_map:
            rpt_coverage[seq_record.id] = len(rpt_map[seq_record.id]), float(len(seq_record))
        else:
            rpt_coverage[seq_record.id] = 0, float(len(seq_record))



    #######################################
    #### Finally: Print results
    wp = sys.stdout if output_file == "-" else open(output_file, "w")


    wp.write("# Consensus coverage (id, num covered bases, length, ratio)\n")
    sum, total = 0.0, 0
    for k,v in sorted(consensus_coverage.items()):
        f = v[0]/v[1] if v[1] > 0 else -1
        sum += f
        total += 1
        wp.write("{id}\t{covered}\t{length}\t{frac}\n".format(id = k, covered = int(v[0]), length = int(v[1]), frac = round(f,4)))
    wp.write("# Average consensus coverage: {avg}\n".format(avg = -1 if total == 0 else round(sum / total, 4)))

    wp.write("# Query covereage (id, num covered bases, length, ratio)\n")
    D = {}
    for id,v in sorted(rpt_coverage.items()):
        fam = id[:id.find("|")]
        f = v[0]/v[1] if v[1] > 0 else 0
        if fam not in D:
            D[fam] = (0,0)
        wp.write("{rpt:<30}\t{covered:<6}\t{length}\t{frac:<6}\n".format(rpt = id, covered = int(v[0]), length = int(v[1]), frac = round(f,4)))
        D[fam] = D[fam][0] + f, D[fam][1] + 1
        
    wp.write("# Average family coverage (family, average coverage) -- not sure this statistic is meaningful\n")
    for fam,T in sorted(D.items()):
        wp.write("{fam}\t{avg}\n".format(fam=fam, avg=round(T[0]/T[1] if T[1] > 0 else -1,4)))
    


if __name__ == "__main__":
    consensus_file = sys.argv[1]
    rm_fa_file = sys.argv[2]
    database_file = sys.argv[3]
    output_file = sys.argv[4]
                
    main(consensus_file, rm_fa_file, database_file, output_file);

