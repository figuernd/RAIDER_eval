#!/software/python/3.3.3/bin/python3.3
import sys
import subprocess
import os
import os.path
import shutil
import argparse
from redhawk import *
import tempfile
import re
import perform_stats
import time
import pickle


flist_start = "START_FILE_LIST"
flist_end = "END_FILE_LIST"
csjobs_start = "START_CHROM_SIM_JOBS"
csjobs_end = "END_CHROM_SIM_JOBS"
tjobs_start = "START_TOOL_JOBS"
tjobs_end = "END_TOOL_JOBS"
rmjobs_start = "START_REPMASK_JOBS"
rmjobs_end = "END_REPMASK_JOBS"
prajobs_start = "START_PRA_JOBS"
prajobs_end = "END_PRA_JOBS"
jobdic_start = "START_JOB_DICT"
jobdic_end = "END_JOB_DICT"
blast_db_start = "START_BLAST_DB"
blast_db_end = "END_BLAST_DB"
stats_start = "START_STATS_JOBS"
stats_end = "END_STATS_JOBS"



def flush_files():
    """If we submitted the evaluation as a PBS job with a set walltime, call this as a debug measure
    to ensure that get updated log files if program quits unexpectedly"""
    logging_fp.flush()
    checkpoint_fp.flush()

def write_flist_to_checkpoint(flist):
    """If we submitted the evaluation as a PBS job with a set walltime, write list of files (either supplied seq
    files or generated chromosome files) to checkpoint file for access upon next run"""
    logging_fp.write("Writing file list to new checkpoint file\n")
    checkpoint_fp.write(flist_start+"\n")
    for f in flist:
        checkpoint_fp.write(f+"\n")
    checkpoint_fp.write(flist_end+"\n")
    flush_files()

def write_jobs_to_checkpoint(jobs, start_id, end_id):
    """If we submitted the evaluation as a PBS job with a set walltime, pickle PBS jobs and store the path
    to the pickle file name in checkpoint file for unpickling upon next run. Note: start_id and end_id are
    indicators when parsing checkpoint file to tell program what kind of object was stored at a location"""
    logging_fp.write("Writing jobs to new checkpoint file\n")
    checkpoint_fp.write(start_id+"\n")
    for ufj in jobs:
        pick_fname = ufj.batch_file_name + ".pickle"
        checkpoint_fp.write(pick_fname+"\n")
        storePBS(ufj, open(pick_fname, "wb"))
    checkpoint_fp.write(end_id+"\n")
    flush_files()

def write_job_dict_to_checkpoint(job_dic, results_dir): 
    """If we submitted the evaluation as a PBS job with a set walltime, pickle list of PBS jobs for each
    type of tool being evaluated and store path to pickle file name in checkpoint file for unpickling
    upon next run."""
    logging_fp.write("Writing job dict to new checkpoint file\n")
    checkpoint_fp.write(jobdic_start + "\n")
    for tool in job_dic:
        pick_name = "{dir}/job_dic.{tname}.pickle".format(dir=results_dir, tname=tool)
        checkpoint_fp.write(pick_name+"\n")
        storePBS(job_dic[tool], open(pick_name, "wb"))
    checkpoint_fp.write(jobdic_end + "\n")
    flush_files()

def write_blast_db_to_checkpoint(blast_db, results_dir): 
    """If we submitted the evaluation as a PBS job with a set walltime, pickle list of PBS jobs for each
    seq_file in blast database and store path to pickle file name in checkpoint file for unpickling
    upon next run."""
    logging_fp.write("Writing job dict to new checkpoint file\n")
    checkpoint_fp.write(blast_db_start + "\n")
    for seq_file in blast_db:
        pick_name = "{dir}/blast_db.{fname}.pickle".format(dir=results_dir, fname=seq_file) #CARLY: worried that fname might have dirs in it...
        checkpoint_fp.write(pick_name+"\n")
        storePBS(blast_db[seq_file], open(pick_name, "wb"))
    checkpoint_fp.write(blast_db_end + "\n")
    flush_files()
    
def save_timed_chrom_sim_jobs(jobs, finished_jobs, flist):
    """If we submitted the evaluation as a PBS job with a set walltime, write list of simulated
    chromosome files that were completed to checkpoint file, and pickle any unfinished 
    chromosome simulation PBS jobs and store paths to pickle file names in checkpoint file"""
    flist.extend([x.sim_output for x in finished_jobs])
    write_flist_to_checkpoint(flist)
    write_jobs_to_checkpoint(jobs - finished_jobs, csjobs_start, csjobs_end)

def save_timed_tool_jobs(jobs, RM_jobs, PRA_jobs):
    """Helper method. Given a list of jobs and repeat masker jobs, pickle jobs and write
    pickled file paths to checkpoint for unpickling upon next run."""
    write_jobs_to_checkpoint(jobs, tjobs_start, tjobs_end)
    write_jobs_to_checkpoint(PRA_jobs, prajobs_start, prajobs_end)
    write_jobs_to_checkpoint(RM_jobs, rmjobs_start, rmjobs_end)
    
def save_timed_PRA_jobs(PRA_jobs):
    """Helper method. Given a list of jobs and job dictionary, pickle jobs and job dict,
    and write pickled file paths to checkpoint for unpickling upon next run."""
    write_jobs_to_checkpoint(PRA_jobs, prajobs_start, prajobs_end)

def save_timed_RM_jobs(jobs, stats_jobs, results_dir, job_dic):
    """Helper method. Given a list of jobs and job dictionary, pickle jobs and job dict,
    and write pickled file paths to checkpoint for unpickling upon next run."""
    write_jobs_to_checkpoint(jobs, rmjobs_start, rmjobs_end)
    write_jobs_to_checkpoint(stats_jobs, stats_start, stats_end)
    write_job_dict_to_checkpoint(job_dic, results_dir)

def recover_job_list(job_type, end_flag):
    logging_fp.write("Reading in {type} jobs from old checkpoint file\n".format(type=job_type))
    job_list = []
    pickname = old_checkpoint_fp.readline().rstrip()
    while pickname != end_flag:
        j = loadPBS(open(pickname, "rb"))
        job_list.append(j)
        pickname = old_checkpoint_fp.readline().rstrip()
    next_step = old_checkpoint_fp.readline().rstrip()
    return job_set, next_step

def recover_job_set(job_type, end_flag):
    logging_fp.write("Reading in {type} jobs from old checkpoint file\n".format(type=job_type))
    job_set = set()
    pickname = old_checkpoint_fp.readline().rstrip()
    while pickname != end_flag:
        j = loadPBS(open(pickname, "rb"))
        job_set.add(j)
        pickname = old_checkpoint_fp.readline().rstrip()
    next_step = old_checkpoint_fp.readline().rstrip()
    return job_set, next_step

def recover_file_list(): 
    # Get file list from old checkpoint file
    file_list = []
    logging_fp.write("Reading file list from old checkpoint file\n")
    fname = old_checkpoint_fp.readline().rstrip()
    while fname != flist_end:
        file_list.append(fname)
        fname = old_checkpoint_fp.readline().rstrip()
    next_step = old_checkpoint_fp.readline().rstrip() 
    return file_list, next_step

def recover_sim_jobs():
    return recover_job_set("chrom sim", csjobs_end)
    
def recover_blast_db():
    logging_fp.write("Reading in blast database from old checkpoint file\n")
    BLAST_DATABASE = {}
    pickname = old_checkpoint_fp.readline().rstrip()
    while pickname != blast_db_end:
        seq_fname = os.path.splitext(os.path.splitext(pickname)[0])[1] #CARLY: still unsure of whether this will work.. maybe need to store basename of seq_file
        BLAST_DATABASE[seq_fname] = loadPBS(open(pickname, 'rb'))
        pickname = old_checkpoint_fp.readline().rstrip()
    next_step = old_checkpoint_fp.readline().rstrip()
    return BLAST_DATABASE, next_step

def recover_tool_jobs():
    return recover_job_list("tool", tjobs_end)

def recover_rm_jobs():
    return recover_job_set("repmask", rmjobs_end)
    
def recover_pra_jobs():
    return recover_job_set("pra", prajobs_end)

def recover_stats_jobs():
    return recover_job_set("stats", stats_end)

def recover_job_dic():
    logging_fp.write("Reading in job dict from old checkpoint file\n")
    old_job_dic = {tool:[] for tool in test_tools}
    pickname = old_checkpoint_fp.readline().rstrip()
    while pickname != jobdic_end:
        tname = os.path.splitext(os.path.splitext(pickname)[0])[1]
        old_job_dic[tname] = loadPBS(open(pickname, 'rb'))
        pickname = old_checkpoint_fp.readline().rstrip()
    flush_files()
    next_step = old_checkpoint_fp.readline().rstrip()


