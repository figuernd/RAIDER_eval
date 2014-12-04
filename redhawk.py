#!/usr/bin/python2.7

############################################################################
# redhawk.py
# Class for submitting / managing redhawk jobs through a python script
# Written by: John Karro, Jens Muller
# Date: Nov. 1, 2011
# Last modified: July 12, 2014
#

import random
import re, string, sys
import subprocess
import getopt
import time
import os
import os.path
import pickle
from argparse import *
import pwd
import datetime
import pickle
import tempfile


epilogue_str = """#!/bin/sh
echo "Redhawk Epilogue Args:" >&2
echo "Job ID: $1" >&2
echo "User ID: $2" >&2
echo "Group ID: $3" >&2
echo "Job Name: $4" >&2
echo "Session ID: $5" >&2
echo "Resource List: $6" >&2
echo "Resources Used: $7" >&2
echo "Queue Name: $8" >&2
echo "Account String: $9" >&2
echo "" >&2
exit 0
"""

bad_ep_str = """
Redhawk Epilogue Args:
Resources Used: cput=00:00:00,mem=0kb,vmem=0kb,walltime=00:00:00
"""

############################################################
# Library-global variables
log_file = "/usr/local/torque/current/var/spool/torque/server_priv/accounting"
uid = os.getuid()
current_user = pwd.getpwuid(uid)[0]  # When run on redhawk with a nohup, os.getlogin() does not work
try:
    subprocess.check_call(["qstat"], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
except:
    pbs_present = False
else:
    pbs_present = True

job_list = set()    # Global list of all jobs that have been submitted and not yet identified as having quit

job_limit = 200 if pbs_present else 2

redhawkStatsRe = re.compile("\s+C\s+[^\s]+\s*$")
redhawkInQueueRe = re.compile("\s+[R|Q]\s+[^\s]+\s*$")

# This contains the default parameter values for __init__.  Doing it this way is no longer necessary,
# but its not worth taking out.
pbs_defaults = {'use_pid':True, 'job_name':None, 'nodes':1, 'ppn':1, 'mem':False, 'walltime':"04:00:00", 'address':None, 'join':False, 'env':None, 'queue':None, 'mail':None, 'output_location':None, 'chdir':None, 'RHmodules':None, 'file_limit':6, 'file_delay':5, 'epilogue_file':None}




############################################################
class PBSError(Exception):
    """Exception class for pbsJobHandler class"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


############################################################
class pbsJobHandler:
    """A pbsJobHandler corresponds to a job launched (or to be launched) on redhawk.  Once the object is created (and provided with a command-line execution command),
       the user can extract various inforamtion about the job (current status, output, etc...) and cleanup files."""
    def __init__(self, batch_file, executable, use_pid = pbs_defaults['use_pid'], job_name = pbs_defaults['job_name'], nodes = pbs_defaults['nodes'], ppn = pbs_defaults['ppn'], 
                 mem = pbs_defaults['mem'], walltime = pbs_defaults['walltime'], address = pbs_defaults['address'], join = pbs_defaults['join'], env = pbs_defaults['env'], 
                 queue = pbs_defaults['queue'], mail = pbs_defaults['mail'], output_location = pbs_defaults['output_location'], chdir = pbs_defaults['chdir'], 
                 RHmodules = pbs_defaults['RHmodules'], file_limit = pbs_defaults['file_limit'], file_delay = pbs_defaults['file_delay'], epilogue_file = pbs_defaults['epilogue_file'],
                 suppress_pbs = None, stdout_file = None, stderr_file = None, arch_type = None):
        """Constructor.  Requires a file name for the batch file, and the execution command.  Optional parmeters include:
           * use_pid: will embded a process id into the batch file name if true.  Default = true.
           * job_name: A name for the redhawk process.  Default = the batch file name.
           * nodes: number of nodes required for the job.   Default = 1.
           * ppn: number of processors needed for the job.  Default = 1.
           * mem: Using 128 Gb machine.  Default = False.
           * walltime: Maximum allowed runtime for the job (hours:minutes:seconds).  Default = 40:00:00.  Max. allowed: 400:00:00.
           * mail = when to send email.  Any combination of:
             b   send mail when job begins
             e   send mail when job ends
             a   send mail when job aborts
           * address: additional email addresses to send notification (comma seperated)
           * join: If true, the stdout and stderr files are combined into one file
           * queue: redhawk queue to run on.  Default: redhawk chooses.
           * output_location: Directory to place output files.  Default: current directory.
           * RHmodules: A list of redhawk modules to be loaded before run (e.g. ['Blast+']).  Default: none.
           * epilogue file: Script needed to track memory usage.  Will overwrite any file of the same name.  By default: <batch_file>.epilogue.py"
           * suppress_pbs: If true, this will run the job with a standard popen, instead of forking it out to the pbs job manager.  (Included so we can run code
                           on other machines for testing.) Will still create all files that would have been created -- emulates redhawk execution as much as possible.
                           Not yes set up to extract resource usage.  If false, will force an attempt to use the pbs job manager.  By default: it will use pbs if
                           possible (specifically: if qstat can be run on the machine)
           * stdout_file: File to receive stdout content.  (Default: <job_name>.o<id>, in output_location directory if specified.)
           * stdin_file: File to receive stderr content. (Default: <job_name.e<id>, in output_location directory if sepcified..)
           * arch_type: Array if specific redhawk architecture to be used (n09, n11, bigmem, 
        """
        if epilogue_file and "/" in epilogue_file:
            raise PBSError("Bad epilogue file name: " + epilogue_file)
        
        self.batch_file_name = batch_file
        if use_pid:
            self.batch_file_name = self.batch_file_name + "." + str(os.getpid())
        self.cmd = executable
        self.jobname = job_name if job_name else batch_file
        self.file_limit = file_limit
        self.file_delay = file_delay
        self.resources = None
        self.status = "unstarted"
        self.suppress_pbs = not pbs_present if suppress_pbs is None else suppress_pbs 
        self.epilogue = os.getcwd() + "/" + "redhawk_epilogue.py"
        f = open(self.batch_file_name, 'w')

        f.write("#!/bin/bash -l\n")
        s="#PBS -N "+ self.jobname +"\n"
        f.write(s)

        self.has_split = False    # Indicate whether the err file has ever been split (added for debugging)
     
        #some defaults:
        self.nodes = nodes
        self.ppn = ppn
        self.mem = mem
        self.walltime = walltime
        self.modules = RHmodules if not self.suppress_pbs else None
        self.output_location = output_location if output_location else "."
        self.arch_type = arch_type if arch_type else []

        s="#PBS -N " + self.jobname + "\n"
        s="#PBS -l nodes="+ str(self.nodes)+":ppn="+str(self.ppn)
        if self.mem:
            s += ":128"
        if self.arch_type:
            s += ":" + ":".join(arch_type)
        s += "\n"
        f.write(s)
        s="#PBS -l walltime="+self.walltime+"\n"
        f.write(s)
        
        if stdout_file:
            self.ofile = self.output_location + "/" + stdout_file
            f.write("#PBS -o %s\n" % self.ofile)
        else:
            self.ofile = None
            f.write("#PBS -o %s\n" % self.output_location)

        if stderr_file:
            self.efile = self.output_location + "/" + stderr_file
            f.write("#PBS -e %s\n" % self.efile)
        else:
            self.efile = None
            f.write("#PBS -e %s\n" % self.output_location)

        if join:
            f.write("#PBS -j oe\n")
            
        if address:    
            s="#PBS -M "+address+"\n"
            f.write(s)

        if queue:
            s="#PBS -q "+queue+"\n"
            f.write(s)

        if env:    
            f.write("#PBS -V\n")

        if mail:
            s="#PBS -m "+mail+"\n"
            f.write(s)

        if chdir:
            s="cd "+chdir+"\n"
            f.write(s)

        else:
            s="cd $PBS_O_WORKDIR\n";
            f.write(s);

        #if timing:
        #    self.timing = self.output_location + "/" + self.batch_file_name + ".timing"
        #    self.cmd = "/usr/bin/time -o " + self.timing + " -f \"%U\" " + self.cmd
        #else:
        #    self.timing = None

        if self.modules != None:
            self.cmd = "; ".join(["module load " + x for x in self.modules]) + "; " + self.cmd

        f.write(self.cmd)

        f.close()
        self.jobid=0
        self.split = self.suppress_pbs   # Set to true when the .e file gets split, unless we are not using pbs (in which case it is irrelevant)
        
        if not os.path.isfile(self.epilogue):
            open(self.epilogue, "w").write(epilogue_str)
            subprocess.call("chmod 500 %s" % (self.epilogue), shell=True)



### submitjob
### Parameters:
###   file is the job script file
###   if preserve is True, don't delete the job script. Delete otherwise
###   optionalFlag is the flag after qsub
###   retry (default set to retry 600 times), the number of times, the job will be submitted in retry
###   seconds between retry (default is 10 seconds)
### return job id if successful
### return -1 if not
    def submit(self, preserve=False, print_qsub = False, job_limit = job_limit, delay=10, user=current_user ):
        """Submit job to redhawk.  Optional parameters:
           * preserve: if False, delete the batch file.  Default = true.
           * job_limit: If the user currently has this many jobs on the batch, wait until one finishes.
           """
        if self.suppress_pbs:   # We are not using the pbs job handler -- direct call to Popen
            self.jobid = random.randint(1000000, 9999999)
            if not self.ofile:
                self.ofile = self.output_location + "/" + self.jobname + ".o" + str(self.jobid)

            if not self.efile:
                self.efile = self.output_location + "/" + self.jobname + ".e" + str(self.jobid)

            if not re.search("[^2]\>", self.cmd):
                self.cmd += " > {ofile}".format(ofile = self.ofile)
            if not re.search("2\>", self.cmd):
                self.cmd += " 2> {efile}".format(efile = self.efile)
            self.p = subprocess.Popen(self.cmd, shell=True)
            

        else:
            if job_limit > 0:
                limit_jobs(limit=job_limit, delay=delay, user=user, use_pbs = not self.suppress_pbs)

            optionalFlag= '-l epilogue=' + self.epilogue
            retry=600
            sleepTimeBetweenRetry=10
            trial=1;
            cmd = "qsub " + optionalFlag + " " + self.batch_file_name

            while (trial < retry): 
                p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # Sleep to ensure the prh file is created before termination
                (output, error) = [x.decode() for x in p.communicate()]
                if p.returncode == 0:
                    break
                trial = trial + 1

            if trial == retry:
                return -1

            t=re.split('\.',output)
            self.jobid=t[0]

        if not preserve:
            os.remove(self.batch_file_name)



        if not self.ofile:
            self.ofile = self.output_location + "/" + self.jobname + ".o" + str(self.jobid)

        if not self.efile:
            self.efile = self.output_location + "/" + self.jobname + ".e" + str(self.jobid)

        self.rfile = self.output_location + "/" + self.jobname + ".r" + str(self.jobid)

        if print_qsub:
            if self.suppress_pbs:
                print("Popen: jobid: {jobid}, process id: {pid}".format(jobid = self.jobid, pid = self.p.pid))
            else:
                print("PBS: jobid: {jobid}".format(jobid = self.jobid))

        self.status = "running"

        job_list.add(self)
        return self

    def submitjob(self, preserve=False, print_qsub = False, job_limit = job_limit, delay=10, user=current_user ):
        """Depricated: replaced with submit()"""
        return self.submit(preserve, print_qsub, job_limit, delay, user)

### isJobRunning
### This is primarily useful for waiting for a _submitted_ job to finish
###    return False if the job is done, completed for sure
###    return True if the job is in Q, R states [ or that PBS/Torque is not available ]
### Prereq is jobid must be a submitted job
    def isJobRunning(self, numTrials = 3, delay = 5):
        """Query of the object represented by the job is still running."""
        if self.status == "finished":
            return False

        if self.suppress_pbs:
            if not self.p.poll():
                return True
            else:
                self.status = "finished"
                job_list.remove(self)
                return False

        # If we are using pbs...
        while True:
            if self.ofile_exists() or self.efile_exists():  #output.find(magicString) >=0 or redhawkStatsRe.search(output):
                while not(self.efile_exists() and self.ofile_exists()):
                    pass
                self.status = "finished"
                job_list.remove(self)
                return False

            cmd = "qstat " + str(self.jobid)
            magicString = 'Unknown Job ID'
            (output, error) = [x.decode() for x in subprocess.Popen(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()]
            if redhawkInQueueRe.search(output):
                return True

            time.sleep(delay)

        raise PBSError("RedHawk error: out of queue, no output file.  OFILE: %s" % (self.ofile))

    
    def wait(self, delay=10):
        """Spin until job completes.
        """
        if self.suppress_pbs:
            self.p.wait()
            self.status = "finished"
        else:
            while self.isJobRunning() == True:
                time.sleep(delay)
        return self.ofile_exists()  

    def wait_on_job(self, delay=10):  
        """Depricated: replace with wait"""
        return self.wait(delay)

    def ofile_name(self):
        """Get the name of the file containing the job stdio output."""
        return self.ofile

    def efile_name(self):
        """Get the name of the file containing the job stderr output."""
        return self.efile

    def rfile_name(self):
        """Get the name of the file containing the rfile output."""
        return self.rfile

    #def timing_name(self):
    #    """Get the name of the timing file."""
    #    return self.timing

    def ofile_exists(self):
        """Does the file contiining the job stdout output exist?"""
        return os.path.isfile(self.ofile)

    def efile_exists(self):
        """Does the file contiining the job stderr output exist?"""
        return os.path.isfile(self.efile)

    def ofile_handle(self):
        """Return a handle to the file containing the job stdout output."""
        if not self.status == "finished":
            raise NameError("redhawk: unfinished ofile check")
        tries = 0
        while not self.ofile_exists() and tries < self.file_limit:
            time.sleep(self.file_delay)
            tries = tries+1
        
        if os.path.isfile(self.ofile_name()):
            return open(self.ofile_name(), "r")

        raise NameError("redhawk: unfound ofile")

    def efile_handle(self):
        """Return a handle to the file containing the job stderr output."""
        if not self.status == "finished":
            raise NameError("redhawk: unfinished efile check")

        tries = 0
        while not self.efile_exists() and tries < self.file_limit:
            time.sleep(self.file_delay)
            tries = tries+1
        
        if os.path.isfile(self.efile_name()):
            return open(self.efile_name(), "r")

        raise NameError("redhawk: unfinished efile check")

    def rfile_handle(self):
        """Return a handle to the file containing the resource description."""
        if self.suppress_pbs:
            raise PBSError("Cannot check resource file on a pbs-suppressed job")
        self.split_efile()
        return open(self.rfile)

    def ofile_string(self):
        """Return the entire contents of the stdout file as a single string."""
        fp = self.ofile_handle()
        if (fp):
            return "\n".join([line.rstrip() for line in fp]) + '\n'
        return None

    def efile_string(self):
        """Return the entire contents of the stderr file as a single string."""
        fp = self.efile_handle()
        if (fp):
            return "\n".join([line.rstrip() for line in fp]) + '\n'
        return None

    # def get_timing(self, delete_timing = True):
    #     """Get the time-generated user runtime and delete the timing file.
    #     (Assumed to be the last line of the timing file.)
    #     By default, erases the timing file."""
    #     if not self.timing:
    #         return None
    #     else:
    #         try:
    #             with open(self.timing) as fp:
    #                 lines = fp.readlines()
    #             if delete_timing:
    #                 os.remove(self.timing)
    #             return float(lines[-1])
    #         except:
    #             sys.stderr.write("Redhawk.runtime: invalid runtime (%s)\n" % (self.timing))
    #             sys.exit(1)
    #             return None
                                
    def erase_files(self, empty_only = False):
        """Erase the stdio and stderr files."""

        try:
            if not empty_only or os.path.getsize(self.ofile_name()) == 0:
                os.remove(self.ofile_name())
        except:
            pass

        try:
            if not empty_only or os.path.getsize(self.efile_name()) == 0:
                os.remove(self.efile_name())
        except:
            pass

        try:
            if not empty_only or path.getsize(self.rfile_name()) == 0:
                os.remove(self.rfile_name())
        except:
            pass

        return None
    
    def get_results(self, resources = False, cleanup=True):
        """Retrieve strings and cleanup."""
        if resources and self.suppress_job:
            raise PBSError("Cannot get resources on a pbs-suppressed job")
        self.wait()
        self.split_efile()
        stdout_str = self.ofile_string()
        stderr_str = self.efile_string()
        if resources:
            T = self.getResources()
        if cleanup:
            self.erase_files()
                
        return (stdout_str, stderr_str, T) if resources else (stdout_str, stderr_str)



    def getResults(self, cleanup=True):
        """Legacy"""
        return self.get_results(cleanup)
	

    def loadResources(self):
        if not self.resources:
            self.wait()

            if self.suppress_pbs:
                self.resources = [-1, -1, -1, -1]
            else:
                try:
                    fp = self.rfile_handle()
                except:
                    self.resources = (-1, -1, -1, -1)
                    return None

                for line in fp:
                    if line.startswith("Resources Used:"):
                        r = re.search("cput=(\d\d):(\d\d):(\d\d),mem=(\d+)kb,vmem=(\d+)kb,walltime=(\d\d):(\d\d):(\d\d)", line)
                        if not r:
                            raise PBSError("Bad resource line: " + line)
                        cpu_time = 60*int(r.group(1)) + 3600*int(r.group(2)) + int(r.group(3))
                        wall_time = 60*int(r.group(6)) + 3600*int(r.group(7)) + int(r.group(8))
                        memory = 1024*int(r.group(4))
                        vmemory = 1024*int(r.group(5))
                        self.resources = (cpu_time, wall_time, memory, vmemory);  
                        break
                



    def getResources(self, cleanup=True):
        """Return cpu_time, wall_time, memory, and virtual memory used"""
        self.loadResources()
        return self.resources

    def cpu_time(self):
        if self.suppress_pbs:
            raise PBSError("Cannot get resources on a pbs-suppressed job")
        self.loadResources()
        return self.resources[0]

    def memory(self):
        if self.suppress_pbs:
            raise PBSError("Cannot get resources on a pbs-suppressed job")
        self.loadResources()
        return self.resources[2]

    def vmemory(self):
        if self.suppress_pbs:
            raise PBSError("Cannot get resources on a pbs-suppressed job")
        self.loadResources()
        return self.resources[3]

    def wait_on_job_limit(self, limit=200, delay=10, user=current_user):
        """Depricated: use the stand-alone function job_limit."""
        sys.stderr.write("pbsJobHandler: wait_on_job_limit method is depricated; please switch to the stand-alone job_limit function.")
        limit_jobs(self, limit, delay, user = user, use_pbs = not self.suppress_pbs)

    def split_efile(self):
        """Split the .e<id> file into a .e<id> and .r<id> file"""
        self.has_split = True
        if not self.split:
            if not self.efile_exists():
                raise PBSError("Call to split_efile() when efile doesn't exist (%s, %s)" % (self.cmd, self.efile))
            self.split = True 

            with open(self.efile) as fp: line = "".join([line for line in fp])
            count = 0
            r = None
            while r is None:
                r = re.search("(.*)Redhawk Epilogue Args:\s*(.+)", line, re.DOTALL)
                if not r:
                    count += 1
                    if count == 1:
                        break
                    time.sleep(5)

            if r:
                first, second = r.group(1,2)
                with open(self.efile, "w") as fp: fp.write(first)
                with open(self.rfile, "w") as fp: fp.write(second)
            else:
                with open(self.rfile, "w") as fp: fp.write(bad_ep_str + "\n")


### Hold until the user has < limit jobs in the circulation
def limit_jobs(limit=job_limit, delay=10, use_pbs=pbs_present, user=current_user):
    """Spin until the user has less < limit jobs in circulation.
        limit = 0 signals no limit."""
    if limit == 0:
        return None
    while 1==1:
        numJobsInQueue = get_number_of_jobs_in_queue(current_user, use_pbs)

        if (numJobsInQueue < limit):
            return None
        time.sleep(delay)


### return the number of jobs in queue, whatever the state is
def get_number_of_jobs_in_queue(user=current_user, using_pbs = pbs_present):
    """Get the number of user jobs currently sitting in the queue or running.  Not yet designed to work with pbs suppression."""
    global job_list
    if not using_pbs:
        job_list = [j for j in job_list if j.isJobRunning()]            
        return len(job_list)

    cmd = "qstat -u "+user + " 2>/dev/null | grep " + user
    output,error = [x.decode() for x in subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()]
    return len([True for line in output.split("\n") if line and not redhawkStatsRe.search(line)])
                                                        
         

def storePBS(pbsList, fp):
    """Pickle and save a PBS job to a specified file pointer."""
    pickle.dump(pbsList, fp)

def loadPBS(fp):
    """Recover a list of pickled jobs from a specified file pointer."""
    return pickle.load(fp)


def relaunch(args = sys.argv, force = False, walltime = "40:00:00", python = "python"):
    """Detects whether program is being run on the head node.  If so, relaunch identical program on a compute node and quit."""
    o, e = [x.decode() for x in subprocess.Popen(["hostname"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()]
    if force or re.match("mualhpcp", o):
        o = pbsJobHandler(batch_file = "relaunch" + str(os.getpid()),  executable = python + " " + " ".join(args), walltime=walltime)
        o, e = o.submit().getResults()
        sys.stdout.write("STDOUT: " + o)
        sys.stderr.write("STDERR: " + e)
        return True
    return False
        
############## Allow for a direct launch of a program
if __name__ == "__main__":
    parser = ArgumentParser(description='Launch a redhawk job.')

    parser.add_argument('command', action = "store", type = str, nargs = '+')

    term = parser.add_argument_group("Input/Output switches")
    term.add_argument('--wait', action = "store_true", dest = "wait", help="wait on job completion", default = False)
    term.add_argument('-o', '--output', action = "store", dest = "output", help = "output file (stdout by default)", default = None)
    term.add_argument('-e', '--error', action = "store", dest = "error", help = "error file (stderr by default)", default = None)
    term.add_argument('-r', '--resources', action = "store", dest = "resources", help = "resource file (same as output by default)", default = None)
    term.add_argument('-S', '--suppress_output', action = "store_true", dest = "suppress", help="Quit without output", default = False)

    settings = parser.add_argument_group("Main job-related settings")
    settings.add_argument('-c', '--create', action = "store_true", dest = "create", help = "Create the batch file and quit.", default = False)
    settings.add_argument('-n', '--nodes', action = "store", type = int, dest = "nodes", help = "Number of nodes", default = 1)
    settings.add_argument('-p', '--ppn', action = "store", type = int, dest = "ppn", help = "Number of processors per node", default = 1)
    settings.add_argument('--big_mem', '--mem', action = "store_true", dest = "mem", help = "Use 128Gb machine", default = False)
    settings.add_argument('-w', '--walltime', action = "store", type = str, dest = "walltime", help = "Reserved walltime", default = "10:00:00")
    settings.add_argument('-m', '--modules', action = "store", type = str, nargs = "+", dest = "RHmodules", help = "required redhawk modules", default = None)
    settings.add_argument('-O', '--output_location', action = "store", type = str, dest = "output_location", help = "Output location", default = None)
    settings.add_argument('-d', '--dir', action = "store", type = str, dest = "target_directory", help = "target directory", default = None)

    settings2 = parser.add_argument_group("Less important job-related settings")
    settings2.add_argument('-b', '--batch', action = "store", type = str, dest = "batch", help="Batch file name", default = "redhawk_run")
    settings2.add_argument('-P', '--pid_off', action = "store_false", dest = "pid", help = "Suppress use of pid in file names", default = True)
    #settings2.add_argument('-T', '--time_off', action = "store_false", dest = "timing", help = "Suppress runtime reporting", default = True)
    settings2.add_argument('-R', '--resources_off', action = "store_false", dest = "print_resources", help = "Suppress resource usage reporting", default = True)
    settings2.add_argument('-K', '--keep_files', action = "store_true", dest = "keep", help = "Keep files generated", default = False)
    settings2.add_argument('-j', '--job_name', dest = "job_name", help = "Redhawk job name", default = None)

    #term = parser.add_argument_group("Alternative jobs (internal option -- not intended for users)")
    #term.add_argument('--post_process', action = 'store_true', dest = 'post_process', help="Process the result", default = False)

    args = parser.parse_args()

    #if args.post_process:
    #    post_process(args.command[0])
    #    sys.exit(0)


    
    p = pbsJobHandler(batch_file = args.batch, executable = " ".join(args.command), job_name = args.job_name, use_pid = args.pid, nodes = args.nodes, ppn = args.ppn, mem = args.mem, 
                      walltime = args.walltime, output_location = args.output_location, chdir = args.target_directory, RHmodules = args.RHmodules)
    if args.create:
        sys.exit(0)
    o = p.submit(preserve = args.keep)

    if args.suppress:
        sys.exit(0)


    o.wait()
    ofp = sys.stdout if not args.output or args.output == '-' else open(args.output, "w")
    efp = sys.stderr if not args.error or args.error == '-' else open(args.error, "w")
    rfp = ofp if not args.resources else (sys.stdout if args.resources == '-' else open(args.resources, "w"))
    out, err = o.get_results(cleanup = not args.keep)
    ofp.write(out)
    ofp.write(err)

    #t = o.get_timing(delete_timing = not args.keep)
    #if args.timing:
    #    rfp.write("Time: " + str(t) + "\n")
    if args.print_resources:
        A = o.getResources()
        rfp.write("CPU Time: %d\n" % (A[0]))
        rfp.write("Wall time: %d\n" % (A[1]))
        rfp.write("Memory: %d\n" % (A[2]))
        rfp.write("Vmemory: %d\n" % (A[3]))
        
        



################################
### Sample code: Running a single job
### exe = "ls *"     # The command we want to run (to run multiple commands, seperate with semi-colons)
### o = pbsJobHandler(batch_file = "batch.txt", executable = exe, mail = "bea");    # Set up the job
### o.submit()                           # Submit the job to a redhawk queue
### if o.isJobRunning(): print "yes"     # Check to see if job is still in the queue or running
### o.wait()                             # "spin" until job is finished
### output_string = o.ofile_string()     # Get the ouput
### o.erase_files()                      # Erase the output files




