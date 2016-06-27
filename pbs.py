from sets import Set
import os
import sys
import time
import imp

class PBS(object):
    
    def __init__(self,job_name,path_local,commands,tm_start,log_func,wall_time,verbose=False):
        
        self.util = imp.load_source('util','util.py')
        
        self.env = imp.load_source('env','env.py')
        
        #name of this job
        self.job_name=job_name #7550
        self.path_master = os.path.join(self.util.CLUSTER_JOB_DIR,path_local)
        #/scratch/local/N on cluster node
        self.path_node = os.path.join(self.util.CLUSTER_NODE_DIR,path_local)
        self.commands  = commands
        self.wall_time = wall_time
        #for log only
        self.log_func = log_func
        self.log_queue = False
        self.log_running = False
        
        self.file_current_cmd = os.path.join(self.path_master,'current.txt')
        self.util.shell_exec('touch %s' % self.file_current_cmd)
        #self.file_node = os.path.join(self.path_master,'node.txt')
        #verbose mode        
        self.verbose = verbose
        #self.job_status = None
        self.pbs_job_id = 0
        self.file_pbs = None
        self.file_error  = None
        self.file_output = None
        self.tm_start = tm_start
        self.status_code={'Q':1, #in the queue
                          'R':2,#running
                          'E':3, #exiting
                          'H':4, #held
                          'S':5, #suspended
                          'T':6, #relocated
                          'W':7, #Waiting
                          'C':8 #completed
                          }
        
        self.required_data_files = self.get_required_data_files()
        
    def run(self):
        """create PBS script and submit to run on cluster nodes
        """
        last_modified_time = os.path.getmtime(self.file_current_cmd)
        file_pbs = self.make_pbs(self.commands)
        msg = self.util.shell_exec('qsub %s' % file_pbs)
        #get the job id returned by qsub
        self.pbs_job_id = str(msg[0].split('.')[0])
        #Wait until job completed or failed.
        completed = False
        
        while not completed:
            
            time.sleep(self.util.PBS_QUERY_INTERVAL)
            
            scode = self.get_status()
            
            if scode==1: #in the queue
                if not self.log_queue:#on/off tag to prevent write more than one status info to the log file.
                    self.log_queue = True
                    self.log_func('Job is in queue')
            
            if scode==2: #running
                if not self.log_running:#on/off tag to prevent write more than one status info to the log file.
                    self.log_running = True
                    self.log_func('Job is loading')
            
            if scode==0 or scode==-1 or scode==None:
                """None ---> The job finished and its ID was removed from the Q before the interval querying.
                    9 ---> A "tmp.done" is found which means the job is successfully finished.
                    -1 ---> A "tmp.fail" is found which means the job is failed.
                """
                completed = True
            
            current_modified_time = os.path.getmtime(self.file_current_cmd)
            if ((current_modified_time - last_modified_time)>0.1):
                last_modified_time = current_modified_time
                sz_current_cmd = file(self.file_current_cmd).read().strip()
                if sz_current_cmd:
                    self.log_func(sz_current_cmd)
    
    def before(self):
        """executed before user's commands"""        
        cmds = []        

        #add all subdirectories to PATH
        #cmd = 'find %s -type d' % self.util.CLUSTER_APP_DIR
        #subdirs = self.util.shell_exec(cmd)
        #if subdirs:
        #    cmds.append('export PATH=$PATH:%s' % ';'.join(subdirs))

        #do not report error if this command fails (sometimes it failed because the file permission) 
        cmds.append('rm -fr %s/* 2>/dev/null || true' % self.util.CLUSTER_NODE_DIR)
        #cmds.append('set -o errtrace')
        
        #make sure all files created in this job can be deleted by other users
        cmds.append('umask 000')
        
        #make  /scratch/local/JOB_NAME
        cmds.append('mkdir -p %s' % self.path_node)
        
        #change to the local disk on cluster node
        cmds.append('cd %s' % self.path_node)
        
        #copy source files from interactive node to cluster node
        #cmds.append('cp -r %s/* %s' % (self.path_master,self.path_node))
        cmds.append('cp -r %s/* .' % self.path_master)
        
        #cmds.append('cp -ur %s/*  %s/' % (CLUSTER_DATA_DIR,os.path.join(CLUSTER_NODE_DIR,self.path_local)))
        
        #only copy required library files
        for f in self.required_data_files:
            #cmds.append('cp -r %s/%s %s' % (self.util.CLUSTER_DATA_DIR,f,self.path_node))
            cmds.append('cp -r %s/%s .' % (self.util.CLUSTER_DATA_DIR,f))
        

        return cmds
        
    def after(self):
        """executed after user's commands"""
        cmds = []
        
        #delete library files
        if self.required_data_files:
            cmds.append('rm -fr %s' % ' '.join(self.required_data_files))

        #copy (update) result files to cluster master node
        cmds.append('cp -ur %s/* %s/' % (self.path_node,self.path_master))

        #delete job folder on cluster node
        cmds.append('rm -fr %s' % self.path_node)
        
        #change working directory on interactive node
        cmds.append('cd %s' % self.path_master)
        
        #make a tag file that indicates the whole processing is done, ready to upload
        cmds.append('touch %s' % 'tmp.done')
        
        return cmds

    def make_pbs(self,commands):
        """create script for PBS"""
        
        #echo 'Failed:' $_ >> stderr.txt
        if self.required_data_files:
            fail_clean = """failclean()
{    
    rm -f %s
    cp -ur %s/* %s/
    rm -fr %s
    cd %s
    touch tmp.fail
    exit 1
}
trap 'failclean' ERR TERM""" % (' '.join(self.required_data_files),self.path_node,self.path_master,self.path_node,self.path_master)
            
        else:
            fail_clean = """failclean()
{
    cp -ur %s/* %s/
    rm -fr %s
    cd %s
    touch tmp.fail
    exit 1
}
trap 'failclean' ERR TERM""" % (self.path_node,self.path_master,self.path_node,self.path_master)
            

        pbs_command = []
        pbs_command.append('#!/bin/bash')
        pbs_command.append('#PBS -A %s' % 'kaplan-em')
        pbs_command.append('#PBS -l nodes=1:ppn=12,walltime=%s' % self.estimate_walltime())
        pbs_command.append('#PBS -N %s' % self.job_name)
        pbs_command.append('#PBS -o %s' % os.path.join(self.path_master,'stdout.txt'))
        pbs_command.append('#PBS -e %s' % os.path.join(self.path_master,'stderr.txt'))
        #pbs_command.append('#PBS -p %d' % self.job_priority)
        
        pbs_command.append(fail_clean)
        #see which node
        pbs_command.append('\n'.join(self.before()))
        #pbs_command.append("uname -a | awk '{print $2}' > %s" % self.file_node)

        #pbs_command.append('\n'.join(self.commands))
        tmp_cmds = []
        for i in range(0,len(self.commands)):
            tmp_cmds.append('echo "Running: %d of %d" > %s' % (i+1,len(self.commands),self.file_current_cmd))
            tmp_cmds.append(self.commands[i])
        pbs_command.append('\n'.join(tmp_cmds))
        
        pbs_command.append('\n'.join(self.after()))
        file_pbs = os.path.join(self.path_master,'tmp.pbs') 
        fh=file(file_pbs,'w')
        fh.write('\n'.join(pbs_command))
        fh.close()
        return file_pbs

    def estimate_walltime(self):
        #estimated_time = (24,0,0) #10 days
        #return '%02d:%02d:%02d' % estimated_time
        #return '%02d:%02d:%02d' % (MAX_WALLTIME_HOURS,0,0)
        return '%02d:%02d:%02d' % (self.wall_time,0,0)

    def get_status(self):
        """
        Status value    Meaning
            E    Exiting after having run
            H    Held
            Q    Queued
            R    Running
            S    Suspended
            T    Being moved to new location
            W    Waiting for its execution time
            C    Recently completed (within the last 5 minutes)
        """
        fdone = os.path.join(self.path_master,'tmp.done')
        ffail = os.path.join(self.path_master,'tmp.fail')
        
        if os.path.exists(fdone):
            return 0
        
        elif os.path.exists(ffail):
            return -1
        
        else:
            cmd = "qstat -u %s | awk '{print $1, $10}' | grep %s | awk '{print $2}'" % (self.util.CLUSTER_USER,str(self.pbs_job_id))
            #cmd = "qstat -u %s | grep %s. | awk '{print $10}'" % (self.util.CLUSTER_USER,str(self.pbs_job_id))
            try:
                jobstatus = self.util.shell_exec(cmd)[0]
                #return self.status_code.get(status,9)
            except Exception,e:
                #print self.util.shell_exec('qstat -u %s' % self.util.CLUSTER_USER)
                print 'Error in JOB QUERY:',e,self.pbs_job_id
                #the job id is not available, that means the job is finished and removed from the queue. 
                return None
            #if none of above code is present, then return "9" which is undefined.
            return self.status_code.get(jobstatus,9)

    def get_required_data_files(self):
        """A command may required some library files, e.g. hg19.nix(Genome Index for HG19)
        ,hg19.fasta, etc.
        
        hg19.bwa.illumina.xxxx
        """
        ret = []
        available_data_files = self.env.DataLibraryFile().get_data_files()
        m = {}.fromkeys(available_data_files)
        for c in self.commands:
            ps = c.split()
            for p in ps:
                if p.endswith('.bwa.illumina') or p.endswith('.bwa.solid'):
                    ret.append('%s.*' % p)
                if m.has_key(p):
                    ret.append(p)
                    if p.endswith(".fasta"): #A.fasta
                        ret.append(p+".fai") #A.fasta.fai
                        ret.append(p[:-len(".fasta")]+".dict") #A.dict
                    if p.endswith(".vcf"): #A.fasta
                        ret.append(p+".idx") #A.fasta.fai
        return list(Set(ret))

