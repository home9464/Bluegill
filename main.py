#!/usr/bin/python
import multiprocessing
import Queue
import time
import os
import sys
import imp

"""main entry
"""

class Main(object):
    def __init__(self):
        self.all_jobs = {}
        self.job_index = 0
        self.util = imp.load_source('util','util.py')
        if self.util.TEST:
            print 'Test mode'
        self.job = imp.load_source('job','job.py')
        self.user_jobs = multiprocessing.Manager().dict()
        self.pending_jobs = Queue.PriorityQueue()
        
    def run(self):
        self.util.make_directories()
        self.util.clean_job_directory()
        print "[%s] Synchronizing ..." % self.util.now(),
        self._sync()
        print '[OK]'
        while True:
            self.util = imp.load_source('util','util.py')
            self.job = imp.load_source('job','job.py')
            self.load()
            time.sleep(self.util.RSYNC_QUERY_INTERVAL)
            self.update_alive_jobs()
            self._sync()
    
    def load(self):
        #temporarily stop receiving new jobs
        #return job_index,alive_jobs
        #job_idx = job_index
        #job_dict = alive_jobs
        
        #clean every time?
        new_jobs,abort_jobs = self.util.sync_query()
        for job_name, create_time in abort_jobs:
            
            owner = self.util.get_job_owner(job_name)
            try:
                pa = self.all_jobs.get(job_name,None)
                if pa:
                    pid = pa[0].pid
                    ji = pa[1]
                    #p.terminate()
                    self.util.shell_exec('kill -9 %s' % pid) #SIGKILL is better than SIGTERM?
                    ps = multiprocessing.Process(target=self.abort,args=(owner,job_name,ji,self.user_jobs))
                    ps.start()
                    #ps.join()
                    #lastly remove this job on Master node
                    #self.util.shell_exec('rm -fr %s' % os.path.join(self.util.CLUSTER_JOB_DIR,str(ji)))
                        
            except Exception,e:
                print e
        
        
        tmp_job_dict = {}  
        tmp_jobs = []

        for job_name,create_time in new_jobs:
            #job name is the input path.
            if not self.all_jobs.has_key(job_name):  #can not accept same job twice
                
                #print self.util.extract_email_from_cmd(job_name)
                job_owner = self.util.get_job_owner(job_name)
                #if self.util.WHITE_BOX:
                #    if not job_owner in self.util.WHITE_BOX: #only accept jobs from users in white box
                #        continue
                    
                #if self.util.BLACK_BOX:
                #    if '*' in self.util.BLACK_BOX or job_owner in self.util.BLACK_BOX:#block any jobs from users in black box
                #        continue
                priority,max_num_jobs,wallhours = self.util.get_user_priority(job_owner)
                tmp_job_dict.setdefault(job_owner,[]).append(job_name)
                myjobs =  self.user_jobs.get(job_owner) #how many jobs are running under this user?
                if myjobs: 
                    _total = len(myjobs) + len(tmp_job_dict[job_owner])
                    if _total < max_num_jobs:
                        tmp_jobs.append((job_name,job_owner,priority,wallhours))
                            
                else:#no other running jobs
                    _total = len(tmp_job_dict[job_owner])
                    if _total < max_num_jobs:
                        tmp_jobs.append((job_name,job_owner,priority,wallhours))
                
                """                    
                job_priority = self.util.get_job_priority(job_owner)
                job_walltime = self.util.get_job_walltime(job_priority)
                if not self.util.PRIORITY_USER.has_key(job_owner):#not a priority user
                    tmp_job_dict.setdefault(job_owner,[]).append(job_name)
                    myjobs =  self.user_jobs.get(job_owner) #how many jobs are running under this user?
                    if myjobs: 
                        _total = len(myjobs) + len(tmp_job_dict[job_owner])
                        if _total <= self.util.NUM_JOBS_MAX:
                            tmp_jobs.append((job_name,job_owner,job_priority,job_walltime))
                            
                    else:#no other running jobs
                        _total = len(tmp_job_dict[job_owner])
                        if _total <= self.util.NUM_JOBS_MAX:
                            tmp_jobs.append((job_name,job_owner,job_priority,job_walltime))
                else:
                    tmp_jobs.append((job_name,job_owner,job_priority,job_walltime))
                """
                
        available_nodes = self.util.NUM_CLUSTER_NODES_AVAILABLE - len(self.all_jobs)
        _job_rank = 0

        for t in tmp_jobs:
            self.pending_jobs.put(t) #push it back
            
        tmp_jobs2 = []
        while not self.pending_jobs.empty():
            job_name,job_owner,job_priority,job_walltime = self.pending_jobs.get()#pop a pending job from the queue
            _job_rank += 1
            #if _job_rank > available_nodes:
            #print 'Rank:', job_name,_job_rank
            #    my_rank = _job_rank-available_nodes 
            #    self.util.accept_pending_job(job_name,my_rank)
            self.util.accept_pending_job(job_name,_job_rank)
            
            tmp_jobs2.append((job_name,job_owner,job_priority,job_walltime)) #push it back

        for t in tmp_jobs2:
            self.pending_jobs.put(t) #push it back
            
        for k in range(0,available_nodes):
            if not self.pending_jobs.empty():
                job_name,job_owner,job_priority,job_walltime = self.pending_jobs.get()#pop a pending job from the queue
                #print 'Run:', job_name
                myjobs =  self.user_jobs.get(job_owner,[])#check out other jobs belong to this user
                """
                if myjobs: #this user has some job in process
                    _priority,_max_num_jobs,_wallhours = self.util.get_user_priority(job_owner)
                    if len(myjobs) >= _max_num_jobs: #and submitted more than 3 jobs
                        job_walltime = self.util.HOURS_LIMITED_WALLTIME
                    #if not self.util.PRIORITY_USER.has_key(job_owner):#this user is not a priority user
                    #    if len(myjobs) >= self.util.NUM_JOBS_UNLIMITED_WALLTIME: #and submitted more than 3 jobs
                    #        job_walltime = self.util.HOURS_LIMITED_WALLTIME
                else: #first job of this job
                    myjobs = []
                """
                self.job_index += 1
                myjobs.append(self.job_index)
                self.user_jobs[job_owner] = myjobs
                #print job_idx,job_owner,job_name,job_walltime
                p = multiprocessing.Process(target=self.begin,args=(job_owner,job_name,self.job_index,job_walltime,self.user_jobs))
                p.start()
                self.all_jobs[job_name] = (p,self.job_index)
    
    def _sync(self):
        self.util.sync_get_app()
        self.util.sync_get_data()
        self.util.sync_get_pipeline()

    def update_alive_jobs(self):
        #eliminate zombie <defunct> child process once it is done 
        ret = {}
        for k,v in self.all_jobs.items():
            if not v[0].is_alive():
                #print '%s is dead' % k
                v[0].terminate()
            else:
                ret[k] = v
        self.all_jobs = ret
        
    def begin(self,owner,name,index,walltime,jobdict):
        #start a new job
        self.job.Job(owner,name,index,walltime,jobdict).start()
    
    
    def abort(self,owner,name,index,jobdict):
        #abort a running job
        self.job.Job(owner,name,index,0,jobdict).abort()
        
if __name__=='__main__':
    Main().run()

