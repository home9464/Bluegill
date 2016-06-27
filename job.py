import re
import sys
import os
import time
import imp
import random

from version import VersionConf

           
class Job:
    
    def __init__(self,job_owner,job_name,job_index,job_walltime,job_dict):

        self.util = imp.load_source('util','util.py')
        self.env = imp.load_source('env','env.py')
        self.gmail = imp.load_source('gmail','gmail.py')
        self.pipelinebuilder = imp.load_source('pipeline','pipeline.py')
        self.pbs = imp.load_source('pbs','pbs.py')
        
        self.tm_start = time.time()
        self.job_name = job_name
        self.job_index = job_index
        self.path_input = job_name #the input path D (not full, JOB_PATH/D) name
        self.path_output = None #by default equals to self.path_input 
        #self.path_local = filter(None,job_path.split(os.sep))[0]+'_'+str(int(time.time()))
        self.path_local = str(job_index)
        self.user_email = None
        self.description = []
        self.user_lab=None
        self.genome_build = 'Unknown'
        self.gnomex_analysis_dir = None
        self.gnomex_analysis_type = None
        self.gnomex_analysis_id = None
        
        #self.wall_time = self.util.MAX_WALLTIME_HOURS
        self.job_owner=job_owner
        self.wall_time = job_walltime
        self.job_dict = job_dict
        #direct output to existing GNomEx Analysis
        #self._analysis_id = None
        self.how_to_send_email = []
        
        self.version = None
        
    def clean_job_dict(self):
        ind = -1
        try:
            myjobs = self.job_dict[self.job_owner]
            for i,j in enumerate(myjobs):
                if j==self.job_index:
                    ind = i
            if not ind == -1:
                myjobs.pop(ind)
                self.job_dict[self.job_owner] = myjobs
        except:
            print 'Error:',self.job_dict
        
    def parse_cmd(self):
        """parse the cmd.txt file from a job folder.
        
        @return [] commands to be executed
        
        @bam,sort
        
        """
        #download the "cmd.txt"
        self.util.sync_get_cmd(self.path_input,self.path_local)
        
        file_command = os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,self.util.FILE_CMD)
        
        total_commands = []
        
        analysis_id = None
        
        command_multiple_lines = []
        
        for line in file(file_command): #parse the "cmd" file to get instructions
            
            line = line.strip()
            if line:
                if line.startswith('#e '):# email
                    pattern_email = re.compile(r"(?:^|\s)[-a-z0-9_.]+@(?:[-a-z0-9]+\.)+[a-z]{2,6}(?:\s|$)",re.IGNORECASE)
                    matched_mail = pattern_email.findall(line.lstrip('#e').strip())
                    if matched_mail and not self.user_email: #can be multiple email-recipients
                        self.user_email = self.gmail.Gmail(matched_mail,self.job_name)
                        
                    #-a  mail will be sent when possible(equal -bcef).
                    #-b  mail will be sent when the job begins execution.
                    #-c  mail will be sent when the job completed successfully.
                    #-e  mail will be sent when the job has exception.
                    #-f  mail will be sent when the job failed.
                    #-n  mail will NOT be sent  (reverse of a)
                    
                    ex = []
                    toks = line.split(' ')
                    for t in toks:
                        t = t.strip()
                        if t.startswith('-'):
                            ex.extend([tt for tt in t[1:]])
                            
                    if not ex:
                        self.how_to_send_email = ['a']
                        
                    else:
                        if 'n' in ex:
                            self.how_to_send_email = []
                        else:
                            self.how_to_send_email = ex
                    
                elif line.startswith('#l '):# lab,  
                    self.user_lab = line.lstrip('#l').strip() 

                elif line.startswith('#a '):# direct output to existing Analysis folder, like #a A401  
                    self.gnomex_analysis_id = line.lstrip('#a').strip()
                    
                elif line.startswith('#o '):# direct output to existing Analysis folder, like #a A401  
                    self.path_output = line.lstrip('#o').strip()
                     
                elif line.startswith('#t '):# 12 hours
                    try:
                        self.wall_time = int(line.lstrip('#t').strip())
                    except:
                        self.wall_time = None
                        
                elif line.startswith('#v'):# description  
                    try:
                        self.version = line.lstrip('#v').strip()
                    except:
                        self.version = None
                        
                    
                elif line.startswith('##'):# description  
                    self.description.append(line.lstrip('##').strip())
                    
                elif line.startswith('#'):# description  
                    self.description.append(line)
                
                else:
                    if line.endswith('\\'):
                        command_multiple_lines.append(line.rstrip('\\').strip())
                        
                    else: #either a single line command or the last line of a multiple line command
                        cx = None
                        if command_multiple_lines:
                            command_multiple_lines.append(line)
                            #cx = ' '.join(command_multiple_lines)
                            cx = ';'.join(command_multiple_lines)
                            command_multiple_lines = []
                        else:
                            cx = line
                            
                        if cx.startswith('@'):
                            #command = filter(None, [s.strip() for s in cx.lstrip('@').split(' ')])
                            command = filter(None, [s.strip() for s in cx.split(' ')])
                            if self.user_lab:
                                for i in range(0,len(command)-1):
                                    if command[i] == '-g':
                                        self.genome_build = self.util.get_genomebuild(command[i+1])
                            #meta_command = command[0]
                            self.gnomex_analysis_type = 'Alignment'
                            #total_commands.extend(self.pipelinebuilder.PipelineBuilder(self.path_input,command,self.version).get_commands())
                            total_commands.append(' '.join(command))
                        else:
                            total_commands.append(cx)
                            
        vconf = VersionConf(self.version)
        mapped_commands = []
        for c in total_commands:
            if c.startswith('@'):
                mapped_commands.extend(self.pipelinebuilder.PipelineBuilder(self.path_input,c,self.version).get_commands())
            else:
                #print '#',self.env.Environment(vconf.get_app_path()).map_command(c)
                mapped_commands.append(self.env.AppLocator(vconf.get_app_path()).locate(c))
        self.version = vconf.get_current_version()
        #print mapped_commands
        
        if not self.user_email:
            err = 'Error: email address was not found. Please use #e MY_EMAIL'
            self.util.shell_exec_remote('echo "%s" > log.txt' % err, job_name=self.job_name, shell=True)
            raise Exception('No Email address provided')

        #if #lab is specified, then redirec the result to GNomEx:
        if not self.path_output:
            if self.user_lab:#create new GNomEx Analysis Folder
                remote_output_path = self.get_output_path()
                self.gnomex_analysis_dir = remote_output_path.split('/')[-1]
                self.util.make_remote_directory(remote_output_path)
                self.util.test_writeable_directory(remote_output_path)
                self.path_output = self.util.LOCAL_SERVER_PATH + remote_output_path.rstrip(os.sep) + os.sep
                if self.user_email:
                    folder_name =  os.path.basename(self.path_output.split(':')[-1].rstrip('/'))
                    #print 'GNomEx Analysis %s was created for %s' % (folder_name,self.job_name)
                    self.log_status('Job %s accepted, GNomEx Analysis %s created' % (self.job_name,folder_name))
                    if 'a' in self.how_to_send_email or 'b' in self.how_to_send_email:
                        self.user_email.send(subject='Job %s accepted, GNomEx Analysis %s created' % (self.job_name,folder_name))
                        #os.path.join(self.util.LOCAL_SERVER_JOB_PATH,self.path_input,folder_name)
                        self.util.shell_exec_remote('touch %s' % folder_name,self.path_input)
            else:
                self.path_output = os.path.join(self.util.LOCAL_SERVER_JOB_PATH,self.path_input).rstrip(os.sep)+os.sep
        else:
            self.path_output = self.util.LOCAL_SERVER_PATH+self.path_output.rstrip(os.sep)+os.sep
         
        #print self.path_output
        return mapped_commands
    
    def start(self):
        print '[%s] <<<<< %s (%d,%s)' % (self.util.now(),self.job_name,self.wall_time,str(self.job_dict[self.job_owner]))

        DRY_RUN =False
        #self.util.shell_exec_remote('[ -f %s ] && echo "1" || echo "0"' % self.util.FILE_TAG_RUNNING,self.path_input)
        msg = self.util.shell_exec_remote('stat -c %%s %s' % self.util.FILE_TAG_DRYRUN,self.path_input)
        if msg[0]=='0':
            DRY_RUN = True

        #remove the b file at first
        self.util.shell_exec_remote('rm -f R* %s %s %s %s %s' % (self.util.FILE_TAG_BEGIN,self.util.FILE_TAG_DRYRUN,self.util.FILE_LOG,self.util.FILE_STDOUT,self.util.FILE_STDERR), self.path_input)
        #create the r file        
        self.util.shell_exec_remote('touch %s' % self.util.FILE_TAG_RUNNING,self.path_input)
        
        self.util.shell_exec('mkdir -p %s' % os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local))
        
        successful = False
        
        try:
            total_commands = self.parse_cmd()
            
            commands = self.env.AppRuntimeMapper().map(total_commands)

            #self.user_email.send(subject = 'HELLO', body_msg='WORLD', body_file = None,additional_subject=None)
            
            if DRY_RUN:
                #print 'DRY RUN'
                #print '\n'.join(commands)
                self.user_email.send(subject = 'Dry run of Job %s completed' % self.job_name, body_msg='\n'.join(commands), body_file = None,additional_subject=None)
                sys.exit(0)
                
            if self.user_email:
                if 'a' in self.how_to_send_email or 'b' in self.how_to_send_email:
                    self.user_email.send(subject = 'Job %s accepted with %d hours of wall time' % (self.job_name,self.wall_time))
            
            self.log_status('Begin at %s' % time.strftime("%A, %b/%d/%Y, %H:%M:%S", time.localtime()))
            
            self.log_status('Pipeline Version: %s' % self.version)
            
            self.log_status('Initialize environments')
            
            self.log_status('Download job data')
            
            files_transfered = self.util.sync_get_job(self.path_input,self.path_local)
            
            if len(commands)>1:
                self.log_status('Following %d commands will be executed:' % len(commands))
            else:
                self.log_status('Following %d command will be executed:' % len(commands))
                
            self.log_status('-----------------------------------------------')
            for i in range(0,len(total_commands)):
                self.log_status('['+str(i+1)+']'+'\t'+total_commands[i].replace(self.util.CLUSTER_APP_DIR,''))
                
            self.log_status('-----------------------------------------------')

            self.pbs.PBS(self.job_name,self.path_local,commands,self.tm_start,self.log_status,self.wall_time).run()
            
            #clean before uploading
            self.util.shell_exec('rm -f %s' % os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,'tmp.fail'))
            self.util.shell_exec('rm -f %s' % os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,'node.txt'))
            self.util.shell_exec('rm -f %s' % os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,'current.txt'))
            self.util.shell_exec('rm -f %s' % os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,'tmp.pbs'))
            
            if os.path.exists(os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,'tmp.done')):
                successful = True
                self.util.shell_exec('rm -f %s' % os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,'tmp.done'))
            else:
                successful = False
            
            #delete all input files
            for f in files_transfered:
                if not f=='log.txt': 
                    self.util.shell_exec('rm -f %s' % os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,f))

            #bug: log.txt was truncated from above?
            self.log_status('Upload result files')
            
            self.util.sync_put_result(os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local)+os.sep,self.path_output)
            
            #self.log_status('Result files can be found at: %s' % self.path_output.split(':')[-1])
            
            #self.log_status('Send notification email')
            
            self.log_status('End at %s' % time.strftime("%a, %b/%d/%Y, %H:%M:%S", time.localtime()))
            
            if self.user_email:
                if successful:
                    analysis_url = None
                    add_subject = None
                    if 'a' in self.how_to_send_email or 'c' in self.how_to_send_email: #All or Complete
                        self.user_email.send(subject = 'Job %s completed' % self.job_name, body_msg=analysis_url, body_file = os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,self.util.FILE_LOG),additional_subject=add_subject)
                else: 
                    #All or Failed
                    if 'a' in self.how_to_send_email or 'f' in self.how_to_send_email:
                        self.user_email.send(subject = 'Job %s failed' % self.job_name, body_file = os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local,self.util.FILE_STDERR))
                        
            #delete the RUNNING tag
            #self.util.shell_exec_remote('rm -f %s' % self.util.FILE_TAG_RUNNING,self.path_input)
        except Exception,e:
            #delete the RUNNING tag
            #self.util.shell_exec_remote('rm -f %s' % self.util.FILE_TAG_RUNNING,self.path_input)
            if str(e)=='record': #job directory is not writable, but it has been processed before.
                pass
            else:
                import traceback
                traceback.print_tb(sys.exc_info()[2])
                print 'Error:',str(e)
                self.log_status(str(e))
                if self.user_email:
                    if 'a' in self.how_to_send_email or 'e' in self.how_to_send_email: #All or Exception
                        self.user_email.send(subject = 'Job %s failed' % self.job_name, body_msg=str(e))
                    
        finally:
            #delete job folder on cluster
            self.util.shell_exec_remote('rm -f %s' % self.util.FILE_TAG_RUNNING,self.path_input)
            self.clean_job_dict()
            self.util.shell_exec('rm -fr %s' % os.path.join(self.util.CLUSTER_JOB_DIR,self.path_local))
            
            #special handling for user "autoaligner" - run a script 'postcmd.sh' after this job is complete.
            
            if self.job_owner == 'autoaligner':
                self.util.shell_exec_remote('.\postcmd.sh',self.path_input)
                
                
            print '[%s] >>>>> %s (%s)' % (self.util.now(),self.path_input,str(self.job_dict[self.job_owner]))
            
            #backup&log this job's cmd.txt,stdout.txt and log.txt
            #tm = time.strftime("%Y_%b_%d_%H_%M_%S", time.gmtime())
            
            #copy the log files?
            #shell_exec('cp -f %s %s' % (os.path.join(CLUSTER_JOB_DIR,self.path_local,FILE_ERROR),os.path.join(CLUSTER_LOG_DIR,self.path_local+'_'+tm+'_'+FILE_ERROR)))
            #shell_exec('cp -f %s %s' % (os.path.join(CLUSTER_JOB_DIR,self.path_local,FILE_STDOUT),os.path.join(CLUSTER_LOG_DIR,self.path_local+'_'+tm+'_'+FILE_STDOUT)))
            #shell_exec('cp -f %s %s' % (os.path.join(CLUSTER_JOB_DIR,self.path_local,FILE_LOG),os.path.join(CLUSTER_LOG_DIR,self.path_local+'_'+tm+'_'+FILE_LOG)))
            #shell_exec('cp -f %s %s' % (os.path.join(CLUSTER_JOB_DIR,self.path_local,FILE_CMD),os.path.join(CLUSTER_LOG_DIR,self.path_local+'_'+tm+'_'+FILE_CMD)))
        
    def abort(self):
        """User initiated job termination
        """
        #delete the "a" file on the source folder
        #print '[%s] ----- %s (%s)' % (self.util.now(),self.job_name,str(self.job_dict[self.job_owner]))
        print '[%s] ----- %s ' % (self.util.now(),self.job_name)
        self.clean_job_dict()
        self.util.shell_exec_remote('rm -f %s %s' % (self.util.FILE_TAG_ABORT,self.util.FILE_TAG_RUNNING),self.path_input)
        #cmd  = "qstat -au %s | awk '{print $1,$4}' | grep ' %s' | awk -F '.' '{print $1}'" % (CLUSTER_USER,self.job_name)
        cmd  = 'showq -w user=%s -n -v |grep %s |cut -d "/" -f1' % (self.util.CLUSTER_USER,self.job_name) 
        #print cmd
        pbs_id = self.util.shell_exec(cmd)
        if pbs_id:
            try:
                self.util.shell_exec('qdel %s' % pbs_id[0])
                time.sleep(10) #leave some time to delete the job in PBS
                
            except Exception,e:
                print 'Failed to qdel %s' % pbs_id[0]
            
            finally:
                self.util.shell_exec('rm -fr %s' % os.path.join(self.util.CLUSTER_JOB_DIR,str(self.job_index)))

    def log_status(self,message):
        """put message into FILE_LOG and upload this file to LOCAL_SERVER.
        User will watch this file to know up-to-date progress & status of his job. 
        """
        source_log = os.path.join(self.util.CLUSTER_JOB_DIR,
                                  self.path_local,
                                  self.util.FILE_LOG)
        dest_path  = self.path_output 
        #dest_path  = os.path.join(LOCAL_SERVER_JOB_PATH,self.path_output)+os.sep
        elapsed = time.time()-self.tm_start
        hours = int(elapsed/3600)
        minutes = int((elapsed-hours*3600)/60)  
        seconds = int(elapsed-hours*3600-minutes*60)
        tm_elapsed = '%03dh:%02dm:%02ds' % (hours, minutes, seconds)
        msg = '%s\t\t%s' % (tm_elapsed,message)
        file(source_log,'a+').write(msg+'\n')
        self.util.rsync(self.util.RYSNC_UPDATE_PARAMS,source_log,dest_path)
    
    """
    def log_stdout(self,message):
        #put standard output and error into FILE_STDOUT and upload this file to LOCAL_SERVER
        file_stdout_cluster = os.path.join(CLUSTER_JOB_DIR,self.path_local,FILE_STDOUT)
        file_stdout_local = os.path.join(LOCAL_SERVER_JOB_PATH,self.path_input)
        #shell_exec('echo "%s">>%s' % (message,file_stdout_cluster))
        file(file_stdout_cluster,'a+').write(message+'\n')
    
        shell_exec('rsync %s %s %s' % (RYSNC_UPDATE_PARAMS,file_stdout_cluster,file_stdout_local))

    def log_stderr(self,message):
        #put standard error and error into FILE_STDOUT and upload this file to LOCAL_SERVER.
        file_stderr_cluster = os.path.join(CLUSTER_JOB_DIR,self.path_local,FILE_STDERR)
        file_stderr_local = os.path.join(LOCAL_SERVER_JOB_PATH,self.path_input)
        #shell_exec('echo "%s">>%s' % (message,file_stdout_cluster))
        file(file_stderr_cluster,'a+').write(message+'\n')
        shell_exec('rsync %s %s %s' % (RYSNC_UPDATE_PARAMS,file_stderr_cluster,file_stderr_local))
    """
