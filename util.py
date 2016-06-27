import subprocess
import os
import re
import time
import math
import operator
from datetime import datetime
import ConfigParser

TEST=False
#TEST = True

config = ConfigParser.RawConfigParser(allow_no_value=True)
config.readfp(file('basic.cnf'))

CLUSTER_USER = config.get("cluster", "CLUSTER_USER")
CLUSTER_NAME = config.get("cluster", "CLUSTER_NAME")
CLUSTER_MASTER_DIR  = config.get("cluster", "CLUSTER_MASTER_DIR")
CLUSTER_MASTER_TMP_DIR = config.get("cluster", "CLUSTER_MASTER_TMP_DIR")
CLUSTER_NODE_DIR = config.get("cluster", "CLUSTER_NODE_DIR")
CLUSTER_NODE_TMP_DIR = config.get("cluster", "CLUSTER_NODE_TMP_DIR")

#can use either the tmp dir on master(larger disk capacity) or node (better performance)
if CLUSTER_NODE_TMP_DIR:
    CLUSTER_TMP_DIR = CLUSTER_NODE_TMP_DIR
elif CLUSTER_MASTER_TMP_DIR:
    CLUSTER_TMP_DIR = CLUSTER_MASTER_TMP_DIR
else:
    CLUSTER_TMP_DIR = '.'

LOCAL_SERVER_USER = config.get("localserver", "LOCAL_SERVER_USER")
LOCAL_SERVER_NAME = config.get("localserver", "LOCAL_SERVER_NAME")
LOCAL_SERVER_DIR = config.get("localserver", "LOCAL_SERVER_DIR")

PATH_APP = config.get("path", "PATH_APP")#'app'
PATH_DATA = config.get("path", "PATH_DATA")#'data'
PATH_LOG = config.get("path", "PATH_LOG")#'log'
PATH_JOB = config.get("path", "PATH_JOB")#'job2'
PATH_PIPELINE = config.get("path", "PATH_PIPELINE")#'conf'

if TEST:
    PATH_JOB = 'job2'
     
#rsync parameters
RSYNC_QUERY_PARAMS      = '-Lre ssh'
RYSNC_DOWNLOAD_PARAMS   = '-Lriue ssh' #L-Softlink, a-Archive, i-  
RYSNC_UPLOAD_PARAMS     = '-Lrue ssh'
RYSNC_UPDATE_PARAMS     = '--append -re ssh'#update-only, transfer-over-ssh

#RSYNC_QUERY_PARAMS      = '-Lae ssh'
#RYSNC_DOWNLOAD_PARAMS   = '-Laiue ssh' #L-Softlink, a-Archive, i-  
#RYSNC_UPLOAD_PARAMS     = '-Laue ssh'
#RYSNC_UPDATE_PARAMS     = '--append -re ssh'#update-only, transfer-over-ssh

RSYNC_QUERY_INTERVAL    = int(config.get("rsync", "RSYNC_QUERY_INTERVAL_SECONDS"))#10 #seconds, query remote data server for new data

FILE_CMD           = config.get("file", "FILE_CMD")#'cmd.txt' #commands of each job
FILE_STDOUT        = config.get("file", "FILE_STDOUT")#'stdout.txt' #standard output and error of each job
FILE_STDERR        = config.get("file", "FILE_STDERR")#'stderr.txt' #standard output and error of each job
FILE_LOG           = config.get("file", "FILE_LOG")#'log.txt' #log information of each job
FILE_TAG_ABORT     = config.get("file", "FILE_TAG_ABORT")#'a' #empty file as "abort a job" place holder
FILE_TAG_BEGIN     = config.get("file", "FILE_TAG_BEGIN")#'b' #empty file as "begin a job" place holder
FILE_TAG_RUNNING   = config.get("file", "FILE_TAG_RUNNING")#'r' #empty file as "running a job" place holder
FILE_TAG_DRYRUN   = config.get("file", "FILE_TAG_DRYRUN")#'d' #empty file as "dry-run a job" place holder

PBS_QUERY_INTERVAL = int(config.get("pbs", "PBS_QUERY_INTERVAL_SECONDS"))#10 #seconds. Query if a job was completed in this interval.
MAX_WALLTIME_HOURS = int(config.get("pbs", "MAX_WALLTIME_HOURS"))#240 #hours
CLUSTER_DOWNTIME= config.get("pbs", "CLUSTER_DOWNTIME")#None

if CLUSTER_DOWNTIME:
    try:
        x = int((datetime(*[int(i) for i in CLUSTER_DOWNTIME.split('/')]) - datetime.now()).total_seconds()/3600)
        if x < MAX_WALLTIME_HOURS:
            MAX_WALLTIME_HOURS = x
    except:
        pass

#PRIORITY_USER = {}
#for u in config.get("priority", "PRIORITY_USER").split(','):
#    PRIORITY_USER[u] = 20

#PRIORITY_USER_WALLTIME=int(config.get("priority", "PRIORITY_USER_WALLTIME"))
#NUM_JOBS_MAX = int(config.get("priority", "NUM_JOBS_MAX"))
#NUM_CLUSTER_NODES_AVAILABLE = int(config.get("priority", "NUM_CLUSTER_NODES_AVAILABLE"))#14
#NUM_JOBS_UNLIMITED_WALLTIME = int(config.get("priority", "NUM_JOBS_UNLIMITED_WALLTIME"))#3
#HOURS_LIMITED_WALLTIME = int(config.get("priority", "HOURS_LIMITED_WALLTIME"))#72

#PRIORITY_USER_WALLTIME=240
#NUM_JOBS_MAX = 6
NUM_CLUSTER_NODES_AVAILABLE = 14
NUM_JOBS_UNLIMITED_WALLTIME = 3
HOURS_LIMITED_WALLTIME = 72

SMTP_SERVER= config.get("misc", "SMTP_SERVER")
RE_PAIRED_END_READ_FILE_1 = config.get("misc", "RE_PAIRED_END_READ_FILE_1")#'(.*)_1\.(.+)'
RE_PAIRED_END_READ_FILE_2 = config.get("misc", "RE_PAIRED_END_READ_FILE_2")#'(.*)_2\.(.+)'

##################################################
CLUSTER_APP_DIR = os.path.join(CLUSTER_MASTER_DIR,PATH_APP)
CLUSTER_JOB_DIR = os.path.join(CLUSTER_MASTER_DIR,PATH_JOB)
CLUSTER_LOG_DIR = os.path.join(CLUSTER_MASTER_DIR,PATH_LOG)
CLUSTER_DATA_DIR = os.path.join(CLUSTER_MASTER_DIR,PATH_DATA)
CLUSTER_PIPELINE_DIR = os.path.join(CLUSTER_MASTER_DIR, PATH_PIPELINE)

##################################################
LOCAL_SERVER_PATH = LOCAL_SERVER_USER+'@'+LOCAL_SERVER_NAME+':'

LOCAL_SERVER_APP_PATH  = LOCAL_SERVER_PATH+os.path.join(LOCAL_SERVER_DIR,PATH_APP)
LOCAL_SERVER_DATA_PATH = LOCAL_SERVER_PATH+os.path.join(LOCAL_SERVER_DIR,PATH_DATA)
LOCAL_SERVER_JOB_PATH  = LOCAL_SERVER_PATH+os.path.join(LOCAL_SERVER_DIR,PATH_JOB)
LOCAL_SERVER_CONF_PATH  = LOCAL_SERVER_PATH+os.path.join(LOCAL_SERVER_DIR,PATH_CONF)



LOCAL_SERVER_APP_DIR  = os.path.join(LOCAL_SERVER_DIR,PATH_APP)
LOCAL_SERVER_DATA_DIR = os.path.join(LOCAL_SERVER_DIR,PATH_DATA)
LOCAL_SERVER_JOB_DIR  = os.path.join(LOCAL_SERVER_DIR,PATH_JOB)
LOCAL_SERVER_CONF_DIR  = os.path.join(LOCAL_SERVER_DIR,PATH_CONF)

##################################################
#for priority, lower is higher
UserPriority = {'bluegill':(0,1024,240), #Admin. (Priority, Max_Jobs, Max_Wallhours)
                 'user1':(1,1024,240), 
                 'user2':(2,1024,240), 
                 'user3':(10,1024,72)}

UserPriorityDefault = (9,3,240)

CLUSTER_R_PATH = '/nfs/bluegill/app/R/lib64/R/bin/R'


##################################################

def get_user_priority(user_name):
    ps = UserPriority.get(user_name,None)
    if ps:
        p,n,h = ps
    else:
        p,n,h = UserPriorityDefault
        
    if h>MAX_WALLTIME_HOURS:
        return p,n,MAX_WALLTIME_HOURS
    else:
        return p,n,h

def extract_email_from_cmd(job_name):
    c = 'cat cmd.txt | grep -E -o \"\\b[a-zA-Z0-9.-_]+@[a-zA-Z0-9.-]+\.[a-zA-Z0-9.-]+\\b\"'
    ret = shell_exec_remote(c,job_name)
    if ret:
        return ret[0]
    else:
        return None

def accept_pending_job(job_name,job_rank):
    #1. remove the b
    #2. make a new rank or change existing rank
    #R14: empty
    px = re.compile('(^R\d+):\sempty')
    c1 = 'rm -f b'
    ret = shell_exec_remote(c1,job_name)
    c2 = 'file R*'
    ret = shell_exec_remote(c2,job_name)
    if ret:
        rx = px.findall(ret[0])
        if rx:
            c3 = 'mv %s R%d' % (rx[0],job_rank)
        else:
            c3 = 'touch R%d' % (job_rank)
        ret = shell_exec_remote(c3,job_name)
    else:
        c3 = 'touch R%d' % (job_rank)
        ret = shell_exec_remote(c3,job_name)

def shell_exec(cmd,shell=True):
    """cmd is a string!"""
    #ALL outputs will be redirected to stdout
    #do not use "time command" to wrap any commands. It is going to disrupt.
    #print cmd
    p = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,shell=shell)
    std_out,std_err = p.communicate()
    if std_out:
        return std_out.strip().split('\n')
    else:
        return []

def shell_exec_remote(cmd,job_name='.',shell=True):
    """execute a command on a remote machine using ssh
    ssh ying@111.111.111.111
    cd PATH_JOB/job_name
    run cmd
    """
    command = """ssh %s@%s 'cd %s; %s'""" % (LOCAL_SERVER_USER,LOCAL_SERVER_NAME,os.path.join(LOCAL_SERVER_JOB_DIR,job_name),cmd)
    #print command
    return shell_exec(command)

def now():
    """
    @return: Wen, Sep/15/2010 09:00:12
    """
    return time.strftime("%A, %b/%d/%Y, %H:%M:%S", time.localtime())



##start: use by GNomEx
def get_job_owner(job_name):
    cmd = """ssh %s@%s 'stat -c %%U %s'""" % (LOCAL_SERVER_USER,LOCAL_SERVER_NAME,os.path.join(LOCAL_SERVER_JOB_DIR,job_name))
    msg = shell_exec(cmd)
    if msg:
        return msg[0]
    else:
        return None

def get_job_priority(job_owner):
    return PRIORITY_USER.get(job_owner,100)

def get_job_walltime(job_priority):
    #return PRIORITY_USER_WALLTIME.get(str(job_priority),24)
    return MAX_WALLTIME_HOURS

def get_associated_lab(job_owner):
    m = {'u0123456':'Bioinformatics Core Facility',}
    #return m.get(job_owner,None)
    return '"Brad Cairns"'

def get_labs():
    labs = {}
    for line in file('labs.txt'):
        line = line.strip()
        try:
            toks = line.split('\t')
            labs[toks[0]]= toks[1]
        except:
            pass
    return labs

def get_organism(genome_build):
    """
Arabidopsis
Aspergillus fumigatus
Asteromyia carbonifera
BioMicro
C. elegans
Canis familiaris
Chicken
Control
Drosophila
E. coli
Emergen
Grape
Human
Mouse
Mycobacterium abscessus
Mycobacterium tuberculosis
Newt
Other
Pig
Planaria
Rat
S. pombe
Sheep
Sodalis
Unknown
Xenopus
Yeast
Yeast Intergenic
Zebrafish
    """
    p = genome_build.lower()
    if p.startswith('hg'):
        organism ='Human'
    elif p.startswith('mm'): 
        organism ='Mouse'
    elif p.startswith('zv'): 
        organism ='Zebrafish'
    elif p.startswith('ce'): 
        organism ='C. elegans'
    elif p.startswith('dm'): 
        organism ='D. melanogaster'
    elif p.startswith('ye'): 
        organism ='Yeast'
    elif p.startswith('sp'): 
        organism ='S. pombe'
    else:
        organism ='Unknown'
    #else:
    #    raise Exception('Can not recognize the organism of genome build %s \nAt now only hgXX (Example: hg19 -> Human), mmXX (for Mouse) and zvXX( for Zebrafish) are supported' % genome_build)
    return organism
 
def get_genomebuild(genome_build):
    p = genome_build.lower()
    ret = 'hg19'
    if p.startswith('hg19'):
        ret ='hg19'
    elif p.startswith('hg18'): 
        ret ='hg18'
    elif p.startswith('mm8'): 
        ret ='mm8'
    elif p.startswith('mm9'): 
        ret ='mm9'
    elif p.startswith('mm10'): 
        ret ='mm10'
    elif p.startswith('zv8'): 
        ret ='zv8'
    elif p.startswith('zv9'): 
        ret ='zv9'
    elif p.startswith('rn4'): 
        ret ='rn4'
    elif p.startswith('ce6'): 
        ret ='ce6'
    elif p.startswith('ce10'): 
        ret ='ce10'
    elif p.startswith('dm3'): 
        ret ='dm3'
    elif p.startswith('ye3'):
        ret ='ye3'
    elif p.startswith('sp07'):
        ret ='sp07'
    else:
        ret ='unknown'
        #raise Exception('Can not recognize the genome build %s \nAt now only hg19, hg18, mm9, mm8, zv9, zv8, ce6, ce10 and ye3 are supported' % genome_build)
    
    return ret 

def get_gnome_desc(genomebuild):
    m = {'ce10':'C_elegans_Oct_2010','hg19':'H_sapiens_Feb_2009','mm9':'M_musculus_Jul_2007'}
    return m.get(genomebuild,genomebuild)
    
def make_directories():
    if not os.path.exists(CLUSTER_APP_DIR):
        shell_exec('mkdir -p %s' % CLUSTER_APP_DIR)

    if not os.path.exists(CLUSTER_JOB_DIR):
        shell_exec('mkdir -p %s' % CLUSTER_JOB_DIR)

    if not os.path.exists(CLUSTER_LOG_DIR):
        shell_exec('mkdir -p %s' % CLUSTER_LOG_DIR)

    if not os.path.exists(CLUSTER_PIPELINE_DIR):
        shell_exec('mkdir -p %s' % CLUSTER_PIPELINE_DIR)
        
    if not os.path.exists(CLUSTER_DATA_DIR):
        shell_exec('mkdir -p %s' % CLUSTER_DATA_DIR)

    if not os.path.exists(CLUSTER_TMP_DIR):
        shell_exec('mkdir -p %s' % CLUSTER_TMP_DIR)

def clean_job_directory():
    if os.path.exists(CLUSTER_JOB_DIR):
        shell_exec('rm -fr %s/*' % CLUSTER_JOB_DIR)

def make_remote_directory(path):
    """create a directory on LOCAL_SERVER"""
    cmd = "ssh %s@%s 'mkdir -p %s'" % (LOCAL_SERVER_USER, LOCAL_SERVER_NAME, path)
    ret = shell_exec(cmd)
    if ret: #error:
        raise Exception('Failed to create directory %s: %s' % (path,ret)) 

def test_writeable_directory(path):
    #The output folder on LOCAL_SERVER must be writable for user "hiseq"."""
    cmd_create_tmp_file = """if [ `touch testfile 2> /dev/null; echo "$?"` -eq 0 ]; then
        echo 0
        rm testfile
        else
        echo 1
        fi"""
    #msg = shell_exec_remote(cmd_create_tmp_file,job_name=path)
    
    msg = shell_exec_remote(cmd_create_tmp_file,job_name=path)
    if len(msg)>1:
        raise Exception("Directory %s does not exist" % path)

    failed = int(msg[0])
    if failed:
        raise Exception("Directory %s is not writable" % path)


def rsync(param,source,dest):
    return shell_exec('rsync %s %s %s' % (param,source,dest))
    
def sync_query():
    #query pending jobs
    cmd = """rsync %s %s | awk '{print $1 "?" $3"#"$4 "?" $5}'""" % (RSYNC_QUERY_PARAMS,LOCAL_SERVER_JOB_PATH+os.sep)
    
    #Link of folder may not work
    #readlink -f FILE
    
    fnames = shell_exec(cmd)
    begin_jobs = []
    abort_jobs = []
    if fnames:
        df = {}
        dt = {}
        for fname in fnames:
            try:
                [prop,ctime,fn] = [f.strip() for f in fname.split('?')]
                if not fn:
                    pass
                elif fn.startswith('.') :
                    pass
                elif fn.startswith('..') :
                    pass
                else:
                    #if prop[0] == 'd': #directory
                    #    pass
                    #if prop[0] == 'l': #link
                    #    pass
                    if prop[0] == '-': #file
                        if not os.path.basename(fn).startswith('.'): #not a hidden file
                            parent_path = os.path.dirname(fn)
                            df.setdefault(parent_path,[]).append(os.path.basename(fn))
                            #Transform time to seconds, for comparing which job should be run firstly. First In First Run.   
                            #'2010/08/20#16:27:40' -> 1282343260.0
                            dt[parent_path] = int(time.mktime(time.strptime(ctime,'%Y/%m/%d#%H:%M:%S')))
            except Exception,e:
                """
                [prop,ctime,fn] = [f.strip() for f in fname.split('@!@')]
                ValueError: need more than 1 value to unpack
                
                rsync: opendir "/mnt/..." failed: Permission denied (13)
                
                """
                pass
                    
        for k,v in df.items():
            
            if (FILE_TAG_ABORT,FILE_TAG_BEGIN) in v:
                return begin_jobs,abort_jobs,
            
            if FILE_TAG_ABORT in v: #the user wants to abort this job after submission.
                abort_jobs.append((k,dt[k]))
                    
            #if FILE_TAG_BEGIN in v or FILE_TAG_RUNNING in v: #the job is ready to run
            if FILE_TAG_BEGIN in v or FILE_TAG_RUNNING in v or FILE_TAG_DRYRUN in v: #the job is ready to run
                #and it has a "cmd.txt" file
                if FILE_CMD in df[k]:
                    begin_jobs.append((k,dt[k]))

                    
        return (sorted(begin_jobs,key=operator.itemgetter(1)),sorted(abort_jobs,key=operator.itemgetter(1)))
    
    else:
        return begin_jobs,abort_jobs


def sync_get_job(path_remote,path_local,files=[]):
    #download all input files for one job
    #@path_remote: A/B/C
    #@files: ['1.txt','2.txt']
    #@local_path: A_123345
    dest = os.path.join(CLUSTER_JOB_DIR,path_local) + os.sep
    msg = []
    if files: #download specified files only
        for f in files:
            #LOCAL_SERVER_JOB_PATH -> '/bluegill/job/'
            source = os.path.join(LOCAL_SERVER_JOB_PATH,path_remote,f)
            
            if f.startswith('/'): #full path to a file
                source = f
                
            rsync(RYSNC_DOWNLOAD_PARAMS, source, dest)
        source = os.path.join(LOCAL_SERVER_JOB_PATH,path_remote,FILE_LOG)
        msg = rsync(RYSNC_DOWNLOAD_PARAMS, source, dest)
        
    else: #if no files are specified, then download all the files under that directory.
        source = os.path.join(LOCAL_SERVER_JOB_PATH,path_remote) + os.sep
        msg = rsync(RYSNC_DOWNLOAD_PARAMS, source, dest)
    
    """
    .d..tp... ./
    >f+++++++ 7550X3_A_1_1.txt.gz
    >f+++++++ 7550X3_A_1_2.txt.gz
    cd+++++++ 1/
    >f+++++++ 1/hello.txt
    cd+++++++ 2/
    >f+++++++ 2/world.txt
    """
    files_transfered = []
    for m in msg:
        if m.startswith('>f'):
            try:
                fn = m.split()[1]
                if not fn in (FILE_LOG,FILE_STDOUT,FILE_STDERR):
                    files_transfered.append(m.split()[1])
            except:
                pass
    return files_transfered
        

def sync_get_cmd(path_remote,path_local):
    #download the "cmd.txt" file for job preprocesing.
    source = os.path.join(LOCAL_SERVER_JOB_PATH, path_remote, FILE_CMD)
    dest = os.path.join(CLUSTER_JOB_DIR,path_local)+os.sep
    rsync(RYSNC_DOWNLOAD_PARAMS, source, dest)

def sync_put_result(from_cluster, to_local):
    #upload result back to local server
    return rsync(RYSNC_UPLOAD_PARAMS, from_cluster, to_local)

def sync_get_app():
    #update APP at cluster
    source = LOCAL_SERVER_APP_PATH + os.sep
    dest = CLUSTER_APP_DIR + os.sep
    return rsync(RYSNC_UPLOAD_PARAMS, source, dest)
    
def sync_get_data():
    #update DATA at cluster
    source = LOCAL_SERVER_DATA_PATH + os.sep
    dest = CLUSTER_DATA_DIR + os.sep
    return rsync(RYSNC_UPLOAD_PARAMS, source, dest)

def sync_get_pipeline():
    #update CONF at cluster
    source = LOCAL_SERVER_CONF_PATH + os.sep
    dest = CLUSTER_PIPELINE_DIR + os.sep
    return rsync(RYSNC_UPLOAD_PARAMS, source, dest)

def update_cmd_params(cmd,reserved,new_params,first_param_index=0):
    """cmd -> ''
    reserved = []
    new_params->'-r " ALL " '
    
    cmd = 'novoalign -r All -d hg19.nix -k -f f1.txt f2.txt -t 60 90 >1.sam 2> &1'
    reserved = ['-f','-d']
    new_params = '-r " Random" -t 30 -d mm9.nix'
    """
    #check 
    #print new_params
    cmd_list = filter(None,[i for i in cmd.split()])
    current_key = None 
    current_value = [] 
    m = {}
    keeped = []
    skip_next = False
    #make dict on original command
    for i in range(first_param_index,len(cmd_list)):
        t = cmd_list[i]
        if skip_next:
            keeped.append(t)
            skip_next = False
            continue
        try:
            #print t
            if t.startswith('-') or t.startswith('--'):
                if current_key: #the previous key has no value, like a on/off param
                    #print current_key,current_value #  
                    m[current_key] = current_value #
                    #params.append((current_key,current_value))
                    current_value = []
                current_key = t
                
            elif t.startswith('|'):
                keeped.append(t)
                if len(t)==1:
                    skip_next = True
                
            elif t.startswith('>'):
                keeped.append(t)
                if len(t)==1:
                    skip_next = True
            elif t.startswith('1>'):
                keeped.append(t)
                if len(t)==2:
                    skip_next = True
            elif t.startswith('2>'):
                keeped.append(t)
                if len(t)==2:
                    skip_next = True
            else:
                current_value.append(t)
        except:
            pass
    if current_key:
        m[current_key] = current_value
        #params.append((current_key,current_value))
    #make dict on new parameters
    n = normalize_params(new_params)
    for r in reserved:
        try: 
            del n[r]
        except:
            pass
    m.update(n)
    ps = []
    for k,v in m.items():
        if v:
            ps.append(str(k)+' '+' '.join(v))
        else:
            ps.append(str(k))
     
    if first_param_index>0:
        new_cmd = ' '.join([' '.join(cmd_list[:first_param_index]),' '.join(ps),' '.join(keeped)])
    else:
        new_cmd = ' '.join([' '.join(ps),' '.join(keeped)])
    return new_cmd

def normalize_params(p,token='"'):
    m = {}
    key = None
    value = []
    quotes = p.split(token)
    if len(quotes)==1: #no token
        pass
    elif len(quotes)%2 == 0: #one token, exception
        raise Exception('Malformatted parameter: %s' % p)
    else:
        pass
    
    string_param = []
    
    p1 = filter(None,[i.strip() for i in p.split()])
    #print p1
    for i in p1:
        #print '#',i
        if i.startswith('-'):
            if not token in i:
                if key:
                    m[key] = value
                    value = []
                key = i
            else:
                raise Exception('Malformatted parameter: %s' % i)

        elif i.startswith(token): #
            if len(i)==1:#could be starting " or ending "
                if string_param: #this is an ending "
                    string_param.append(token)
                    value.append(''.join(string_param))
                    string_param = []
                else: #this is an starting "
                    string_param.append(token)
            else: #"XXXX
                #print i
                if i.endswith(token): #"XXXX"
                    string_param.append(i)
                    value.append(''.join(string_param))
                    string_param = []
                else: #just "xxxx
                    string_param.append(i)
                
        elif i.endswith(token): #xxx"
            if string_param:
                string_param.append(i)
                value.append(''.join(string_param))
                string_param = []
            else:
                raise Exception('Malformatted parameter: %s' % i)
        else: #general 
            if string_param:
                string_param.append(str(i))
            else:
                value.append(i)
    m[key] = value
    special = {}.fromkeys(['>','<','|'])
    for k,v in m.items():
        for vv in v:
            if vv.startswith(token):
                pass
            else:
                for c in vv:
                    if special.has_key(c):
                        raise Exception('Illegal parameter %s' % vv)
    return m

def get_next_jobdir(basename, basepath='.'):
    ds = []
    for root, dirs, files in os.walk(basepath):
        for d in dirs:
            try:
                #if d==basename:
                #    return d+str(1)
                #else: 
                if d.startswith(basename):
                    ds.append(int(d.lstrip(basename)))
            except:
                pass
            
    sds = sorted(ds,reverse=True)
    if sds:
        job_folder_index = sds[0]+1
    else:
        job_folder_index = 1
    return basename+str(job_folder_index)
    #make new job dir
    #myfolder = os.path.join(PATH_JOB,uname,str(job_folder_index))
    #if not os.path.exists(myfolder):
    #    shell_exec('mkdir -p %s' % myfolder)
    #else:
    #    shell_exec('rm -fr %s/*' % myfolder)

#def make_command(versionNumber,pipelineName):
#print get_version('v1','@snpindel -g hg19')
        
"""
def get_free_nodes():
    free_nodes=0
    nodes = range(241,253)
    total_nodes = len(nodes)
    our_nodes=['em'+str(i) for i in nodes]
    msg = shell_exec('qnodes')
    nextline_is_what_we_want = False
    current_node=''
    for line in msg:
        try:
            if nextline_is_what_we_want:
                status = line.strip().split('=')[1].strip()
                if status=='free':
                    free_nodes += 1 
                nextline_is_what_we_want = False
             
            for n in our_nodes:
                if line.strip().startswith(n):
                    current_node = n
                    nextline_is_what_we_want = True
        except:
            return 'Nodes status are unknown'
    #print '%d of %d nodes are free now' % (free_nodes,total_nodes)
    return '%d nodes are free now' % free_nodes

def check_security(fn_cmd):
    #command = ['ls','-l']
    #check the command safety to avoid things like "rm -fr /home/*"
    #1. I am not a root user, so blocking commands like "mkfs" are not necessary.
    #2. User can only operate her/his own folder
    #
    
    limited_path_operations = {}.fromkeys(['.','..','/','$'])
    limited_commands={}.fromkeys(['sh','bash','csh','ksh','zsh','tcsh'])
    command = []
    for line in file(fn_cmd):
        
        if not line[:-1].strip(): #blank line
            continue
        
        if not line.startswith('@') and not line.startswith('#'):
            if line.startswith('!'):
                command = filter(None, [s.strip() for s in line[1:].split(' ')])
            else:
                command = filter(None, [s.strip() for s in line.split(' ')])
                
            if limited_commands.has_key(command[0]):
                raise Exception('Illegal command: %s' % command)
            else: #command seems OK, but parameters especially path must be checked.
                for c in command:
                    if limited_path_operations.has_key(c[0]): #the first letter of any commands and parameters
                        raise Exception('Illegal command: %s' % command)


def map_job_path(job_name,realpath=True):
    if realpath:
        return '/'.join(job_name.split('_'))
    else:
        return '_'.join(job_name.split('/')) #'A/B/C' ---> 'A_B_C', 'A'--->'A'

def tarball(path,postfix='.sam',gz=True):
        #merge all files under path which has .fq postfix.
        files = []
        for root, dirs, files in os.walk('.'):
            for f in files:
                if f.endswith(postfix):
                    files.append(os.path.join(root,f))
        if gz:
            opt = '-zcf'
        else:
            opt = '-cf'
        cmd = 'tar %s %s %s' % (opt,'hello.tar.gz', '*.'+postfix) 
        shell_exec(cmd)
        return files

def query_files(path_remote,suffix='.txt.gz'):
    ret = []
    source = os.path.join(LOCAL_SERVER_JOB_DIR,path_remote) + os.sep
    cmd = "rsync %s %s | grep %s | awk '{print $5}' " % (RSYNC_QUERY_PARAMS,source,suffix)
    files = shell_exec(cmd)
    for line in files: 
        ret.append(line)
    return ret

"""


