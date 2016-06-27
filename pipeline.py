import os
import imp
from inputmgr import InputManager  
from appmgr import AppManager  

class PipelineBuilder:
    
    def __init__(self,job_source,commmand,version=None):
        self.util = imp.load_source('util','util.py')
        self.app = imp.load_source('app','app.py')
        self.env = imp.load_source('env','env.py')

        #self.metacmds = filter(None, [s.strip() for s in commmand.lstrip('@').split(' ')])
        self.metacmds = filter(None, [s.strip() for s in commmand.split(' ')])
        #print self.metacmds
        self.reserved_files = [self.util.FILE_CMD,self.util.FILE_STDOUT,self.util.FILE_STDERR,self.util.FILE_LOG]
        self.job_source=job_source
        self.version=version

    def get_commands(self):
        
        bool_reserve_intermediate_files = False 
        bool_make_bam=False
        bool_gzip_sam=False
        aligner = None
        tmp = []
        navs = []
        nested=False
        
        pipeline_name = self.metacmds[0]
        
        for i in range(1,len(self.metacmds)):
            if self.metacmds[i].startswith("["):
                if self.metacmds[i].endswith("]"): #single param like [-k]
                    navs.append('#'+self.metacmds[i][1:-1])
                else:
                    tmp.append('#'+self.metacmds[i][1:])
                    nested = True
            else:
                if nested:
                    if self.metacmds[i].endswith("]"):#multiple param like [-k -r None]
                        nested = False
                        tmp.append(self.metacmds[i][:-1])
                        navs.append(' '.join(tmp))
                        tmp = []
                    else:
                        tmp.append(self.metacmds[i])
                else:
                    navs.append(self.metacmds[i])
        
        kcmds = {}
        key = None
        values = []
        for t in navs:
            if t.startswith('-'):
                if key:
                    kcmds[key] = values 
                    values = []
                key = t
            elif t.startswith('#'):
                values.append(t[1:])
            else:
                values.append(t)
                
        kcmds[key] = values
        
        #if not kcmds.has_key('-g'):
        #    raise Exception("No genome : use -g")
        #else:
        #    genome=kcmds['-g'][0]
        
        #platform = kcmds.get('-p','illumina')
        
        #if kcmds.has_key('-r'):
        #    bool_reserve_intermediate_files = True
            
        #if kcmds.has_key('-bam'):
        #    bool_make_bam = True
            
        #input        
        inputs = []

        case_control_group = None
        grouping_file = kcmds.get('-s',None)
        if grouping_file:
            case_control_group = self.util.shell_exec_remote('cat %s' % grouping_file[0],self.job_source)
            
        bool_multiple_bams = False
        if not kcmds.has_key('-i'):
            raise Exception('No input files : use "-i INPUT"')
        else:
            im = InputManager(self.job_source,kcmds['-i'][0],case_control_group)

        total_commands = []
        for finput in im.iter_inputs():
            total_commands.extend(AppManager(self.version,pipeline_name,finput,kcmds).get_command())
        return total_commands
