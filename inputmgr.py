"""Manage (discover, pair) input files in a pipeline.

When user specified some input files in a pipeline command file, this module will 
try to find if specified input files do exist.
"""
import os,sys,re,imp

class InputManager(object):
    
    def __init__(self, path, inputs, case_control=None):
        """
        @param inputs: [string]. Examaple, ['A,B,C,D']

        """
        self.util = imp.load_source('util','util.py')
        self.job_path = path
        self.inputs = inputs
        self.paired = []
        if case_control:
            #print case_control_group
            self._pair(case_control) 
        self.reserved_files = ['cmd.txt','stdout.txt','stderr.txt','log.txt']

    #def _get_path(self,file_name):
    #    #cmd = 
    #    ms = shell_exec('file -b /home/ying/test/BIN')[0]
    #    print ms.strip('symbolic link to').strip().strip('`').strip("'")
   
    def _find(self,file_name):
        """list files on job folder that match pattern.
        pattern can be any regular expression like "*.txt.gz"
        """
        msg = self.util.shell_exec_remote('ls -lpmR1 %s' % file_name, self.job_path) #FOLDER,'*.txt'
        ret = []
        path = None
        for m in msg:
            if m:
                if m.startswith('ls: '): #this file/pattern does not exist
                    raise Exception('Cannot find input "%s"' % file_name)
                
                elif m.endswith(':'):
                    path = m.strip(':')
                    continue
                
                elif m.endswith('/'): #directory
                    continue
                
                else:
                    if not m in self.reserved_files:
                        if path:
                            ret.append(os.path.join(path,m))
                        else:
                            ret.append(m)
        
        #check out if all input files have same suffix                    
        if len({}.fromkeys([k.split('.')[-1] for k in ret]))>1:
            #raise Exception('All input files that match pattern "%s" must have same suffix. Result [%s]' % (pattern,' '.join(v)))
            raise Exception('All input files that match pattern "%s" must have same file extensions. Result [%s]' % (file_name,' '.join(ret)))
        if ret:
            return ret
        else:
            raise Exception('No input files match pattern %s' % file_name)

    
    def _pair(self, case_control):
        """
        A_1.txt.gz    A_2.txt.gz,A_3.txt.gz
        A_4.txt.gz,A_5.txt.gz A_1.txt.gz
        """
        for line in case_control:
            line = line.strip()
            #print 'LINE:',line
            if line and not line.startswith('#'):
                k = []
                j = [t.strip() for t in line.split('\t')]
                k.append(j[0])
                if not len(j)==2:
                    k.append('')
                else:
                    k.append(j[1])
                #    raise Exception('Two columns( control case) data format is required for paired analysis')
                self.paired.append((k[0].split(','),k[1].split(',')))
        #print self.paired
    
    def iter_inputs(self):
        #@return [X]
        #if len(x) >1: #case-control
        #if len(x) ==1: #single
        
        ret = []
        inputs = self.inputs.split(',')
        
        if not inputs:
            inputs = ['.'] #all files under current folder
                
        #path_of_inputs = []
        total = []
        for i in inputs:
            total.extend(self._find(i))
        
        if not self.paired:
            yield [total]#[[total]]
            
        else:      
            for i in self.paired:
                #print i
                case = []
                control = []
                control_samples = i[0]
                case_samples = i[1]
                for j in total:
                    for k1 in control_samples:
                        if os.path.basename(j) == k1:
                            control.append(j)
                    for k2 in case_samples:
                        if os.path.basename(j) == k2:
                            case.append(j)
                #if not control:
                #    raise Exception('No control samples? Check out your grouping file specified by "-s"')
                #if not case:
                #    raise Exception('No cases samples? Check out your grouping file specified by "-s"')
                #print control,case
                yield [control, case] #[[control],[case]]
