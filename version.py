"""version controlled pipelines
"""
import os,sys,imp,re
from sets import Set

class VersionConf:
    """parse the version configuration file.
     
    """
    def __init__(self, version=None):
        self.util = imp.load_source('util', 'util.py')
        self.conf = {} 
        self.path_app = []
        self.PATH_PIPELINE = self.util.CLUSTER_CONF_DIR
        self.version = version
        if not version:
            self.version = self.get_current_version()
        #self.fn = '%s/%s.txt' % (self.PATH_PIPELINE,self.version)
        self.fn = os.path.join(self.PATH_PIPELINE,'%s.txt' % self.version)
        self._parse()
        
    def _parse(self):
        m = {}
        p = {}
        for line in file(self.fn):
            line = line.strip()
            if line:
                if not line.startswith('#'):
                    
                    if line.startswith('@'):
                        pipelineName,pipelineSteps = line.split(':')
                        p[pipelineName] = [i.strip() for i in pipelineSteps.strip().split(',')] 
                        
                    else:
                        vs = line.split(':')
                        Index = vs[0].strip() 
                        AppName = vs[1].strip() 
                        AppPath = vs[2].strip()
                        self.path_app.append(AppPath)
                        Genome = vs[3].strip()
                        AppDesc = vs[4].strip()
                        AppParameters = ':'.join(vs[5:]).strip()
                        
                        m.setdefault(Index,[]).append((Index,AppName,AppPath,Genome,AppDesc,AppParameters))
                            
        for pipelineName,commands in p.items():
            for i in commands:
                j = m.get(i,None)
                if j:
                    for k in j:
                        #print k
                        #self.conf.setdefault(pipelineName,[]).append((k[0],k[1],k[2],k[5]))
                        
                        self.conf.setdefault(pipelineName,[]).append((k[0],k[1],k[2],k[3],k[5]))
                        
    def get_current_version(self):
        """get the current version of pipelines
        """
        myversion = None
        if self.version:
            return self.version
        
        else: 
            v = []
            px = re.compile('^(\d+)\.txt')
            py = re.compile('^v(\d+)')
            for root, dirs, files in os.walk(self.PATH_PIPELINE):
                for f in files:
                    m = px.findall(f)
                    n = py.findall(f)
                    if m:   
                        v.append(int(m[0]))
                    if n:   
                        myversion = int(n[0])
                        
            if v:
                if myversion in v: #find specified version
                    return str(myversion)
                else: #use the most recent version
                    return str(sorted(v, reverse=True)[0])
            else:
                raise Exception('Can not find version configuration file.')
    
    def get_app_path(self):
        """get the path to APP
        """
        ret = []
        for i in self.path_app:
            j = os.path.join(self.util.LOCAL_SERVER_DIR, self.util.PATH_APP, i)
            ret.append(j)
        return list(Set(ret))
    
    def get_data_path(self):
        """get the path to DATA
        """
        return os.path.join(self.util.CLUSTER_DATA_DIR, self.get_current_version())

    
    def get_pipeline(self, pipelineName=None):
        #print self.conf.items()
        return self.conf.get(pipelineName, None)
