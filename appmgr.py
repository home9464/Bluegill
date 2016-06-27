"""
"""
import os,sys,imp,re,inspect,time
from version import VersionConf

class AppManager(object):
    def __init__(self,version,pipeline_name,inputs,global_params={}):
        """
        ver(int): requested version
        genome(str): genome, hg19, mm9, etc.
        platform(str): illumina, solid, etc.
        links([], steps):
        
        inputs:
        1 A_1.txt.gz, A_2.txt.gz
        2 .fq.gz
        3 .bam
        4 .sam
        5 .sam.gz
        5 .vcf
        """
        self.util = imp.load_source('util','util.py')
        
        #add to commands
        self.commands = []
        
        vobj = VersionConf(version)
        
        self.apppath = vobj.get_app_path()
         
        self.version = vobj.get_current_version()
        
        self.datapath = vobj.get_data_path() 
        
        self.links = vobj.get_pipeline(pipeline_name)
        #print self.links
        if not self.links:
            raise Exception('Pipeline [%s] was not found in version %s' % (pipeline_name,version))
        
        
        self.platform = 'illumina'
        if global_params.has_key('-p'):
            self.platform = global_params['-p'][0]
        
        
        if not global_params.has_key('-g'):
            raise Exception("No genome : use -g")
        else:
            self.genome = global_params['-g'][0]
        
        #change the products in-place
        #[(group1),(group2),(group3),...]
        #print 'INPUT:',inputs
        
        self.products = inputs #[[]], it is a nested list , "list of list" 
        
        self.control_case = False
        if len(inputs)==2: #case-control
            self.control_case = True
           
        self.sys_params = global_params

        self.delete_intermediate_file = True
        if self.sys_params.has_key('-r'):
            self.delete_intermediate_file = False

        #self.gatk_call_variants_with_multiple_bams = False
        #if self.sys_params.has_key('-m'):
        #    self.gatk_call_variants_with_multiple_bams = True

        self.KNOWN_INDEL_INTERVALS  = '%s.indel.intervals' % self.genome

    def _normalize_app_name(self,appname):
        """
        GenomeAnalysisTK.jar -T UnifiedGenotyper
        bwa aln
        VarScan.jar somaticFilter
        """
        ret = [i for i in appname.split(' ') if i.strip()]
          
    def get_command(self):
        t = {}
        k = []
        for i in self.links:
            t.setdefault(i[0],[]).append(i[1:])
            if k:
                if k[-1] == i[0]:
                    continue  
            k.append(i[0])
        #print t['9']
        for j in k: #for each link
            cmds = t[j]
            appName, appPath, Genome,appParams = cmds[0]
            app = appName.split(' ')[0].replace('.','_') #app.jar -> app_jar
            getattr(AppManager,app)(self,cmds)
            
        return self.commands
    
    def get_product(self):
        return self.products
    
    def pair_fastq(self,files):
        """@return: [[A_1,A_2],[B_1,B_2],...]"""
        k1 = {}
        fs = sorted(files)
        if len(fs)==1:
            return [fs]
        for f in fs:
            #m1 =  re.match('(.*)_1\.(.+)',f)
            #m2 =  re.match('(.*)_2\.(.+)',f)
            m1 =  re.match(self.util.RE_PAIRED_END_READ_FILE_1,f)
            m2 =  re.match(self.util.RE_PAIRED_END_READ_FILE_2,f)
            if m1 or m2:
                if m1:
                    k1.setdefault(''.join(m1.groups()),[]).append(f) 
                if m2:
                    k1.setdefault(''.join(m2.groups()),[]).append(f)
            else: #any other files that do not end with '1.x' or '2.x'.
                fn = '.'.join(f.split('.')[:-1])
                k1.setdefault(fn,[]).append(f)
        return k1.values()

    
    def get_file_type(self):
        for i in self.products:#i is []
            for k in i:
                j = k.lower()
                if j.endswith('.txt'):
                    return 'FASTQ'
                elif j.endswith('.txt.gz'):
                    return 'FASTQ'
                elif j.endswith('.fq'):
                    return 'FASTQ'
                elif j.endswith('.fq.gz'):
                    return 'FASTQ'
                elif j.endswith('.fastq'):
                    return 'FASTQ'
                elif j.endswith('.fastq.gz'):
                    return 'FASTQ'
                elif j.endswith('.sam'):
                    return 'SAM'
                elif j.endswith('.sam.gz'):
                    return 'SAM'
                elif j.endswith('.bam'):
                    return 'BAM'
                elif j.endswith('.vcf'):
                    return 'VCF'
                elif j.endswith('.annovar'):
                    return 'ANNOVAR'
                elif j.endswith('.pileup'):
                    return 'PILEUP'
                elif j.endswith('.mpileup'):
                    return 'MPILEUP'
                elif j.endswith('.gvf'):
                    return 'GVF'
                elif j.endswith('.cdr'):
                    return 'CDR'
                else:
                    raise Exception("File format was not recognized")
            
    ########################
    #Novoalign
    ########################
    def _get_app(self,conf,required_file_type):
        """
        myname = inspect.stack()[1][3] # get parent's name -> who is calling me?'
        """
        if not self.products:
            raise Exception("Can not find any input files")
        
        t = self.get_file_type()
        
        if required_file_type:
            if not t in required_file_type.split(','):
                raise Exception('Application [%s] Only accept file format [%s]' % (conf[0][0],required_file_type)) 
         
        m = []
        for c in conf: #if there is grouped apps, the len(conf) > 1, otherwise len(conf)==1 
            appname,apppath,genome,appparams = c
            if genome:
                if not self.genome in genome.split(','):
                    raise Exception('[%s] can not be used for genome [%s] in version [%s], available genomes are %s' % (appname,self.genome,self.version,genome))
            x = '%s %s' % (os.path.join(apppath.replace(self.util.LOCAL_SERVER_APP_DIR,self.util.CLUSTER_APP_DIR),appname),appparams)
            m.append(x)
        if not m:
            raise Exception('Failed to get necessary application ')
        return m

    def novoalign(self,conf):
        """predefined_parameters -> novoalign,/PATH/TO/novoalign,parameters
        """
        m = self._get_app(conf,'FASTQ')
        app = m[0]
        
        sam_header= "-o SAM $'@RG\\tID:bluegill\\tPL:illumina\\tLB:libtmp\\tSM:%s\\tCN:INSTITUTION'"
        
        current = []
        for p in self.products:
            pair = self.pair_fastq(p)#i is []
            tmp = []
            genome_index_file = '%s.nov.%s.nix' % (self.genome,self.platform)
            for i in pair:#i is []
                cc = None
                if len(i)==1: #single-end sequence read file
                    f = i[0]
                    output = f.rstrip('.gz')
                    output = output.rstrip('.gzip')
                    output = output.rstrip('.txt')
                    output = output.rstrip('.fastq')
                    output = output.rstrip('.fq')
                    output = output.rstrip('.fasta')
                    output = output.rstrip('.fa')
                    cmd_string = "%s %s -d %s -f %s | gzip >%s.sam.gz"
                    cc = cmd_string %(app,sam_header%output,genome_index_file,f,output)
                
                elif len(i)==2: #pair-end sequence read file
                    f1,f2 = i
                    output = os.path.commonprefix((f1,f2)).strip('_')
                    if not output:
                        output = 'tmp'
                    cmd_string = "%s %s -d %s -f %s %s| gzip >%s.sam.gz"
                    cc = cmd_string %(app,sam_header%output,genome_index_file,f1,f2,output)

                    
                else:
                    raise Exception('Two many input files: %d' % len(self.products))
            
                if cc:
                    params = self.sys_params.get('-novoalign',None)
                    if params:
                        reserved_params = ['-d','-f','-o']
                        cc = self.util.update_cmd_params(cc,reserved_params,params[0],first_param_index=1)
                        self.commands.append('%s' % cc)
                    else:
                        self.commands.append(cc)
                    tmp.append(output+'.sam.gz')
                else:
                    pass
            if tmp:
                current.append(tmp)
        self.products = current 

    def bwa(self,conf):
        """
        3:bwa aln:/bluegill/app/bwa/1.0/::COMMENT:-t 24 
        3:bwa sampe:/bluegill/app/bwa/1.0/::COMMENT:-P

        3:bwa aln:/bluegill/app/bwa/1.0/::COMMENT:-t 24 
        3:bwa samse:/bluegill/app/bwa/1.0/::COMMENT:-P
        """
        m = self._get_app(conf,'FASTQ')
        bwaaln = m[0]
        bwase = m[1]
        bwape = m[2]
        #print m
        current = []

        sam_header = "-r '@RG\\tID:bluegill\\tPL:illumina\\tLB:libtmp\\tSM:%s\\tCN:INSTITUTION'"
        for p in self.products:
            pair = self.pair_fastq(p)#i is []
            tmp = []
            genome_index_file = '%s.bwa.%s' % (self.genome,self.platform)
            for i in pair:#i is []
                cmd_aln = ['','']
                cmd_sam = None
                if len(i)==1: #single-end sequence read file
                    f = i[0]
                    output = f.split('.')[0]
                    cmd_aln[0] = "%s %s %s > %s.sai" % (bwaaln,self.genome,f,output)
                    cmd = "%s %s %s %s.sai %s | gzip > %s.sam.gz" % (bwase,sam_header%output,genome_index_file,output,f,output)
                
                elif len(i)==2: #pair-end sequence read file
                    f1,f2 = i
                    output = os.path.commonprefix((f1,f2)).strip('_')
                    if not output:
                        output = 'tmp'
                    cmd_aln[0] = "%s %s %s > %s.1.sai" % (bwaaln,genome_index_file,f1,output)
                    cmd_aln[1] = "%s %s %s > %s.2.sai" % (bwaaln,genome_index_file,f2,output)
                    cmd = "%s %s %s %s.1.sai %s.2.sai %s %s | gzip > %s.sam.gz" % (bwape,sam_header%output,genome_index_file,output,output,f1,f2,output)
                
                else:
                    raise Exception('Two many input files: %d' % len(self.inputs))
            
                #if self.user_defined_params:
                user_params_aln = self.sys_params.get('-bwaaln',None)
                if user_params_aln:
                    reserved_params_aln = []
                    if cmd_aln[0]:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_aln[0],reserved_params_aln,user_params_aln[0],first_param_index=2))
                    if cmd_aln[1]:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_aln[1],reserved_params_aln,user_params_aln[0],first_param_index=2))
                        
                else: #no user-defined parameters
                    if cmd_aln[0]:
                        self.commands.append(cmd_aln[0])
                    if cmd_aln[1]:
                        self.commands.append(cmd_aln[1])
                
                if len(cmd_aln)==1: #bwa samse
                    user_params_samse = self.sys_params.get('-bwasamse',None)
                    reserved_params_samse = ['-r']
                    if user_params_samse:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_sam,reserved_params_samse,user_params_samse[0],first_param_index=2))
                    else:
                        self.commands.append(cmd)
                    self.commands.append('rm -fr *.sai')
                    
                else:
                    user_params_sampe = self.sys_params.get('-bwasampe',None)
                    reserved_params_sampe = ['-r']
                    if user_params_sampe:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_sam,reserved_params_sampe,user_params_sampe[0],first_param_index=2))
                        
                    else:
                        self.commands.append(cmd)
                    self.commands.append('rm -fr *.sai')
                tmp.append('%s.sam.gz' % output)
            if tmp:
                current.append(tmp)
        self.products = current 
    
    ########################
    #Picard
    ########################
    def picard(self,conf,appk,required_inputfiletype,outputfilesuffix='.bam',b_change_product=True):
        """
        """
        m = self._get_app(conf,required_inputfiletype)
        app = m[0]
        appn = os.path.basename(app.split(' ')[0])

        suffix = None
        if appk:
            if outputfilesuffix:
                suffix = '.%s.%s' % (appk,outputfilesuffix)
        else:
            if outputfilesuffix:
                suffix = '.%s' % outputfilesuffix
        current = []
        for p in self.products:
            tmp = []
            for i in p:
                px = re.compile(suffix)
                if not px.findall(i): #no need to sort again
                    fout = None
                    if i.endswith('.sam'):
                        fbase = i.rstrip('.sam')
                    elif i.endswith('.sam.gz'):
                        fbase = i.rstrip('.sam.gz')
                    elif i.endswith('.bam'):
                        fbase = i.rstrip('.bam')
                    else:
                        fbase = 'tmp'
                    if suffix:
                        fout = fbase+suffix
                    else:
                        fout = fbase
                    
                    self.commands.append('%s INPUT=%s OUTPUT=%s' % (app,i,fout))
                    if self.delete_intermediate_file:
                        if appk: #BuildBamIndex
                            self.commands.append('rm -f %s' % fout.rstrip(suffix)+'.ba*')
                        
                    if b_change_product:
                        tmp.append(fout)
                    else:
                        tmp.append(i)
                
                else: #no need to do again. SortBam will not work on "sort.bam" 
                    tmp.append(i)
            current.append(tmp)
        self.products=current  

    ########################
    #Picard
    ########################
    def SortSam_jar(self,conf):
        appk = 'sort'
        filetype='BAM,SAM'
        outsuffix = 'bam'
        self.picard(conf,appk,filetype,outsuffix)
        
    ########################
    #Picard
    ########################
    def FixMateInformation_jar(self,conf):
        appk = 'mate'
        filetype='BAM,SAM'
        outsuffix = 'bam'
        self.picard(conf,appk,filetype,outsuffix)

    ########################
    #Picard
    ########################
    def MarkDuplicates_jar(self,conf):
        appk = 'dup'
        filetype='BAM'
        outsuffix = 'bam'
        m = self._get_app(conf,filetype)
        app = m[0]
        appn = os.path.basename(app.split(' ')[0])

        suffix = '.%s.%s' % (appk,outsuffix)
                        
        current = []
        for p in self.products:
            tmp = []
            for i in p:
                px = re.compile(suffix)
                if not px.findall(i): #no need to sort again
                    fout = None
                    if i.endswith('.sam'):
                        fbase = i.rstrip('.sam')
                    elif i.endswith('.sam.gz'):
                        fbase = i.rstrip('.sam.gz')
                    elif i.endswith('.bam'):
                        fbase = i.rstrip('.bam')
                    else:
                        fbase = 'tmp'
                    
                    fout = fbase+suffix
                    
                    self.commands.append('%s INPUT=%s OUTPUT=%s M=%s.dupmetrics' % (app,i,fout,fbase))
                    if self.delete_intermediate_file:
                        if appk: #BuildBamIndex
                            self.commands.append('rm -f %s' % fout.rstrip(suffix)+'.ba*')
                    tmp.append(fout)
                else: #no need to do again. SortBam will not work on "sort.bam" 
                    tmp.append(i)
            current.append(tmp)
        self.products = current
        
    
    ########################
    #Picard
    ########################
    def BuildBamIndex_jar(self,conf):
        appk = None
        filetype='BAM'
        outsuffix = 'bai'
        #no changes on existing products
        self.picard(conf,appk,filetype,outsuffix,False)
        
    ########################
    #Picard
    ########################
    def CalculateHsMetrics_jar(self,conf):
        appk = None
        filetype='BAM'
        outsuffix = 'matrics'
        #no changes on existing products
        self.picard(conf,appk,filetype,outsuffix,False)

    ########################
    #Picard
    ########################
    def CollectInsertSizeMetrics_jar(self,conf):
        appk = None
        filetype='BAM'
        outsuffix = 'InsertSizeMetrics'
        #no changes on existing products
        self.picard(conf,appk,filetype,outsuffix,False)

    ########################
    #Picard
    ########################
    def CollectMultipleMetrics_jar(self,conf):
        appk = None
        filetype='BAM'
        outsuffix = None
        #no changes on existing products
        self.picard(conf,appk,filetype,outsuffix,False)

    ########################
    #Picard
    ########################
    def CollectRnaSeqMetrics_jar(self,conf):
        """
        REF_FLAT=File    Gene annotations in refFlat form. Format described here: http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat Required.
        """
        appk = None
        filetype='BAM'
        outsuffix = 'RnaSeqMetrics'
        #no changes on existing products
        self.picard(conf,appk,filetype,outsuffix,False)
        
    ########################
    #Picard
    ########################
    def ReorderSam_jar(self,conf):
        
        #appk = 'order'
        #filetype='SAM,BAM'
        #outsuffix = 'bam'
        #self.picard(conf,appk,filetype,outsuffix)

        appk = 'order'
        filetype='BAM,SAM'
        outsuffix = 'bam'
        m = self._get_app(conf,filetype)
        app = m[0]
        appn = os.path.basename(app.split(' ')[0])

        suffix = '.%s.%s' % (appk,outsuffix)
                        
        current = []
        for p in self.products:
            tmp = []
            for i in p:
                px = re.compile(suffix)
                if not px.findall(i): #no need to sort again
                    fout = None
                    if i.endswith('.sam'):
                        fbase = i.rstrip('.sam')
                    elif i.endswith('.sam.gz'):
                        fbase = i.rstrip('.sam.gz')
                    elif i.endswith('.bam'):
                        fbase = i.rstrip('.bam')
                    else:
                        fbase = 'tmp'
                    
                    fout = fbase+suffix
                    
                    self.commands.append('%s INPUT=%s OUTPUT=%s REFERENCE=%s.fasta' % (app,i,fout,self.genome))
                    if self.delete_intermediate_file:
                        if appk: #BuildBamIndex
                            self.commands.append('rm -f %s' % fout.rstrip(suffix)+'.ba*')
                    tmp.append(fout)
                else: #no need to do again. SortBam will not work on "sort.bam" 
                    tmp.append(i)
            current.append(tmp)
        self.products = current
        
        
    ########################
    #Picard
    ########################
    def CleanSam_jar(self,conf):
        appk = 'clean'
        filetype='SAM'
        outsuffix = 'sam'
        self.picard(conf,appk,filetype,outsuffix)

    ########################
    #GATK
    ########################

    def GenomeAnalysisTK_jar(self,conf):
        
        #m = self._get_app(conf,'FASTQ')
        m = self._get_app(conf,None)
        app = m[0]
        t = re.compile('^.+?\s+-T\s+(.+?)\s.*').findall(app)
        if not t:
            raise Exception('-T AppName is required to run GATK')
        tx = t[0]
        
        current = []
        if tx=='UnifiedGenotyper':
            t = self.get_file_type()
            if not t == 'BAM':
                raise Exception('BAM file is required to run GATK-UnifiedGenotyper')
            #reserved_params = ['-I','-o','--dbsnp','-glm']
            reserved_params = ['-I','-o','-R']
            if self.control_case:
                for p in self.products:
                    tmp = []
                    _prefix = os.path.commonprefix(p).strip('_').strip()
                    if _prefix:
                        fvcf = '%s.vcf' % _prefix
                    else:
                        fvcf = 'tmp.vcf'
                    cmd_params = '-R %s.fasta -I %s -o %s' % (self.genome,' -I '.join(p),fvcf)
                    user_params = self.sys_params.get('-UnifiedGenotyper',None)
                    if user_params:
                        cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                    self.commands.append('%s %s' % (app,cmd_params))
                    tmp.append(fvcf)
                    current.append(tmp)
            else:
                for p in self.products:
                    tmp = []
                    for fbam in p:
                        #fvcf = '%s' % fbam.split('.')[0]+'.vcf'
                        fvcf = fbam.rstrip('.bam')+'.vcf'
                        #cmd_params = '-I %s -o %s --dbsnp %s' % (fbam,fvcf,self.dbsnp)
                        cmd_params = '-R %s.fasta -I %s -o %s' % (self.genome,fbam,fvcf)
                        user_params = self.sys_params.get('-UnifiedGenotyper',None)
                        if user_params:
                            cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                        self.commands.append('%s %s' % (app,cmd_params))
                        tmp.append(fvcf)
                    current.append(tmp)
            self.products = current

        if tx=='HaplotypeCaller':
            t = self.get_file_type()
            if not t == 'BAM':
                raise Exception('BAM file is required to run GATK-HaplotypeCaller')
            #reserved_params = ['-I','-o','--dbsnp','-glm']
            reserved_params = ['-I','-o','-R']

            """
            if self.gatk_call_variants_with_multiple_bams:
                for p in self.products:
                    tmp = []
                    _prefix = os.path.commonprefix(p).strip('_').strip()
                    if _prefix:
                        fvcf = '%s.vcf' % _prefix
                    else:
                        fvcf = 'tmp.vcf'
                
                    #cmd_params = '-I %s -o %s --dbsnp %s' % (' -I '.join(self.products),fn_vcf_output,self.dbsnp)
                    cmd_params = '-R %s.fasta -I %s -o %s' % (self.genome,' -I '.join(p),fvcf)
                    user_params = self.sys_params.get('-HaplotypeCaller',None)
                    if user_params:
                        cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                    self.commands.append('%s %s' % (app,cmd_params))
                    tmp.append(fvcf)
                    current.append(tmp)
            """
            if self.control_case:
                for p in self.products:
                    tmp = []
                    _prefix = os.path.commonprefix(p).strip('_').strip()
                    if _prefix:
                        fvcf = '%s.vcf' % _prefix
                    else:
                        fvcf = 'tmp.vcf'
                    cmd_params = '-R %s.fasta -I %s -o %s' % (self.genome,' -I '.join(p),fvcf)
                    user_params = self.sys_params.get('-HaplotypeCaller',None)
                    if user_params:
                        cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                    self.commands.append('%s %s' % (app,cmd_params))
                    tmp.append(fvcf)
                    current.append(tmp)
            
            else:
                for p in self.products:
                    tmp = []
                    for fbam in p:
                        #fvcf = '%s' % fbam.split('.')[0]+'.vcf'
                        fvcf = fbam.rstrip('.bam')+'.vcf'
                        #cmd_params = '-I %s -o %s --dbsnp %s' % (fbam,fvcf,self.dbsnp)
                        cmd_params = '-R %s.fasta -I %s -o %s' % (self.genome,fbam,fvcf)
                        user_params = self.sys_params.get('-HaplotypeCaller',None)
                        if user_params:
                            cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                        self.commands.append('%s %s' % (app,cmd_params))
                        tmp.append(fvcf)
                    current.append(tmp)
            self.products = current
        
        if tx== 'RealignerTargetCreator':
            #-T RealignerTargetCreator -R ucsc.hg19.fasta -o hg19.1000G_biallelic.indels.intervals --known 1000G_biallelic.indels.hg19.vcf
            pass
        
        if tx=='IndelRealigner':
            t = self.get_file_type()
            if not t == 'BAM':
                raise Exception('BAM file is required to run GATK-IndelRealigner')
            px = re.compile(".realign.")
            reserved_params = ['-I','-o','-R']
            for p in self.products:
                tmp = []
                for fbam in p:
                    if not px.findall(fbam):
                        fbam_output = fbam.rstrip('bam')+'realign.bam'
                        cmd_params = '-R %s.fasta -targetIntervals %s -I %s -o %s' % (self.genome,self.KNOWN_INDEL_INTERVALS,fbam,fbam_output)
                        self.commands.append('%s %s' % (app,cmd_params))
                    
                        if self.delete_intermediate_file:
                            self.commands.append('rm -f %s' % fbam.rstrip('.bam')+'.ba*')
                        tmp.append(fbam_output)
                    else:
                        tmp.append(fbam_input)
                current.append(tmp)
            self.products = current
        
        if tx=='BaseRecalibrator':
            for p in self.products:
                for fbam in p:
                    fbase = fbam.rstrip('bam')
                    foutput = fbase+'grp'
                    #app = 'GenomeAnalysisTK.jar -T IndelRealigner -R %s -targetIntervals %s -known %s' % (self.genome_fasta,self.hg19_indel_interval,self.hg19_known_indel_vcf)
                    #params = '-R %s.fasta -knownSites %s -I %s -o %s' % (self.genome,self.dbSNP,fbam,foutput)
                    params = '-R %s.fasta -I %s -o %s' % (self.genome,fbam,foutput)
                    self.commands.append('%s %s' % (app,params))
        
        if tx=='PrintReads':
            current = []
            px = re.compile(".recal.")
            for p in self.products:
                tmp = []
                for fbam in p:
                    if not px.findall(fbam):
                        fbase = fbam.rstrip('.bam')
                        foutput = fbase+'.recal.bam'
                        params = '-R %s.fasta -I %s -BQSR %s.grp -o %s' % (self.genome,fbam,fbase,foutput) 
                        self.commands.append('%s %s' % (app,params))
                        if self.delete_intermediate_file:
                            self.commands.append('rm -f %s' % fbam.rstrip('.bam')+'.ba*')
                            self.commands.append('rm -f %s' % fbam.rstrip('.bam')+'.grp')
                        tmp.append(foutput)
                    else:
                        tmp.append(fbam)
                current.append(tmp)
            self.products = current
        
        if tx=='VariantRecalibrator':
            t = self.get_file_type()
            if not t == 'VCF':
                raise Exception('VCF file is required to run GATK-VariantRecalibrator')
            
            for p in self.products:
                tmp = []
                for fvcf in p: 
                    frecal = fvcf.rstrip('.vcf')+'.recal'
                    ftranches = fvcf.rstrip('.vcf')+'.tranches'
                    #cmd_params = '-resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s.hapmap.vcf \
                    #-resource:omni,known=false,training=true,truth=false,prior=12.0 %s.omni.vcf \
                    #-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 %s.dbsnp.vcf \
                    #-an QD -an DP -an ReadPosRankSum -an MQRankSum -an FS -an MQ --maxGaussians 4 \
                    #-input %s -recalFile %s -tranchesFile %s' % (self.reference_genome,self.reference_genome,self.reference_genome,fvcf,frecal,ftranches)
                    cmd_params = '-R %s.fasta -input %s -recalFile %s -tranchesFile %s' % (self.genome,fvcf,frecal,ftranches)
                    reserved_params = ['-input','-recalFile','-R','-tranchesFile']
                    user_params = self.sys_params.get('-VariantRecalibrator',None)
                    if user_params:
                        cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                    self.commands.append('%s %s' % (app,cmd_params))
                    tmp.append(fvcf)
                current.append(tmp)
            self.products = current


        if tx=='ApplyRecalibration':
            ts_filter_level = '99.0'
            
            for p in self.products:
                tmp = []
                for f in p:
                    if f.endswith('.vcf'):
                        frecal = f.rstrip('.vcf')+'.recal'
                        ftranches = f.rstrip('.vcf')+'.tranches'
                        fout = f.rstrip('.vcf')+'.tmp.vcf'
                        reserved_params = ['-R','-input','-o','-recalFile','-tranchesFile']
                        cmd_params = '-recalFile %s -tranchesFile %s --ts_filter_level %s -input %s -o %s' % (frecal,ftranches,ts_filter_level,f,fout)
                        user_params = self.sys_params.get('-ApplyRecalibration',None)
                        if user_params:
                            cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                        self.commands.append('%s %s' % (app,cmd_params))
                        if self.delete_intermediate_file:
                            self.commands.append('rm -f %s' % frecal)
                            self.commands.append('rm -f %s' % ftranches)
                        tmp.append(fout)
                current.append(tmp)
            self.products = current   
                
        if tx=='SelectVariants':
            for p in self.products:
                tmp = []
                for f in p:
                    fout = f.rstrip('.tmp.vcf')+'.recal.vcf'
                    reserved_params = ['-R','--variant','-o']
                    cmd_params = '--variant %s -o %s' % (f,fout)
                    user_params = self.sys_params.get('-SelectVariants',None)
                    if user_params:
                        cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                    self.commands.append('%s %s' % (app,cmd_params))
                    if self.delete_intermediate_file:
                        self.commands.append('rm -f %s' % f)
                    tmp.append(fout)
                current.append(tmp)
            self.products = current
               
        if tx=='DepthOfCoverage':
            for p in self.products:
                for f in p:
                    fout = f.rstrip('bam')+'coverage'
                    params = '-R %s.fasta -I %s -o %s' % (self.genome,f,fout)
                    self.commands.append('%s %s' % (app,params))
        

        if tx=='VariantFiltration':
            """    
            $$$
            #For SNPs
            java -jar /bluegill/app/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /bluegill/data/hg19.fasta --variant x1.bwa.gatk.snp.vcf -o x1.bwa.gatk.snp.filter.vcf -cluster 3 -window 10 \
            --filterExpression "MQ < 40.0" --filterName "MQ<40.0" \
            --filterExpression "QD < 2.0" --filterName "QD<2.0" \
            --filterExpression  "FS > 60.0" --filterName "FS>60.0" \
            --filterExpression  "HaplotypeScore > 13.0" --filterName "HaplotypeScore>13.0" \
            -filterExpression  "MQRankSum < -12.5" --filterName "MQRankSum<-12.5" \
            --filterExpression  "ReadPosRankSum < -8.0" --filterName "ReadPosRankSum<-8.0"
        
            #For Indels
            java -jar /bluegill/app/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /bluegill/data/hg19.fasta --variant x1.bwa.gatk.indel.vcf -o x1.bwa.gatk.indel.filter.vcf -cluster 3 -window 10 \
            --filterExpression "QD < 2.0" --filterName "QD<2.0" \
            --filterExpression "FS > 200.0" --filterName "FS>200.0" \
            --filterExpression "InbreedingCoeff < -0.8" --filterName "InbreedingCoeff<-0.8" \
            --filterExpression "ReadPosRankSum < -20.0" --filterName "ReadPosRankSum<-20.0"

            AA ancestral allele
            AC allele count in genotypes, for each ALT allele, in the same order as listed
            AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
            AN total number of alleles in called genotypes
            BQ RMS base quality at this position
            CIGAR cigar string describing how to align an alternate allele to the reference allele
            DB dbSNP membership
            DP combined depth across samples, e.g. DP=154
            END end position of the variant described in this record (esp. for CNVs)
            H2 membership in hapmap2
            MQ RMS mapping quality, e.g. MQ=52
            MQ0 Number of MAPQ == 0 reads covering this record
            NS Number of samples with data
            SB strand bias at this position
            SOMATIC indicates that the record is a somatic mutation, for cancer genomics
            VALIDATED validated by follow-up experiment
        
        
        
            MQ: RMS mapping quality, e.g. MQ=52
        
            HaplotypeScore: Consistency of the site with two (and only two) segregating haplotypes.
                        Higher scores are indicative of regions with bad alignments, often leading to artifactual SNP and indel calls.

 
            QD(QualByDepth): Variant confidence (given as (AB+BB)/AA from the PLs) / unfiltered depth.
                    Low scores are indicative of false positive calls and artifacts.
                        
        
            FS(FisherStrand): Phred-scaled p-value using Fisher's Exact Test to detect strand bias in the reads.
                    the variation being seen on only the forward or only the reverse strand.  More bias is indicative of false positive calls.
        
        
            HRun(HomopolymerRun): Largest contiguous homopolymer run of the variant allele in either direction on the reference.
        
            VQSLOD: Only present when using Variant quality score recalibration. 
                    Log odds ratio of being a true variant versus being false under the trained gaussian mixture model.
        
            "QD < 2.0", "MQ < 40.0", "FS > 60.0", "HaplotypeScore > 13.0", "MQRankSum < -12.5", "ReadPosRankSum < -8.0".
        
            """
            for p in self.products:
                tmp = []
                for fin in p:
                    self.commands.append('vcf2snp.py %s' % fin)
                    vcf_snp = fin.rstrip('vcf')+'snp.vcf' 
                    vcf_indel = fin.rstrip('vcf')+'indel.vcf'

                    vcf_snp_filter = vcf_snp.rstrip('vcf')+'filter.vcf' 
                    vcf_indel_filter = vcf_indel.rstrip('vcf')+'filter.vcf' 
            
                    #-cluster 3 -window 10 --filterName "SNP" --filterExpression "QD<2.0||MQ<40.0||FS>60.0||HaplotypeScore>13.0" 
                    params_snp = '--variant %s -o %s ' % (vcf_snp,vcf_snp_filter)
                    self.commands.append('%s %s' % (app,params_snp))
                    vcf_snp_pass = vcf_snp_filter.rstrip('vcf')+'pass.vcf' 
                    self.commands.append("""awk '/^#/ {print $0; next} $7 == "PASS" {print $0}' %s  > %s """ % (vcf_snp_filter,vcf_snp_pass))


                    #--filterName "INDEL" --filterExpression "QD<2.0||FS>200.0"
                    params_indel = '--variant %s -o %s  ' % (vcf_indel,vcf_indel_filter)
                    self.commands.append('%s %s' % (app,params_indel))
                    vcf_indel_pass = vcf_indel_filter.rstrip('vcf')+'pass.vcf' 
                    self.commands.append("""awk '/^#/ {print $0; next} $7 == "PASS" {print $0}' %s  > %s """ % (vcf_indel_filter,vcf_indel_pass))
                    
            
                    vcf_filter = fin.rstrip('vcf')+'filter.vcf' 
                    self.commands.append('cat %s %s > %s' % (vcf_snp_pass,vcf_indel_pass,vcf_filter))
                    self.commands.append('rm -f %s' % (' '.join([vcf_snp,vcf_indel,vcf_snp_filter,vcf_indel_filter,vcf_snp_pass,vcf_indel_pass])))
                    
                    tmp.append(vcf_filter)
                current.append(tmp)
            self.products = current

        if tx=='ClipReads':
            #-T ClipReads -I my.bam -I your.bam -o my_and_your.clipped.bam -R Homo_sapiens_assembly18.fasta \
            #-XF seqsToClip.fasta -X CCCCC -CT "1-5,11-15" -QT 10
            current = []
            px = re.compile(".clip.")
            for p in self.products:
                tmp=[]
                for fbam in p:
                    if not px.findall(fbam):
                        fout = fbam.rstrip('.bam')+'.clip.bam'
                        params = '-I %s -o %s' % (fbam,fout)  
                        self.commands.append('%s %s' % (app,params))
                        if self.delete_intermediate_file:
                            self.commands.append('rm -f %s' % fbam)
                        tmp.append(fout)
                    else:
                        tmp.append(fbam)
                current.append(tmp)       
            self.products = current

        if tx=='ReduceReads':
            #http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_compression_reducereads_ReduceReads.html
            #-T ReduceReads -R ref.fasta -I myData.bam -o myData.reduced.bam
            current = []
            px = re.compile(".reduce.")
            for p in self.products:
                tmp=[]
                for fbam in p:
                    if not px.findall(fbam):
                        fout = fbam.rstrip('.bam')+'.reduce.bam'
                        params = '-R %s.fasta -I %s -o %s' % (self.genome,fbam,fout)  
                        self.commands.append('%s %s' % (app,params))
                        if self.delete_intermediate_file:
                            self.commands.append('rm -f %s' % fbam)
                        tmp.append(fout)
                    else:
                        tmp.append(fbam)
                current.append(tmp)       
            self.products = current

        if tx=='VariantEval':
            #app = 'GenomeAnalysisTK.jar  -T VariantEval -EV TiTvVariantEvaluator -R %s.fasta -D %s -eval %s -l INFO' % (self.genome,self.dbsnp,fvcf)
            t = self.get_file_type()
            if not t == 'VCF':
                raise Exception('VCF file is required to run GATK-VariantEval')

            for p in self.products:
                for fvcf in p:
                    #app = 'GenomeAnalysisTK.jar  -T VariantEval -EV TiTvVariantEvaluator -R %s.fasta -D %s -eval %s -l INFO' % (self.genome,self.dbsnp,fvcf)
                    ks = ('TiTvVariantEvaluator','CountVariants','VariantQualityScore')
                    for k in ks: 
                        cmd_params = '-EV %s -R %s.fasta -l INFO -eval %s -o %s' % (k,self.genome,fvcf,fvcf.rstrip('vcf')+k)
                        self.commands.append('%s %s' % (app,cmd_params))
            
    ########################
    #VAAST
    ########################
    def vaast_converter(self,conf):
        m = self._get_app(conf,'VCF')
        app = m[0]
        current = []
        for p in self.products:
            tmp = []
            for fvcf in p:
                output = fvcf.rstrip('vcf')+'gvf'
                cmd = '%s %s -s -b %s > %s' % (app,fvcf,self.genome,output)
                self.commands.append(cmd)
                tmp.append(output)
            current.append(tmp)
        self.products = current 

    ########################
    #VAAST
    ########################
    def VAT(self,conf):
        gene_reference = '%s.refgene.gff3' % self.genome
        m = self._get_app(conf,'GVF')
        app = m[0]
        current = []
        for p in self.products:
            tmp = []
            for fgvf in p:
                output = fgvf.rstrip('gvf')+'vat.gvf'
                cmd = '%s -f %s -a %s.fasta  %s > %s' % (app,gene_reference,self.genome,fgvf,output)
                self.commands.append(cmd)
                tmp.append(output)
            current.append(tmp)
        self.products = current 

    ########################
    #VAAST
    ########################
    def VST(self,conf):
        """
        -o 'I(U(0..3),C(4..5))'

        Set operations to perform on GVF files. Files within operations
        are referenced by their 0-based index into the list of files
        given.  Operations can be nested arbitrarily deep.  The operations
        must be quoted or escaped.
        
      U - Union: All variants.
      I - Intersection: Variants shared by all files.
      C - Compliment: Variants unique to the first file.
      D - Difference: Variants unique to any file.
      S - Shared: Variants shared by n files. S(n,0..2);
          The value for n can be a positive integer, in which case all
          variants present in at least n files will be retained.  In
          addition, the value of n can be an integer preceded by one
          of the following comparison operators:
            = - Exactly n files share the variant.
            > - Greater than n files share the variant.
            < - Less than n files share the variant.
        """
        #appname,apppath,genome,appparams = conf[0]
        #print conf[0]
        #conf = [('VarScan.jar mpileup2snp', '/bluegill/version/app/varscan/2.3.3/', 'hg19', '--output-vcf 1')]
        #current = []
        #mainapp,subapp = appname.split(' ')

        m = self._get_app(conf,'GVF')
        app = m[0]
        
        current = []

        
        if not self.control_case:
            raise Exception('Application [VST] require both control and case samples')
        
        #print self.products
        
        _control,_case = self.products
        #print 'PRODUCTS:',_control,_case
        #_normal = _control[0] 
        #_tumor = _case[0]
        
        tmp = []
        fout = os.path.commonprefix(_control).strip('_')
        if not fout:
            fout = 'tmp.'+str(time.time())
            time.sleep(1)
        foutput = fout+'.cdr'
        
        if len(_control) >1:
            cmd = "%s  -o 'I(0..%d)' -b %s %s > %s" % (app,len(_control)-1,self.genome,' '.join(_control),foutput)
        else: 
            cmd = "%s  -o 'I(0)' -b %s %s > %s" % (app,self.genome,' '.join(_control),foutput)   
        self.commands.append(cmd) 
        tmp.append(foutput)
        current.append(tmp)
                     
        tmp = []
        fout = os.path.commonprefix(_case).strip('_')
        if not fout:
            fout = 'tmp.'+str(time.time())
            time.sleep(1)
        foutput = fout+'.cdr'
        if len(_case) >1:
            cmd = "%s  -o 'I(0..%d)' -b %s %s > %s" % (app,len(_case)-1,self.genome,' '.join(_case),foutput)
        else: 
            cmd = "%s  -o 'I(0)' -b %s %s > %s" % (app,self.genome,' '.join(_case),foutput)   
        self.commands.append(cmd) 
        tmp.append(foutput)
        current.append(tmp)
        
        self.products = current 
        
    ########################
    #VAAST
    ########################
    def VAAST(self,conf):
        if not self.control_case:
            raise Exception('Application [VAAST] require both control and case samples')
        
        gene_reference = '%s.refgene.gff3' % self.genome
        
        current = []
        tmp=[]
        _control, _case = self.products
        bg_cdr =  _control[0]
        tg_cdr =  _case[0]
        output = os.path.commonprefix((bg_cdr,tg_cdr)).rstrip('_').rstrip('.cdr')        
        if not output:
            output = 'tmp'
        foutput = output+'.vaast'
        m = self._get_app(conf,'CDR')
        app = m[0]
        #cmd_params = "-iht r -lh n -fast_gp -d 1e4 -r 0.00035 -m lrt -k %s -o %s %s %s 2>%s.log" % (self.genome,fn_vaast_output,bg_cdr,tg_cdr,output)
        cmd_params = "-m lrt -o %s %s %s %s 2>vaast.log" % (foutput,gene_reference,bg_cdr,tg_cdr)
        self.commands.append('%s %s' % (app,cmd_params))
        tmp.append(foutput)
        current.append(tmp)
        self.products = current 
        #print self.products
    ########################
    #ANNOVAR
    ########################
    def convert2annovar_pl(self,conf):
        #convert2annovar.pl %s -format vcf4 -includeinfo > %s
        m = self._get_app(conf,'VCF')
        app = m[0]
        current = []
        for p in self.products:
            tmp = []
            for fvcf in p:
                fout = fvcf.rstrip('.vcf')+'.annovar'
                self.commands.append('%s %s > %s' % (app,fvcf,fout))
                tmp.append(fout)
            current.append(tmp)
        self.products = current   

    ########################
    #ANNOVAR
    ########################
    def annotate_variation_pl(self,conf):
        #annotate_variation.pl:/bluegill/app/annovar/0.1/::COMMENT:-geneanno -dbtype refseq
        #
        m = self._get_app(conf,'ANNOVAR')
        app = m[0]
        gene_model = 'refgene'
        pt = re.compile('^.+?\s+-dbtype\s+(\w+)')
        m = pt.findall(app)
        if m:
            gene_model = m[0].lower()
            
        fn_gene = ''
        fn_seq = ''
        fn_ref = ''
        path_libfiles='.'
        
        if gene_model=='refgene':
            fn_gene = '%s_refGene.txt' % self.genome
            fn_seq = '%s_refGeneMrna.fa' % self.genome
            fn_ref = 'refLink.txt.gz'
            
        elif gene_model=='knowngene':
            fn_gene = '%s_knownGene.txt' % self.genome
            fn_seq = '%s_knownGeneMrna.fa' % self.genome
            fn_ref = '%s_kgXref.txt' % self.genome
            
        elif gene_model=='ensgene':
            fn_gene = '%s_ensGene.txt' % self.genome
            fn_seq = '%s_ensGeneMrna.fa' % self.genome
            #fn_ref = ''
            
        else:
            raise Exception('Unsupported gene model [%s], available options are [refGene,knownGene or ensGene]' % gene_model)
        
        self.commands.append('cp -L %s .' % os.path.join(self.util.CLUSTER_DATA_DIR,fn_gene))
        
        self.commands.append('cp -L %s .' % os.path.join(self.util.CLUSTER_DATA_DIR,fn_seq))
        
        if fn_ref:
            self.commands.append('cp -L %s .' % os.path.join(self.util.CLUSTER_DATA_DIR,fn_ref))
            
        current = []
        for p in self.products:
            tmp = []
            for fannovar in p: 
                cmd_params = '--buildver %s' % self.genome 
                user_params = self.sys_params.get('-annovar',None)
                if user_params:
                    cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                self.commands.append('%s %s %s .' % (app,cmd_params,fannovar))
                tmp.append(fannovar)
                
            if self.delete_intermediate_file:
                self.commands.append('rm -fr %s %s %s' % (fn_gene,fn_seq,fn_ref))
            current.append(tmp)   
        self.products = current   

    ########################
    #samtools
    ########################
    def samtools(self,conf):
        
        appname,apppath,genome,appparams = conf[0]
        mainapp,subapp = appname.split(' ')
        
        if subapp=='mpileup':
            current = []
            m = self._get_app(conf,'BAM')
            app = m[0]
            cmd_params = '-f %s.fasta' % self.genome
            for p in self.products:
                tmp = []
                pname = os.path.commonprefix(p)
                if not pname:
                    pname = 'samtools.tmp'
                else:
                    pname = pname.strip('_')
                    #raise Exception('No Common Name')
                pname = '%s.%s' % (pname,subapp)
                self.commands.append('%s %s %s > %s' % (app,cmd_params,' '.join(p),pname))
                tmp.append(pname)
                current.append(tmp)
            self.products = current
               
        if subapp=='merge':
            current = []
            m = self._get_app(conf,'BAM')
            app = m[0]
            for p in self.products:
                tmp = []
                pname = os.path.commonprefix(p)
                if not pname:
                    raise Exception('No Common Name')
                fout = '%s.bam' % pname
                self.commands.append('%s %s %s' % (app,fout,' '.join(p)))
                tmp.append(fout)
                current.append(tmp)
            self.products = current   
        
    ########################
    #USeq
    ########################
    def MergePairedSamAlignments_jar(self,conf):
        """
        
        MergePairedSamAlignments.jar -f /Novo/Run7/ -c -s /Novo/STPParsedBams/run7.bam -d 10000
        
        -f The full path file or directory containing raw xxx.sam(.gz/.zip OK)/.bam file(s)
          paired alignments that are sorted by query name (standard novoalign output).
          Multiple files will be merged.

        Default Options:
            -s Save file, defaults to that inferred by -f. If an xxx.sam extension is provided,
                  the alignments won't be sorted by coordinate and saved as a bam file.
            -a Maximum alignment score (AS:i: tag). Defaults to 120, smaller numbers are more
                  stringent. Approx 30pts per mismatch for novoalignments.
            -q Minimum mapping quality score, defaults to 13, larger numbers are more stringent.
                  Set to 0 if processing splice junction indexed RNASeq data.
            -r The second paired alignment's strand is reversed. Defaults to not reversed.
            -d Maximum acceptible base pair distance for merging, defaults to 5000.
            -m Don't cross check read mate coordinates, needed for merging repeat matches. Defaults
                  to checking.
            -o Merge all proper paired alignments. Defaults to only merging those that overlap.
            -k Skip merging paired alignments. Defaults to merging. Useful for testing effect of
                  merging on downstream analysis.
        
        java -Xmx2048M -jar /bluegill/version/app/useq/8.5.0/Apps/MergePairedSamAlignments.jar -f ../test3/X1_3.sam.gz -s 2.bam -d 10000
         
         
        If input is a BAM, make sure it was sorted by "queryname", not "coordinate"  
        Error, your bam file appears sorted by coordinate. Sort by query name and restart.                                   
        """ 
        
        m = self._get_app(conf,'BAM,SAM')
        app = m[0]
        current = []
        for p in self.products:
            tmp = []
            for f in p:
                fout='%s.mpsa.bam' % f 
                self.commands.append('%s -f %s -s %s' %  (app,f,fout))
                tmp.append(fout)
            current.append(tmp)   
        self.products = current   

    ########################
    #USeq
    ########################
    def Sam2USeq_jar(self,conf):
        """
        Example: java -Xmx1500M -jar pathTo/USeq/Apps/Sam2USeq -f /Data/SamFiles/ -r
             -v H_sapiens_Feb_2009 -b ccdsExons.bed.gz
             
              /home/Annotation/Human/Hg19/CCDSExonsRptMskedMerged.bed.zip
              
        Required Options:
        -f Full path to a bam or a sam file (xxx.sam(.gz/.zip OK) or xxx.bam) or directory
              containing such. Multiple files are merged.
        -v Versioned Genome (ie H_sapiens_Mar_2006, D_rerio_Jul_2010), see UCSC Browser,
              http://genome.ucsc.edu/FAQ/FAQreleases.
        """
        m = self._get_app(conf,'BAM,SAM')
        app = m[0]
        #current = []  
        for p in self.products:
            tmp = []
            for f in p:
                #fout = 
                self.commands.append('%s -f %s' %  (app,f))
                #tmp.append(fsnp)
                #tmp.append(findel)
            #current.append(tmp)  
        #self.products = current
    
    def RNASeq_jar(self,conf):
        """ UCSC RefFlat/RefSeq gene table
        cp /home/Genomes/Mouse/Mm10/GeneAnnotations/mm10EnsGenes_JustNormal.ucsc.gz  mm10.EnsGenes_JustNormal.ucsc.gz        
        
        @align -novoalign [-o SAM -r All 50 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ] -g mm10EnsTransRad46bpNum100kMin10SplicesChrNormPhiXAdaptr -i *txt.gz

        #Delete the fastq data to create space on the cluster node
        rm *txt.gz

        #Move alignments to appropriate directories
        mkdir T C
        mkdir T/1 T/2 T/3 C/1 C/2 C/3
        mv 9434X1* C/1/
        mv 9434X2* C/2/
        mv 9434X3* C/3/
        mv 9434X4* T/1/
        mv 9434X5* T/2/
        mv 9434X6* T/3/

        #Launch RNASeq app
        java -Xmx2G -jar pathTo/USeq/Apps/RNASeq -v D_rerio_Dec_2008 -t /Data/PolIIMut/ -c /Data/PolIIWT/ -s /Data/Results/MutVsWT -g /Anno/zv8Genes.ucsc
        
The RNASeq application is a wrapper for processing RNA-Seq data through a variety of
USeq applications. It uses the DESeq package for calling significant differential
expression.  3-4 biological replicas per condition are strongly recommended. See 
http://useq.sourceforge.net/usageRNASeq.html for details constructing splice indexes,
aligning your reads, and building a proper gene (NOT transcript) table.

The pipeline:
   1) Converts raw sam alignments containing splice junction coordinates into genome
          coordinates outputting sorted bam alignemnts.
   2) Makes relative read depth coverage tracks.
   3) Scores known genes for differential exonic and intronic expression using DESeq
         and alternative splicing with a chi-square test.
   4) Identifies unannotated differentially expressed transfrags using a window
         scan and DESeq.

Use this application as a starting point in your transcriptome analysis.

Options:
-s Save directory, full path.
-t Treatment alignment file directory, full path.  Contained within should be one
       directory per biological replica, each containing one or more raw
       SAM (.gz/.zip OK) files.
-c Control alignment file directory, ditto.  
-n Data is stranded. Only analyze reads from the same strand as the annotation.
-v Genome version (e.g. H_sapiens_Feb_2009, M_musculus_Jul_2007), see UCSC FAQ,
      http://genome.ucsc.edu/FAQ/FAQreleases.
-g UCSC RefFlat or RefSeq gene table file, full path. Tab delimited, see RefSeq Genes
       http://genome.ucsc.edu/cgi-bin/hgTables, (uniqueName1 name2(optional) chrom
       strand txStart txEnd cdsStart cdsEnd exonCount (commaDelimited)exonStarts
       (commaDelimited)exonEnds). Example: ENSG00000183888 C1orf64 chr1 + 16203317
       16207889 16203385 16205428 2 16203317,16205000 16203467,16207889 . NOTE:
       this table should contain only ONE composite transcript per gene (e.g. use
       Ensembl genes NOT transcripts). Use the MergeUCSCGeneTable app to collapse
       transcripts to genes. See the RNASeq usage guide for details.
-r Full path to R, defaults to '/usr/bin/R'. Be sure to install Ander's DESeq
       (http://www-huber.embl.de/users/anders/DESeq/) R library.

Advanced Options:
-m Combine replicas and run single replica analysis using binomial based statistics,
       defaults to DESeq and a negative binomial test.
-a Maximum alignment score. Defaults to 120, smaller numbers are more stringent.
-d Minimum FDR threshold for filtering windows, defaults to 0.5
-o Don't delete overlapping exons from the gene table.
-e Print verbose output from each application.

Example: java -Xmx2G -jar pathTo/USeq/Apps/RNASeq -v D_rerio_Dec_2008 -t 
      /Data/PolIIMut/ -c /Data/PolIIWT/ -s
      /Data/Results/MutVsWT -g /Anno/zv8Genes.ucsc 
         
        """
        
        #appname,apppath,genome,appparams = conf[0]
        
        if not self.control_case:
            raise Exception('Application [RNASeq.jar] require both control and case samples')
        
        _control,_case = self.products
        #print _control,_case
        #fout = os.path.commonprefix((_control[0],_case[0])).strip('_')
        #if not fout:
        #    fout = 'tmp'

        m = self._get_app(conf,'BAM,SAM')
        app = m[0] #not a grouped app
        #gene_table = '%s.ensgenes.ucsc' % self.genome        
        #if not os.path.exists(os.path.join(self.util.CLUSTER_DATA_DIR,gene_table)):
        #    error_msg = 'The default gene table %s was not found, please choose one genetable (ends with .ucsc) from %s and use "-rnaseq [-g MY_GENE_TABLE.ucsc]"' % (gene_table,os.path.join(self.util.LOCAL_SERVER_DIR,self.util.PATH_DATA))
        #    raise Exception(error_msg)
        name_control = 'T'
        name_case = 'C'
        name_result = 'R'
        
        tx = str(time.time())
        dir_control=name_control+tx
        dir_case=name_case+tx
        dir_result = name_result+tx

        _control,_case = self.products
        #files_control = _control[0] 
        #files_case = _case[0]
        
        cmd = 'mkdir %s %s %s' % (dir_control,dir_case,dir_result)
        self.commands.append(cmd)
        
        cmd = 'mv %s %s/' % (' '.join(_control),dir_control)
        self.commands.append(cmd)
        
        cmd = 'mv %s %s/' % (' '.join(_case),dir_case)
        self.commands.append(cmd)
        
        #gv = 'M_musculus_Dec_2011'
        
        #cmd_params = "-s %s -v %s -t %s -c %s -r %s -g %s" % (dir_control,gv,dir_control,dir_case,self.util.CLUSTER_R_PATH,gene_table)
        
        cmd_params = "-t %s -c %s -s %s -r %s" % (dir_control,dir_case,dir_result,self.util.CLUSTER_R_PATH)

        #cmd_params = "-s . -v %s -t %s -c %s -r %s -g %s" % (self.gv,tt,tc,self.R_PATH,gene_table)
        #reserved_params = {'-v','-t','-c','-r'}
        #user_params = self.sys_params.get('-rnaseq',None)
        #if user_params:
        #    cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
        #print 'APP:',app
        self.commands.append('%s %s' % (app,cmd_params))

    ########################
    #VarScan
    ########################
    def VarScan_jar(self,conf):
        #print self.products
        appname,apppath,genome,appparams = conf[0]
        #print conf[0]
        #conf = [('VarScan.jar mpileup2snp', '/bluegill/version/app/varscan/2.3.3/', 'hg19', '--output-vcf 1')]
        current = []
        mainapp,subapp = appname.split(' ')
        if subapp=='pileup2snp' or subapp=='pileup2indel' or subapp=='pileup2cns':
            #m = self._get_app(conf,'PILEUP')
            xapp = '%s' % os.path.join(apppath.replace(self.util.LOCAL_SERVER_APP_DIR,self.util.CLUSTER_APP_DIR),appname)
            for p in self.products:
                tmp = []
                for fbam in p:
                    fout = fbam.rstrip('.pileup')
                    self.commands.append('%s %s %s' % (xapp,fbam,appparams))
                    tmp.append(fout)
                current.append(tmp)
            self.products = current   

        elif subapp=='mpileup2snp' or subapp=='mpileup2indel' or subapp=='mpileup2cns':
            tmp = []
            #m = self._get_app(conf,'PILEUP')
            xapp = '%s' % os.path.join(apppath.replace(self.util.LOCAL_SERVER_APP_DIR,self.util.CLUSTER_APP_DIR),appname)
            for p in self.products:
                tmp = []
                for fbam in p:
                    fout = fbam.rstrip('.mpileup')
                    self.commands.append('%s %s %s' % (xapp,fbam,appparams))
                    tmp.append(fout)
                current.append(tmp)
            self.products = current
               
        elif subapp=='somatic':
            """
            java -jar VarScan.jar somatic [normal_pileup] [tumor_pileup] [output] OPTIONS
                normal_pileup - The SAMtools pileup file for Normal
                tumor_pileup - The SAMtools pileup file for Tumor
                output - Output base name for SNP and indel output
            Two inputs are mandatory
            """
            tmp = []
            if not self.control_case:
                raise Exception('Application [VarScan.jar somatic] require both control and case samples')
            _control,_case = self.products
            normal_pileup = _control[0] 
            tumor_pileup = _case[0]
            fout = os.path.commonprefix((normal_pileup,tumor_pileup)).strip('_')
            if not fout:
                fout = 'tmp'
                #raise Exception('No Common Name')
            #m = self._get_app(conf,'MPILEUP')
            xapp = '%s' % os.path.join(apppath.replace(self.util.LOCAL_SERVER_APP_DIR,self.util.CLUSTER_APP_DIR),appname)
            self.commands.append('%s %s %s %s %s' % (xapp,_control[0],_case[0],fout,appparams))
            #tmp.append('%s.indel.vcf' % fout)
            #tmp.append('%s.snp.vcf' % fout)
            #self.commands.append('cat %s.snp.vcf %s.indel.vcf > %s.vcf ' % (fout,fout,fout))
            self.commands.append("awk '/^#/{next;} {print $0}' %s.indel.vcf >> %s.snp.vcf; mv %s.snp.vcf %s.vcf; rm -f %s.indel.vcf" % (fout,fout,fout,fout,fout))
            tmp.append('%s.vcf' % fout)
            current.append(tmp)
            self.products = current
        else:
            raise Exception('%s was not implemented yet' % subapp)

    ########################
    #utility methods
    ########################

    def mergelanes(self):
        """
        merge multiple lane BAMs into one sample BAM.
        """
        current = []
        for p in self.products:
            tmp = []
            if len(p)==1:
                tmp.append(p)
            else:
                m = {}
                for i in p:
                    j = i.split('/')
                    if len(j)>1:
                        m.setdefault(j[0],[]).append(i)
                    
                for k,v in m.items():
                    common_prefix = os.path.commonprefix(v).strip('_').strip()
                    if not common_prefix:
                        common_prefix = 'merge'
                    if common_prefix.endswith(os.sep): #only the folder is same
                        common_prefix += 'merge'
                        
                    fout = '%s.bam' % common_prefix
                    cmd = 'samtools merge %s %s' % (fout,' '.join(v))
                    self.commands.append(AppLocator(self.apppath).locate(cmd))
                    if delete_intermediate_file:
                        self.commands.append('rm -f %s' % ' '.join([i.rstrip('.bam')+'.ba*' for i in v]))                    
                    tmp.append(fout)
            current.append(tmp)  
        self.products = current 
        #self.SortSam()

    def vcf2snpindel(self):
        """
        python /bluegill/app/vcf2snp.py x1.vcf
        
        #snp
        java -jar /bluegill/app/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /bluegill/data/hg19.fasta --variant x1.snp.vcf  -o x1.snp.filter.vcf -cluster 3 -window 10 --filterName "SNP" --filterExpression "QD<2.0||MQ<40.0||FS>60.0||HaplotypeScore>13.0"
        awk '/^#/ {print $0; next} $7 == "PASS" {print $0}' x1.snp.filter.vcf > x1.snp.filter.pass.vcf
        
        #indel
        java -jar /bluegill/app/gatk/GenomeAnalysisTK.jar -T VariantFiltration -R /bluegill/data/hg19.fasta --variant x1.indel.vcf  -o x1.indel.filter.vcf --filterName "INDEL" --filterExpression "QD<2.0||FS>200.0"
        awk '/^#/ {print $0; next} $7 == "PASS" {print $0}' x1.indel.filter.vcf > x1.indel.filter.pass.vcf
        
        cat x1.snp.filter.pass.vcf x1.indel.filter.pass.vcf > x1.pass.vcf
        """
        current = []  
        for p in self.products:
            tmp = []
            for fvcf in p:
                fsnp = fvcf.rstrip('vcf')+'snp.vcf' 
                findel = fvcf.rstrip('vcf')+'indel.vcf' 
                self.commands.append('vcf2snpindel %s' %  fvcf)
                tmp.append(fsnp)
                tmp.append(findel)
            current.append(tmp)  
        self.products = current
        

if __name__=='__main__':
    version = '1'
    #vobj = VersionConf(ver)
    #print vobj.get_current_version() 
    #print vobj.get_app_path() 
    #links = vobj.get_pipeline('@align -bam')
    #print links[0]
    #print links[1]

    #pipeline_name = '@align'
    #pipeline_name = '@align -bam'
    pipeline_name = '@snpindel -g hg19'
    #pipeline_name = '@pair'
    #pipeline_name = '@snpindel -g hg19 -bwa'
    genome = 'hg19'
    platform = 'illumina'
    inputs=['A_1.txt.gz','A_2.txt.gz']
    params={}
    ax = AppManager(version,pipeline_name,inputs,genome,platform,params)
    print '\n'.join(ax.get_command())
    #print ax.get_product()


    