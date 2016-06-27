"""This module defined the applications to be executed"""

import os
import re
import inspect
import imp
import time

class App(object):            
    def __init__(self,
                 inputs=None,
                 params=None,
                 reserve_intermediate_files=False):
        """
        @param inputs: [], list of inputs 
        @param params: {}, the parameters passed to this app
        @param reserve_intermediate_files: Boolean, if true then reserver all intermediate files
        """
        self.util = imp.load_source('util','util.py')
        self.inputs=inputs
        self.params=params #dict
        self.reserve_intermediate = reserve_intermediate_files
        self.commands = []
        self.outputs=[]
        self.reserved_parameter_keys = []

        
    def get_output(self):
        return self.outputs
    
    def get_command(self):
        return self.commands

    def get_param(self,dict_params):
        p = []
        for k,v in dict_params.items():
            if v:
                p.append(str(k)+' '+str(v))
            else:
                p.append(str(k))
        return ' '.join(p)

    def __str__(self):
        return self.get_command()
    
##################################################
#            Mapping Application                 #
##################################################
class Aligner(App):
    #def __init__(self,platform,reference_genome,sequence_files,params=None,reserve_intermediate_files=False,gzip_output=False):
    def __init__(self,platform,reference_genome,sequence_files,params=None,reserve_intermediate_files=False):
        #supported mode
        App.__init__(self,sequence_files,params,reserve_intermediate_files)
        self.platform=platform
        self.genome = reference_genome
        #SAM header tags
        self.READGROUP='READGROUP1' #Read Group
        self.PLATFORM='illumina'
        self.LIBRARY='libtmp'
        self.SAMPLE_NAME='sample'
        self.CN ='INSTITUTION'
        self.sam_header='@RG\\tID:%s\\tPL:%s\\tLB:%s\\tSM:%s\\tCN:%s' % (self.READGROUP,self.PLATFORM,self.LIBRARY,self.SAMPLE_NAME,self.CN)
        #self.gzip=gzip_output
        self.paired = self.pair_files(sequence_files)
        
    def pair_files(self,files):
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

class BWA(Aligner):
    #def __init__(self,platform,reference_genome,sequence_files,params=None,reserve_intermediate_files=False,gzip_output=False):
    def __init__(self,platform,reference_genome,sequence_files,params=None,reserve_intermediate_files=False):
        #Aligner.__init__(self,platform,reference_genome,sequence_files,params,reserve_intermediate_files,gzip_output)
        Aligner.__init__(self,platform,reference_genome,sequence_files,params,reserve_intermediate_files)
        self.genome = reference_genome+'.bwa.'+platform
        #self.inputs = sequence_files
        #self.paired = self.pair_files(sequence_files)
        
    def get_command2(self):
        output='tmp'
        if not self.inputs:
            raise Exception("Can not find any input files")

        for i in self.paired:#i is []
            cmd_aln = ['','']
            cmd_sam = None
            if len(i)==1: #single-end sequence read file
                f = i[0]
                #already aligned before?
                if f.endswith('.sam') or f.endswith('.sam.gz') or f.endswith('.bam') or f.endswith('.vcf'):
                    self.outputs = self.inputs
                    return [] #no commands!
                else:
                    output = f.split('.')[0]
                    
                    if self.platform=='illumina':
                        #cmd_aln[0] = "bwa aln -t 24 -I %s %s > %s.sai"
                        cmd_aln[0] = "bwa aln -t 24 %s %s > %s.sai"
                        #cmd_sam = "bwa samse -r '%s' %s %s.sai %s> %s.sam"
                        cmd_sam = "bwa samse -r '%s' %s %s.sai %s | gzip > %s.sam.gz"
                            
                         
                    elif self.platform=='solid':
                        #cmd_aln[0] = "bwa aln -c -t 24 -I %s %s > %s.sai"
                        cmd_aln[0] = "bwa aln -c -t 24 %s %s > %s.sai"
                        #cmd_sam = "bwa samse -r '%s' %s %s.sai %s> %s.sam"
                        cmd_sam = "bwa samse -r '%s' %s %s.sai %s | gzip > %s.sam.gz"

                    else:
                        raise Exception('%s is not supported' % self.platform)
                    
                    cmd_aln[0] = cmd_aln[0] % (self.genome,f,output)
                    cmd_sam  = cmd_sam % (self.sam_header,self.genome,output,f,output)
                    
            elif len(i)==2: #pair-end sequence read file
                
                f1,f2 = i
                if f1.endswith('.sam') or f1.endswith('.bam') or f1.endswith('.vcf'):
                    self.outputs = self.inputs
                    return []
                    
                else:
                    output = os.path.commonprefix((f1,f2)).strip('_')
                    if not output:
                        output = 'tmp'
                    
                    if self.platform=='illumina':
                        #cmd_aln[0] = "bwa aln -t 24 -I %s %s > %s.1.sai"
                        #cmd_aln[1] = "bwa aln -t 24 -I %s %s > %s.2.sai"
                        cmd_aln[0] = "bwa aln -t 24 %s %s > %s.1.sai"
                        cmd_aln[1] = "bwa aln -t 24 %s %s > %s.2.sai"
                        #cmd_sam = "bwa sampe -P -r '%s' %s %s.1.sai %s.2.sai %s %s > %s.sam"
                        cmd_sam = "bwa sampe -P -r '%s' %s %s.1.sai %s.2.sai %s %s | gzip > %s.sam.gz"
                         
                    elif self.platform=='solid':
                        #cmd_aln[0] = "bwa aln -c -t 24 -I %s %s > %s.1.sai"
                        #cmd_aln[1] = "bwa aln -c -t 24 -I %s %s > %s.2.sai"
                        cmd_aln[0] = "bwa aln -c -t 24 %s %s > %s.1.sai"
                        cmd_aln[1] = "bwa aln -c -t 24 %s %s > %s.2.sai"
                        #cmd_sam = "bwa sampe -P -r '%s' %s %s.1.sai %s.2.sai %s %s > %s.sam"
                        cmd_sam = "bwa sampe -P -r '%s' %s %s.1.sai %s.2.sai %s %s | gzip > %s.sam.gz"

                    else:
                        raise Exception('%s is not supported' % self.sz_mode)
                    
                    cmd_aln[0] = cmd_aln[0] % (self.genome,f1,output)
                    cmd_aln[1] = cmd_aln[1] % (self.genome,f2,output)
                    cmd_sam  = cmd_sam % (self.sam_header,self.genome,output,output,f1,f2,output)
            else:
                raise Exception('Two many input files: %d' % len(self.inputs))
            
            if self.params:
                user_params_aln = self.params.get('-bwaaln',None)
                if user_params_aln:
                    reserved_params_aln = []
                    if cmd_aln[0]:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_aln[0],reserved_params_aln,user_params_aln[0],first_param_index=2))
                    if cmd_aln[1]:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_aln[1],reserved_params_aln,user_params_aln[0],first_param_index=2))
                else:
                    if cmd_aln[0]:
                        self.commands.append(cmd_aln[0])
                    if cmd_aln[1]:
                        self.commands.append(cmd_aln[1])
                
                if len(cmd_aln)==1: #bwa samse
                    user_params_samse = self.params.get('-bwasamse',None)
                    reserved_params_samse = []
                    if user_params_samse:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_sam,reserved_params_samse,user_params_samse[0],first_param_index=2))
                    else:
                        self.commands.append(cmd_sam)
                    self.commands.append('rm -fr *.sai')
                else:
                    user_params_sampe = self.params.get('-bwasampe',None)
                    reserved_params_sampe = []
                    if user_params_sampe:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_sam,reserved_params_sampe,user_params_sampe[0],first_param_index=2))
                    else:
                        self.commands.append(cmd_sam)
                    self.commands.append('rm -fr *.sai')
                    
            #self.outputs.append('%s.sam' % output)
            self.outputs.append('%s.sam.gz' % output)
            
            #if self.platform=='solid':
            #    self.commands.append("cscode.py %s.sam %s.tmp" %(output,output))
            #    self.commands.append("mv %s.tmp %s.sam" %(output,output))
            
        #print 'OUTPUT:',self.outputs
        return self.commands

    def get_command(self):
        """
        bwa mem
        """
        output='tmp'
        if not self.inputs:
            raise Exception("Can not find any input files")

        for i in self.paired:#i is []
            cmd_sam = None
            if len(i)==1: #single-end sequence read file
                f = i[0]
                #already aligned before?
                if f.endswith('.sam') or f.endswith('.sam.gz') or f.endswith('.bam') or f.endswith('.vcf'):
                    self.outputs = self.inputs
                    return [] #no commands!
                else:
                    output = f.split('.')[0]
                    
                    if self.platform=='illumina':
                        cmd_sam = "bwa mem -r '%s' %s %s | gzip > %s.sam.gz"
                            
                    else:
                        raise Exception('%s is not supported' % self.platform)
                    
                    cmd_sam  = cmd_sam % (self.sam_header,self.genome,f,output)
                    
            elif len(i)==2: #pair-end sequence read file
                
                f1,f2 = i
                if f1.endswith('.sam') or f1.endswith('.bam') or f1.endswith('.vcf'):
                    self.outputs = self.inputs
                    return []
                    
                else:
                    output = os.path.commonprefix((f1,f2)).strip('_')
                    if not output:
                        output = 'tmp'
                    
                    if self.platform=='illumina':
                        cmd_sam = "bwa mem -r '%s' %s %s %s | gzip > %s.sam.gz"
                         
                    else:
                        raise Exception('%s is not supported' % self.sz_mode)
                    cmd_sam  = cmd_sam % (self.sam_header,self.genome,f1,f2,output)
                    
            else:
                raise Exception('Two many input files: %d' % len(self.inputs))
            
            if self.params:
                user_params_bwamem = self.params.get('-bwamem',None)
                if user_params_bwamem:
                    reserved_params_bwamem = []
                    if cmd_aln[0]:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_sam,reserved_params_bwamem,user_params_bwamem[0],first_param_index=2))
                    if cmd_aln[1]:
                        self.commands.append('%s' % self.util.update_cmd_params(cmd_sam,reserved_params_bwamem,user_params_bwamem[0],first_param_index=2))
                else:
                    if cmd_sam:
                        self.commands.append(cmd_sam)
                
            self.outputs.append('%s.sam.gz' % output)
        return self.commands

class Bowtie2(Aligner):
    pass

class Bowtie1(Aligner):
    pass
    
class Novocraft(Aligner):
    """a command in a pipeline, it must have at least one input, and it should have at least one output. 
    1. the output will be used as input for next command, for example, SortSam's output .bam file 
    2. the output will NOT be directly used as input for next command, for example, BuildBamIndex's output .bai file,
    which is necessary to run same commands, but you do not need to explicitly specify a .bai file in the command.
    """
    #def __init__(self,platform,reference_genome,sequence_files,params=None,reserve_intermediate_files=False,gzip_output=False):
    def __init__(self,platform,reference_genome,sequence_files,params=None,reserve_intermediate_files=False):
        """platform -> "solid","illumina" or "bisulphite"
        hg19.nov.solid.nix
        hg19.nov.illumina.nix
        hg19.nov.bisulphite.nix
        """
        #Aligner.__init__(self,platform,reference_genome,sequence_files,params,reserve_intermediate_files,gzip_output)
        Aligner.__init__(self,platform,reference_genome,sequence_files,params,reserve_intermediate_files)
        self.genome = reference_genome+'.nov.'+platform+'.nix'
        #self.inputs = sequence_files

    def get_command(self):
        output='tmp'
        if not self.inputs:
            raise Exception("Can not find any input files")

        for i in self.paired:#i is []
            cc = None
            if len(i)==1: #single-end sequence read file
                f = i[0]
                #already aligned before?
                if f.endswith('.sam') or f.endswith('.sam.gz') or f.endswith('.bam') or f.endswith('.vcf'):
                    self.outputs.append(f)

                    #no commands!
                else:
                    #output = f.split('.')[0]
                    output = f.rstrip('.gz')
                    output = output.rstrip('.gzip')
                    output = output.rstrip('.txt')
                    output = output.rstrip('.fastq')
                    output = output.rstrip('.fq')
                    output = output.rstrip('.fasta')
                    output = output.rstrip('.fa')
                    
                    if self.platform=='solid':
                        #cmd_string = "novoalignCS -r None -k -o SAM $'%s' -d %s -f %s >%s.sam 2>%s.log"
                        #cmd_string = "novoalignCS -r Random -k -o SAM $'%s' -d %s -f %s >%s.sam"
                        cmd_string = "novoalignCS -k -o SAM $'%s' -d %s -f %s | gzip >%s.sam.gz"

                        
                    elif self.platform=='illumina':
                        #cmd_string = "novoalign -r None -k -o SAM $'%s' -d %s -f %s >%s.sam 2>%s.log"
                        #cmd_string = "novoalign -r Random -k -o SAM $'%s' -d %s -f %s >%s.sam"
                        cmd_string = "novoalign -k -o SAM $'%s' -d %s -f %s | gzip >%s.sam.gz"
                        
                    elif self.platform=='bisulphite':
                        #cmd_string = "novoalign -r None -k -b 4 -o SAM $'%s' -d %s -f %s >%s.sam 2>%s.log"
                        #cmd_string = "novoalign -r Random -k -b 4 -o SAM $'%s' -d %s -f %s >%s.sam"
                        #if self.gzip:
                        cmd_string = "novoalign -k -b 4 -o SAM $'%s' -d %s -f %s | gzip >%s.sam.gz"
                        
                    else:
                        raise Exception('%s is not supported' % self.platform)
                    
                    cc = cmd_string %(self.sam_header,self.genome,f,output)
                
            elif len(i)==2: #pair-end sequence read file
                
                f1,f2 = i
                if f1.endswith('.sam') or f1.endswith('.sam.gz') or f1.endswith('.bam') or f1.endswith('.vcf'):
                    self.outputs.append(f1)
                    self.outputs.append(f2)
                    
                else:
                    output = os.path.commonprefix((f1,f2)).strip('_')
                    if not output:
                        output = 'tmp'
                    
                    if self.platform=='solid':
                        #cmd_string = "novoalignCS -r None -k -o SAM $'%s' -d %s -f %s %s >%s.sam 2>%s.log"
                        #cmd_string = "novoalignCS -r Random -k -o SAM $'%s' -d %s -f %s %s >%s.sam"
                        #if self.gzip:
                        cmd_string = "novoalignCS -k -o SAM $'%s' -d %s -f %s %s | gzip >%s.sam.gz"

                    elif self.platform=='illumina':
                        #cmd_string = "novoalign -r None -k -o SAM $'%s' -d %s -f %s %s >%s.sam 2>%s.log"
                        #cmd_string = "novoalign -r Random -k -o SAM $'%s' -d %s -f %s %s >%s.sam"
                        #if self.gzip:
                        cmd_string = "novoalign -k -o SAM $'%s' -d %s -f %s %s | gzip >%s.sam.gz"

                    elif self.platform=='bisulphite':
                        #cmd_string = "novoalign -r None -k -b 4 -o SAM $'%s' -d %s -f %s >%s.sam 2>%s.log"
                        #cmd_string = "novoalign -r Random -k -b 4 -o SAM $'%s' -d %s -f %s %s >%s.sam"
                        #if self.gzip:
                        cmd_string = "novoalign -k -b 4 -o SAM $'%s' -d %s -f %s %s | gzip >%s.sam.gz"

                    else:
                        raise Exception('%s is not supported' % self.sz_mode)
                    
                    cc = cmd_string %(self.sam_header,self.genome,f1,f2,output)
                    
            else:
                raise Exception('Two many input files: %d' % len(self.inputs))
            
            if cc:
                if self.params:
                    user_params = self.params.get('-novoalign',None)
                    if user_params:
                        reserved_params = ['-d','-f']
                        cc = self.util.update_cmd_params(cc,reserved_params,user_params[0],first_param_index=1)
                        self.commands.append('%s' % cc)
                        #self.commands.append(update_cmd_params(cc,self.reserved_parameter_keys,user_params))
                    else:
                        self.commands.append(cc)
                #self.outputs.append(output+'.sam')
                self.outputs.append(output+'.sam.gz')
                #if self.platform=='solid':
                #    self.commands.append("cscode.py %s.sam %s.tmp" %(output,output))
                #    self.commands.append("mv %s.tmp %s.sam" %(output,output))
            else:
                pass
            
        return self.commands


##################################################
#            GATK                                #
##################################################

class Genotyper(App):
    def __init__(self,genome,alignment_files,params=None,reserve_intermediate_files=False):
        App.__init__(self,inputs=alignment_files,params=params,reserve_intermediate_files=reserve_intermediate_files)
        self.genome = genome

        
class GATK(Genotyper):
    """GATK v2
    """
    def __init__(self,analysis_type,genome,alignment_files,params=None,reserve_intermediate_files=False,bool_multiple_bams=False):
        """
        -CountCovariates []
        -DepthOfCoverage []
        -UnifiedGenotyper []
        -VariantFiltration []
        """
        Genotyper.__init__(self,genome,alignment_files,params,reserve_intermediate_files)
        self.analysis_type = analysis_type #recal or genotype
        #self.dbsnp = self.genome+'.dbsnp_132.vcf'
        
        my_dbsnp = self.params.get('-dbsnp',None)
        if my_dbsnp:
            self.dbsnp = my_dbsnp[0]
        else:        
            self.dbsnp = self.genome+'.dbsnp.vcf'
                    
        self.genome_fasta = self.genome+'.fasta'

        #hg19 specific files
        #self.hg19_indel_interval = 'hg19.1000G_biallelic.indels.intervals'
        
        self.hg19_omni_vcf = 'hg19.1000G_omni2.5.vcf'
        
        self.hg19_hapmap_vcf = 'hg19.hapmap_3.3.vcf'
        
        #self.hg19_known_snp_vcf = ['hg19.1000G_omni2.5.vcf','hg19.hapmap_3.3.vcf',self.dbsnp]
        
        #self.hg19_known_indel_vcf = ['hg19.1000G_omni2.5.vcf','hg19.1000G_phase1.indels.vcf']
        
        self.hg19_known_indel_vcf = 'hg19.1000G_phase1.indels.vcf'
        
        self.hg19_indel_intervals = 'hg19.indel.intervals'
        
        self.call_mode='SNP' #SNP, INDEL or BOTH
        
        #current files to be processed
        self.products = []
        self.products = self.inputs  #'A.sam' or ['A.sam','B.sam']
        #self.validation_stringency='LENIENT'
        self.validation_stringency='SILENT'
        
        self.multiple_bams = bool_multiple_bams
        
    def get_output(self):
        return self.products
        
    def get_command(self):
        """Best Practice Variant Detection with the GATK v4, for release 2.0
        
        for each sample
            lanes.bam <- merged lane.bams for sample
            dedup.bam <- MarkDuplicates(lanes.bam)
            realigned.bam <- realign(dedup.bam) [with known sites included if available]
            recal.bam <- recal(realigned.bam)
            sample.bam <- recal.bam
        """
        for fx in self.products:
            if fx.endswith('.vcf'):
                return []
                    
        #SAM->BAM
        self.FixMateInformation()

        #BAM->BAM
        self.MarkDuplicates()

        #BAM->BAM
        #self.ReduceReads() #not compatible with HaplotypeCaller
        
        #BAM->BAM
        self.SortSam()
        
        if self.analysis_type == 'dedup':
            return self.commands
         
        #BAM->BAM
        #print 'INPUTS',self.inputs
        if self.genome=='hg19': #realign around indels
            self.IndelRealigner()
            self.BuildBamIndex()
            #if len(self.inputs) > 1:
            #    if not self.multiple_bams:
                    #if type(self.inputs).__name__=='list':#merge these bam files then recalibrate
            #        self.MergeSamFiles()
            #        self.BuildBamIndex()
            self.MergeLanes()
            self.BaseRecalibrator()
            #self.TableRecalibration()
            self.SortSam()
            
        if self.analysis_type == 'recal':
            return self.commands
        
        #BAM -> BAM
        if self.params.has_key('-DepthOfCoverage'):
            self.DepthOfCoverage()
            
        #BAM -> VCF
        elif self.analysis_type=='snp':
            self.UnifiedGenotyper('SNP')
            #self.HaplotypeCaller('SNP')
            
        #BAM -> VCF
        elif self.analysis_type=='indel':
            self.UnifiedGenotyper('INDEL')
            #self.HaplotypeCaller('INDEL')
            
        #BAM -> VCF
        elif self.analysis_type=='snpindel':
            self.UnifiedGenotyper('BOTH')
            #self.HaplotypeCaller('BOTH')
            
        #BAM -> VCF
        elif self.analysis_type=='annot':
            self.UnifiedGenotyper('BOTH')
            #self.HaplotypeCaller('BOTH')
            
        else:
            raise Exception('Analysis "%s" is not supported' % self.analysis_type)

        if self.genome=='hg19': #realign around indels
            self.VariantRecalibrator()
        
        #VCF -> VCF
        #if self.params.has_key('-VariantFiltration') or self.params.has_key('-vf'):
        #    self.VariantFiltration()

        #VCF -> VCF
        #if self.params.has_key('-VariantRecalibrator') or self.params.has_key('-vr'):
        #    self.VariantRecalibrator()
         
        #VCF -> VCF
        if self.params.has_key('-VariantEval'):
            self.VariantEvalWalker()
            
        #if self.params.has_key('-TiTvVariantEvaluator') or self.params.has_key('-tt'):
        #    self.TiTvVariantEvaluator()
             
        #VCF -> VCF
        #if self.params.has_key('-CountVariants') or self.params.has_key('-cv'):
        #    self.CountVariants()
        
        #VCF -> VCF
        #if self.params.has_key('-VariantQualityScore') or self.params.has_key('-vq'):
        #    self.VariantQualityScore()
             
        return self.commands
      
    def _VarianceEval(self):
        """http://www.broadinstitute.org/gsa/wiki/index.php/VariantEval
        
        This method can NOT be used directly, instead it MUST be called in other Evaluation Modules 
        """
        for fvcf in self.products:
            myname = inspect.stack()[1][3] # get parent's name -> who is calling me?'
            freport = fvcf.rstrip('vcf')+myname
            app = 'GenomeAnalysisTK.jar  -T VariantEval -EV %s -R %s -D %s -eval %s -l INFO' % (myname,self.genome_fasta,self.dbsnp,fvcf)
            cmd_params = '-o %s' % freport
            reserved_params = []
            user_params = self.params.get('-%s' % myname,None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])            
            self.commands.append('%s %s' % (app,cmd_params))
        
    def TiTvVariantEvaluator(self):
        self._VarianceEval()
        
    def CountVariants(self):
        self._VarianceEval()
        
    def VariantQualityScore(self):
        self._VarianceEval()
    
    def VariantEvalWalker(self):
        """-R ref.fasta -T VariantEval -o output.eval.gatkreport --eval:set1 set1.vcf  --eval:set2 set2.vcf  [--comp comp.vcf]
        
        Transition/Transversion (Ti/Tv) ratio
        
        Whole Genome: 2.0-2.1
        Exome: 3.0-3.5
        
        """
        app = 'GenomeAnalysisTK.jar -T VariantEval -R %s' % (self.genome_fasta)
        _inputs = []
        for i,fvcf in enumerate(self.products):
            _inputs.append('--eval:set%s %s' % (i+1,fvcf))
        _output = 'output.eval.gatkreport'
        cmd_params = '-o %s %s' % (_output,' '.join(_inputs))
        reserved_params = []
        user_params = self.params.get('-VariantEval',None)
        if user_params:
            cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])            
        self.commands.append('%s %s' % (app,cmd_params))
        
    
    def IndelRealigner(self):
        """
        java -Xmx4g -jar ~/tool/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ucsc.hg19.fasta 
        -o hg19.1000G_biallelic.indels.intervals --known 1000G_biallelic.indels.hg19.vcf
                
        java -Xmx4g -Djava.io.tmpdir=/path/to/tmpdir \
          -jar /path/to/GenomeAnalysisTK.jar \
          -I <lane-level.bam> \
          -R <ref.fasta> \
          -T IndelRealigner \
          -targetIntervals <intervalListFromStep1Above.intervals> \
          -o <realignedBam.bam> \
          --known /path/to/indel_calls.vcf
          --consensusDeterminationModel KNOWNS_ONLY \
          -LOD 0.4
        """
        tmp = []
        for fn_bam_input in self.products:
            px = re.compile(".realign.")
            if not px.findall(fn_bam_input):
                fn_bam_output = fn_bam_input.rstrip('bam')+'realign.bam'
                #app = 'GenomeAnalysisTK.jar -T IndelRealigner -R %s -targetIntervals %s -known %s' % (self.genome_fasta,self.hg19_indel_interval,self.hg19_known_indel_vcf)
                app = 'GenomeAnalysisTK.jar -T IndelRealigner -R %s -targetIntervals %s' % (self.genome_fasta,self.hg19_indel_intervals)
                cmd_params = '-I %s -o %s' % (fn_bam_input,fn_bam_output)
                self.commands.append('%s %s' % (app,cmd_params))
                if not self.reserve_intermediate:
                    self.commands.append('rm -f %s' % fn_bam_input.rstrip('.bam')+'.ba*')
                tmp.append(fn_bam_output)
            else:
                tmp.append(fn_bam_input)
        self.products = tmp   


    def BaseRecalibrator(self):
        """ 
        -T BaseRecalibrator \
        -I my_reads.bam \
        -R resources/Homo_sapiens_assembly18.fasta \
        -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
        -knownSites another/optional/setOfSitesToMask.vcf \
        -o recal_data.grp

        -T PrintReads -R reference.fasta -I input.bam -BQSR recalibration_report.grp -o output.bam

        """
        tmp = []
        for fn_bam_input in self.products:
            px = re.compile(".recal.")
            if not px.findall(fn_bam_input):
                fbase = fn_bam_input.rstrip('bam')
                fn_grp_output = fbase+'grp'
                #app = 'GenomeAnalysisTK.jar -T IndelRealigner -R %s -targetIntervals %s -known %s' % (self.genome_fasta,self.hg19_indel_interval,self.hg19_known_indel_vcf)
                app = 'GenomeAnalysisTK.jar -T BaseRecalibrator -R %s -knownSites %s' % (self.genome_fasta,self.dbsnp)
                cmd_params = '-I %s -o %s' % (fn_bam_input,fn_grp_output)
                self.commands.append('%s %s' % (app,cmd_params))
                #if not self.reserve_intermediate:
                #    self.commands.append('rm -f %s' % fn_bam_input.rstrip('.bam')+'.ba*')
            
                fn_bam_output = fbase+'recal.bam'
                app = 'GenomeAnalysisTK.jar -T PrintReads -R %s' % self.genome_fasta
                cmd_params = '-I %s -BQSR %s -o %s' % (fn_bam_input,fn_grp_output,fn_bam_output) 
                reserved_params = ['-I','-BQSR','-o']
                user_params = self.params.get('-PrintReads',None)
                if user_params:
                    cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])            
                self.commands.append('%s %s' % (app,cmd_params))
        
                if not self.reserve_intermediate:
                    self.commands.append('rm -f %s' % fbase+'ba*')
                    self.commands.append('rm -f %s' % fn_grp_output)
                tmp.append(fn_bam_output)
            else:
                tmp.append(fn_bam_input)
        self.products = tmp   
 
    #deprecated in GATK2.0
    def TableRecalibration(self):
        """re-calibrate quality score
        
        java -Xmx4g -jar GenomeAnalysisTK.jar \
           -R resources/Homo_sapiens_assembly18.fasta \
           -knownSites bundle/hg18/dbsnp_132.hg18.vcf \
           -knownSites another/optional/setOfSitesToMask.vcf \
           -I my_reads.bam \
           -T CountCovariates \
           -cov ReadGroupCovariate \
           -cov QualityScoreCovariate \
           -cov CycleCovariate \
           -cov DinucCovariate \
           -recalFile my_reads.recal_data.csv        
        """
        tmp = []
        
        knownSites = '  '.join(['-knownSites '+i for i in (self.hg19_omni_vcf,self.hg19_hapmap_vcf,self.dbsnp)])
        #knownSites = '-knownSites %s -knownSites %s -knownSites %s' % (self.hg19_omni_vcf,self.hg19_hapmap_vcf,self.dbsnp)
        for fbam in self.products:
            fbase = fbam.rstrip('bam')
            fn_csv_output = fbase+'recal.csv'
            app = 'GenomeAnalysisTK.jar -T CountCovariates -R %s' % self.genome_fasta
            reserved_params = ['-I','-recalFile']
            cmd_params = '-I %s -recalFile %s -cov QualityScoreCovariate %s' % (fbam,fn_csv_output,knownSites)
            user_params = self.params.get('-CountCovariates',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,cmd_params))
            
            fn_bam_output = fbase+'recal.bam'
            app = 'GenomeAnalysisTK.jar -T TableRecalibration -R %s' % self.genome_fasta
            cmd_params = '-I %s -recalFile %s -o %s' % (fbam,fn_csv_output,fn_bam_output) 
            reserved_params = ['-I','-recalFile','-o']
            user_params = self.params.get('-TableRecalibration',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])            
            self.commands.append('%s %s' % (app,cmd_params))
        
            if not self.reserve_intermediate:
                self.commands.append('rm -f %s' % fbase+'ba*')
                self.commands.append('rm -f %s' % fn_csv_output)
            tmp.append(fn_bam_output)
            
        self.products = tmp   

    def DepthOfCoverage(self):
        """
        - no suffix: per locus coverage
        - _summary: total, mean, median, quartiles, and threshold proportions, aggregated over all bases
        - _statistics: coverage histograms (# locus with X coverage), aggregated over all bases
        - _interval_summary: total, mean, median, quartiles, and threshold proportions, aggregated per interval
        - _interval_statistics: 2x2 table of # of intervals covered to >= X depth in >=Y samples
        - _gene_summary: total, mean, median, quartiles, and threshold proportions, aggregated per gene
        - _gene_statistics: 2x2 table of # of genes covered to >= X depth in >= Y samples
        - _cumulative_coverage_counts: coverage histograms (# locus with >= X coverage), aggregated over all bases
        - _cumulative_coverage_proportions: proprotions of loci with >= X coverage, aggregated over all bases        
        """
        app = 'GenomeAnalysisTK.jar -T DepthOfCoverage -R %s' % self.genome_fasta
        reserved_params = ['-I','-o']
        fout = 'coverage'
        for fin in self.products:
            if fin.endswith('.bam'):
                fout = fin.rstrip('bam')+'coverage'
            else:
                raise Exception('DepthOfCoverage: Only BAM is accepted')
            
            cmd_params = '-I %s -o %s -baseCounts -omitIntervals -omitLocusTable -omitSampleSummary' % (fin,fout)
            user_params = self.params.get('-DepthOfCoverage',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,cmd_params))
            
    def UnifiedGenotyper(self,genotype):
        """java -jar GenomeAnalysisTK.jar \
               -R resources/Homo_sapiens_assembly18.fasta \
               -T UnifiedGenotyper \
               -I sample1.bam [-I sample2.bam ...] \
               --dbsnp dbSNP.vcf \
               -o snps.raw.vcf \
               -stand_call_conf [50.0] \
               -stand_emit_conf 10.0 \
               -dcov [50] \
               [-L targets.interval_list]
        """
        self.call_mode = genotype
        tmp = []
        app = 'GenomeAnalysisTK.jar -T UnifiedGenotyper -R %s -glm %s' % (self.genome_fasta,self.call_mode)
        reserved_params = ['-I','-o','--dbsnp','-glm']
        
        #print self.multiple_bams,self.products
        
        if self.multiple_bams:
            _prefix = os.path.commonprefix(self.products).strip('_').strip()
            if _prefix:
                fn_vcf_output = '%s.vcf' % _prefix
            else:
                fn_vcf_output = 'tmp.vcf'
            cmd_params = '-I %s -o %s --dbsnp %s' % (' -I '.join(self.products),fn_vcf_output,self.dbsnp)
            user_params = self.params.get('-UnifiedGenotyper',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,cmd_params))
            tmp.append(fn_vcf_output)
            
        else:
            for fin in self.products:
                fn_vcf_output = '%s' % fin.split('.')[0]+'.vcf'
                cmd_params = '-I %s -o %s --dbsnp %s' % (fin,fn_vcf_output,self.dbsnp)
                user_params = self.params.get('-UnifiedGenotyper',None)
                if user_params:
                    cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                self.commands.append('%s %s' % (app,cmd_params))
                tmp.append(fn_vcf_output)
                
        self.products = tmp
        
    def HaplotypeCaller(self,genotype):
        """
        $$$
        -T HaplotypeCaller -R reference/human_g1k_v37.fasta
             -I sample1.bam [-I sample2.bam ...] \
             --dbsnp dbSNP.vcf \
             -stand_call_conf [50.0] \
             -stand_emit_conf 10.0 \
             [-L targets.interval_list]
             -o output.raw.snps.indels.vcf
             
            -D dbSNP137.vcf 
            -A DepthOfCoverage 
            -A HaplotypeScore 
            -A MappingQualityRankSumTest 
            -A FisherStrand 
            -A ReadPosRankSumTest 
            -A QualByDepth 
            -et NO_ET 
            -K mykey 
            -L 2
        """
        tmp = []
        app = 'GenomeAnalysisTK.jar -T HaplotypeCaller -R %s' % (self.genome_fasta)
        #_annotation = ' -A '.join(('DepthOfCoverage','HaplotypeScore','MappingQualityRankSumTest','FisherStrand','ReadPosRankSumTest'))
        #_annotation = ' -A '.join(('DepthOfCoverage','MappingQualityRankSumTest','FisherStrand','ReadPosRankSumTest'))
        reserved_params = ['-I','-o','--dbsnp','-D']
        _prefix = os.path.commonprefix(self.products).strip('_').strip()
        if not _prefix:
            _prefix = 'out'
        _output_vcf = '%s.vcf' % _prefix
        _input_bam = '  '.join(['-I '+i for i in self.products])
        cmd_params = '%s -o %s --dbsnp %s' % (_input_bam,_output_vcf,self.dbsnp)
        user_params = self.params.get('-HaplotypeCaller',None)
        if user_params:
            cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
        self.commands.append('%s %s' % (app,cmd_params))
        tmp.append(_output_vcf)
        self.products = tmp
    
    def VariantAnnotator(self):
        """
        GT:
        0/0 - the sample is homozygous reference
        0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
        1/1 - the sample is homozygous alternate
        
        
        java -Xmx2g -jar GenomeAnalysisTK.jar \
           -R ref.fasta \
           -T VariantAnnotator \
           -I input.bam \
           -o output.vcf \
           -A DepthOfCoverage
           --variant input.vcf \
           --dbsnp dbsnp.vcf
        """
        
        fn_vcf_output = 'output.vcf'
        app = 'GenomeAnalysisTK.jar -T VariantAnnotator -R %s' %  self.genome_fasta
        cmd_params = '-I input.bam -o output.vcf --variant input.vcf --dbsnp dbsnp.vcf'
        return fn_vcf_output

    def VariantFiltration(self):
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
        
        tmp = []
        for fin in self.products:
            self.commands.append('vcf2snp.py %s' % fin)
            vcf_snp = fin.rstrip('vcf')+'snp.vcf' 
            vcf_indel = fin.rstrip('vcf')+'indel.vcf'
            
            vcf_snp_filter = vcf_snp.rstrip('vcf')+'filter.vcf' 
            vcf_indel_filter = vcf_indel.rstrip('vcf')+'filter.vcf' 
            
            app = 'GenomeAnalysisTK.jar -T VariantFiltration -R %s' % self.genome_fasta
            
            
            params_snp = '--variant %s -o %s -cluster 3 -window 10 --filterName "SNP" --filterExpression "QD<2.0||MQ<40.0||FS>60.0||HaplotypeScore>13.0" ' % (vcf_snp,vcf_snp_filter)
            reserved_params = ['--variant','-o']
            user_params = self.params.get('-vfs',None)
            if not user_params:
                user_params = self.params.get('-VariantFiltrationSnp',None)
                
            if user_params:
                params_snp = self.util.update_cmd_params(params_snp,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,params_snp))
            vcf_snp_pass = vcf_snp_filter.rstrip('vcf')+'pass.vcf' 
            self.commands.append("""awk '/^#/ {print $0; next} $7 == "PASS" {print $0}' %s  > %s """ % (vcf_snp_filter,vcf_snp_pass))


            params_indel = '--variant %s -o %s --filterName "INDEL" --filterExpression "QD<2.0||FS>200.0" ' % (vcf_indel,vcf_indel_filter)
            reserved_params = ['--variant','-o']
            user_params = self.params.get('-vfi',None)
            if not user_params:
                user_params = self.params.get('-VariantFiltrationIndel',None)
            if user_params:
                params_indel = self.util.update_cmd_params(params_indel,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,params_indel))
            vcf_indel_pass = vcf_indel_filter.rstrip('vcf')+'pass.vcf' 
            self.commands.append("""awk '/^#/ {print $0; next} $7 == "PASS" {print $0}' %s  > %s """ % (vcf_indel_filter,vcf_indel_pass))
            
            vcf_filter = fin.rstrip('vcf')+'filter.vcf' 
            self.commands.append('cat %s %s > %s' % (vcf_snp_pass,vcf_indel_pass,vcf_filter))
            self.commands.append('rm -f %s' % (' '.join([vcf_snp,vcf_indel,vcf_snp_filter,vcf_indel_filter,vcf_snp_pass,vcf_indel_pass])))
            tmp.append(vcf_filter)
            
        self.products = tmp
            
    def VariantRecalibrator(self,ts_filter_level='99.0'):
        """
        http://www.broadinstitute.org/gsa/gatkdocs/release/org_broadinstitute_sting_gatk_walkers_variantrecalibration_VariantRecalibrator.html
        
        java -Xmx4g -jar /bluegill/app/gatk/GenomeAnalysisTK.jar -T VariantRecalibrator -R /bluegill/data/hg19.fasta \
                       -input $x \
                       -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /bluegill/data/hg19.hapmap_3.3.vcf \
                       -resource:omni,known=false,training=true,truth=false,prior=12.0 /bluegill/data/hg19.1000G_omni2.5.vcf \
                       -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 /bluegill/data/hg19.dbsnp_132.vcf \
                       -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
                       -recalFile $x.recal \
                       -tranchesFile $x.tranches \
                       -rscriptFile $x.plots.R

        java -Xmx4g -jar /bluegill/app/gatk/GenomeAnalysisTK.jar -T ApplyRecalibration -R /bluegill/data/hg19.fasta \
        -input $x \
        --ts_filter_level 99 \
        -tranchesFile $x.tranches \
        -recalFile $x.recal \
        -o $x.2                   

        awk '$1 ~ "#" {print $0; next} $3 ~ "rs" {next} $7 ~ "TruthSensitivityTranche99.90to100.00" {print $0}' $x.2 > $x.3
        """
        tmp = []
        #self.hg19_omni_vcf,self.hg19_hapmap_vcf,self.dbsnp
        #DP,FS,MQ,MQRankSum,QD,ReadPosRankSum,

        for fin in self.products:
            fn_recal = fin.rstrip('vcf')+'recal'
            fn_tranches = fin.rstrip('vcf')+'tranches'
            #app =  'GenomeAnalysisTK.jar -T VariantRecalibrator -R %s -mode %s' % (self.genome_fasta,self.call_mode)
            app =  'GenomeAnalysisTK.jar -T VariantRecalibrator -R %s -mode BOTH' % self.genome_fasta
            
            cmd_params = '-input %s -recalFile %s -tranchesFile %s\
           -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s \
           -resource:omni,known=false,training=true,truth=false,prior=12.0 %s \
           -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 %s \
           -an QD -an DP -an ReadPosRankSum -an MQRankSum -an FS -an MQ --maxGaussians 4' % (fin,fn_recal,fn_tranches,self.hg19_hapmap_vcf,self.hg19_omni_vcf,self.dbsnp)
           
            reserved_params = ['-input','-o','-R']
            user_params = self.params.get('-VariantRecalibrator',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,cmd_params))

            #ApplyRecalibration
            fout1 = fin.rstrip('vcf')+'tmp.vcf'
            #app = 'GenomeAnalysisTK.jar -T ApplyRecalibration -R %s -mode %s' % (self.genome_fasta,self.call_mode)
            app = 'GenomeAnalysisTK.jar -T ApplyRecalibration -R %s -mode BOTH' % self.genome_fasta
            reserved_params = ['-R','-input','-o','-recalFile','-tranchesFile']
            cmd_params = '-input %s -o %s -recalFile %s -tranchesFile %s --ts_filter_level %s' % (fin,fout1,fn_recal,fn_tranches,ts_filter_level)
            user_params = self.params.get('-ApplyRecalibration',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,cmd_params))
            #fout2 = fin.rstrip('vcf')+'recal.vcf'
            #cmd_awk= """awk '$1 ~ "#" {print $0; next} $3 ~ "rs" {next} $7 ~ "TruthSensitivityTranche99.90to100.00" {print $0}' %s > %s""" % (fout1,fout2)
            
            #sensitivity_keyword = 'VQSRTrancheBOTH99.90to100.00'
            #ts_filter_level
            #cmd_awk= """awk '$1 ~ "#" {print $0; next} $3 ~ "rs" {next} $7 ~ "%s" {print $0}' %s > %s""" % (sensitivity_keyword,fout1,fout2)
            #self.commands.append(cmd_awk)
            #self.commands.append('rm -fr *.tmp.vcf')
            #tmp.append(fout2)
            fout2 = fin.rstrip('vcf')+'recal.vcf'
            app = 'GenomeAnalysisTK.jar -T SelectVariants -R %s' % self.genome_fasta
            reserved_params = ['-R','--variant','-o']
            cmd_params = '--variant %s -o %s -select "VQSLOD>4.0" ' % (fout1,fout2)
            user_params = self.params.get('-SelectVariants',None)
            if user_params:
                print '######',user_params
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,cmd_params))
            self.commands.append('rm -f %s' % fn_recal)
            self.commands.append('rm -f %s' % fn_tranches)
            self.commands.append('rm -f %s' % fout1)
            tmp.append(fout2)
            
        self.products = tmp

    def split_snp_indel(self):
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
        
        tmp = []
        for fin in self.products:
            output_snp = fin.rstrip('vcf')+'snp.vcf' 
            output_indel = fin.rstrip('vcf')+'indel.vcf' 
            self.commands.append('vcf2snp.py %s' %  fin)
            tmp.append(output_snp)
            tmp.append(output_indel)
        self.products = tmp
        
    def ClipReads(self):
        #http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_ClipReads.html
        #use "-ClipReads" for options
        #-T ClipReads -I my.bam -I your.bam -o my_and_your.clipped.bam -R Homo_sapiens_assembly18.fasta \
        #-XF seqsToClip.fasta -X CCCCC -CT "1-5,11-15" -QT 10
        
        tmp = []
        app = 'GenomeAnalysisTK.jar -T ClipReads -R %s' %  self.genome_fasta
        reserved_params = ['-I','-o']
        for fbam in self.products:
            px = re.compile(".clip.")
            if not px.findall(fbam):
                fout = fbam.rstrip('.bam')+'.clip.bam'
                cmd_params = '-I %s -o %s' % (fbam,fout)  
                user_params = self.params.get('-ClipReads',None)
                if user_params:
                    cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                self.commands.append('%s %s' % (app,cmd_params))
                if not self.reserve_intermediate:
                    self.commands.append('rm -f %s' % fbam)
                tmp.append(fout)
            else:
                tmp.append(fbam)
                
        self.products = tmp

    def CallableLoci(self):
        #-T CallableLociWalker -I my.bam -summary my.summary -o my.bed
        tmp = []
        app = 'GenomeAnalysisTK.jar -T CallableLociWalker -R %s' %  self.genome_fasta
        reserved_params = ['-I','-summary','-o']
        for fbam in self.products:
            fsummary = fbam.rstrip('.bam')+'.callableloci.summary'
            fbed = fbam.rstrip('.bam')+'.callableloci.bed'
            cmd_params = '-I %s -summary %s -o %s' % (fbam,fsummary,fbed)  
            user_params = self.params.get('-CallableLoci',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,cmd_params))
            tmp.append(fbam)
        self.products = tmp
        
    def ReduceReads(self):
        #http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_compression_reducereads_ReduceReads.html
        #-T ReduceReads -R ref.fasta -I myData.bam -o myData.reduced.bam
        tmp = []
        app = 'GenomeAnalysisTK.jar -T ReduceReads -R %s' %  self.genome_fasta
        reserved_params = ['-I','-o']
        for fn_bam_input in self.products:
            px = re.compile(".reduce.")
            if not px.findall(fn_bam_input):
                fout = fn_bam_input.rstrip('.bam')+'.reduce.bam'
                cmd_params = '-I %s -o %s' % (fn_bam_input,fout)  
                user_params = self.params.get('-ReduceReads',None)
                if user_params:
                    cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                self.commands.append('%s %s' % (app,cmd_params))
                tmp.append(fout)
            else:
                tmp.append(fn_bam_input)
        self.products = tmp

    def CombineVariants(self):
        #http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html            
        #-T CombineVariants -R ref.fasta --variant input1.vcf --variant input2.vcf -o output.vcf -genotypeMergeOptions UNIQUIFY
        app = 'GenomeAnalysisTK.jar -T CombineVariants -R %s' %  self.genome_fasta
        reserved_params = ['--variant','-o']
        if len(self.products)==1:
            return
        else:
            #cmd_params = '-I %s -o %s' % (fbam,fout)  
            cmd_params = ''  
            common_prefix = os.path.commonprefix(self.products).strip('_').strip()
            if not common_prefix:
                common_prefix = 'merge'
            fout = '%s.vcf' % common_prefix
            finput = '--variant '.join(self.products)
            user_params = self.params.get('-CombineVariants',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s' % (app,cmd_params))
            if not self.reserve_intermediate:
                self.commands.append('rm -f %s' % ' '.join(self.products))
                
            self.products = [fout]
                    
    #PICARD
    def SortSam(self):
        tmp = []
        for fsam in self.products:
            px = re.compile(".sort.bam")
            if not px.findall(fsam):
                fout = None
                if fsam.endswith('.sam'):
                    fbase = fsam.rstrip('.sam')
                elif fsam.endswith('.sam.gz'):
                    fout = fsam.rstrip('.sam.gz')
                elif fsam.endswith('.bam'):
                    fbase = fsam.rstrip('.bam')
                else:
                    raise Exception('SortSam: Only SAM or BAM is accepted')
                fout = fbase+'.sort.bam'
                self.commands.append('SortSam.jar INPUT=%s OUTPUT=%s CREATE_INDEX=true SO=coordinate COMPRESSION_LEVEL=5 MAX_RECORDS_IN_RAM=5000000 VERBOSITY=ERROR QUIET=true TMP_DIR=%s VALIDATION_STRINGENCY=%s' % (fsam,fout,self.util.CLUSTER_TMP_DIR,self.validation_stringency))
                if not self.reserve_intermediate:
                    self.commands.append('rm -f %s' % fout.rstrip('.sort.bam')+'.ba*')
                tmp.append(fout)
            else:
                tmp.append(fsam)
                
        self.products = tmp
    
    def BuildBamIndex(self):
        for fbam in self.products:
            self.commands.append("BuildBamIndex.jar INPUT=%s VERBOSITY=ERROR QUIET=true VALIDATION_STRINGENCY=%s" % (fbam,self.validation_stringency)) 
    
    def MergeLanes(self,same_sample=False):
        """
        merge multiple lane BAMs into one sample BAM.
        """
        tmp = []
        if len(self.products)==1:
            return
        else:
            m = {}
            for i in self.products:
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

                self.commands.append('samtools merge %s %s' % (fout,' '.join(v)))
                if not self.reserve_intermediate:
                    self.commands.append('rm -f %s' % ' '.join([i.rstrip('.bam')+'.ba*' for i in v]))
                    
                tmp.append(fout)
            self.products = tmp 
            self.SortSam()
       
    def FixMateInformation(self):
        tmp = []
        for fbam in self.products:
            px = re.compile(".mate.")
            if not px.findall(fbam):
                fout = None
                if fbam.endswith('.sam'):
                    fout = fbam.rstrip('.sam')+'.mate.bam'
                elif fbam.endswith('.sam.gz'):
                    fout = fbam.rstrip('.sam.gz')+'.mate.bam'
                elif fbam.endswith('.bam'):
                    fout = fbam.rstrip('.bam')+'.mate.bam'
                else:
                    raise Exception('FixMateInformation: Only SAM or BAM is accepted')
                self.commands.append('FixMateInformation.jar INPUT=%s OUTPUT=%s TMP_DIR=%s SO=coordinate VERBOSITY=ERROR QUIET=true CREATE_INDEX=true VALIDATION_STRINGENCY=%s' % (fbam,fout,self.util.CLUSTER_TMP_DIR,self.validation_stringency))
                if not self.reserve_intermediate:
                    self.commands.append('rm -f %s' % fbam)
                tmp.append(fout)
            else:
                tmp.append(fbam)
        self.products = tmp
    
    def MarkDuplicates(self):
        tmp = []
        for fbam in self.products:
            px = re.compile(".dup.")
            if not px.findall(fbam):
                fout = None
                if fbam.endswith('.sam'):
                    fout = fbam.rstrip('.sam')+'.dup.bam'
                elif fbam.endswith('sam.gz'):
                    fout = fbam.rstrip('.sam.gz')+'.dup.bam'
                elif fbam.endswith('.bam'):
                    fout = fbam.rstrip('.bam')+'.dup.bam'
                else:
                    raise Exception('MarkDuplicates: Only SAM or BAM is accepted')
                self.commands.append('MarkDuplicates.jar INPUT=%s OUTPUT=%s M=%s.duplicate VERBOSITY=ERROR QUIET=true REMOVE_DUPLICATES=true ASSUME_SORTED=true TMP_DIR=%s VALIDATION_STRINGENCY=%s' % (fbam,fout,fbam,self.util.CLUSTER_TMP_DIR,self.validation_stringency))
                if not self.reserve_intermediate:
                    self.commands.append('rm -f %s' % fbam.rstrip('.bam')+'.ba*')
                tmp.append(fout)
            else:
                tmp.append(fbam)
                
        self.products = tmp

    #end of picard tools
    #######################################################################

class Annotator(App):
    def __init__(self,reference_genome,background_variance_files,target_variance_files,params=None,reserve_intermediate_files=False):
        App.__init__(self,[background_variance_files,target_variance_files],params,reserve_intermediate_files)
        self.genome = reference_genome
        self.bg_vcfs = background_variance_files
        self.tg_vcfs = target_variance_files

class snpEff(Annotator):
    #http://snpeff.sourceforge.net/manual.html
    #http://downloads.sourceforge.net/project/snpeff/snpEff_latest_core.zip?r=http%3A%2F%2Fsnpeff.sourceforge.net%2Fdownload.html&ts=1344880466&use_mirror=voxel
    #http://sourceforge.net/projects/snpeff/files/databases/v3.0/snpEff_v3.0_hg19.zip
    def __init__(self,reference_genome,background_variance_files,target_variance_files,params=None,reserve_intermediate_files=False):
        Annotator.__init__(self,reference_genome,background_variance_files,target_variance_files,params,reserve_intermediate_files)
        self.data_dir = 'snpEff'
        self.config_file = 'snpEff.config'
        #self.gene_model = 'knownGene'
        
    def annotate(self):
        #copy the config file for snpEff
        cmd = 'cp %s .' % os.path.join(self.util.CLUSTER_DATA_DIR,self.config_file)
        self.commands.append(cmd)
        
        #copy the data with dir to current directory
        cmd = 'cp -R %s .' % os.path.join(self.util.CLUSTER_DATA_DIR,self.data_dir)
        
        self.commands.append(cmd)
        
        reserved_params = []
        merge_type = self.params.get('-m',None)

        if not self.bg_vcfs: #no background, just annotate any input files use ANNOVAR
            if merge_type:
                common_prefix = os.path.commonprefix(self.tg_vcfs).strip()
                if not common_prefix:
                    common_prefix = 'merge'
                    
                if common_prefix.endswith(os.sep): #only the folder is same
                    common_prefix += 'merge'
                
                fvcf_merge = '%s.vcf' % common_prefix
                self.commands.append('vcfop.py -i %s > %s' % (' '.join(self.tg_vcfs),fvcf_merge))
                app = 'snpEff.jar %s %s' % (self.genome,fvcf_merge)
                cmd_params = '-o txt' 
                user_params = self.params.get('-snpeff',None)
                if user_params:
                    cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                    
                fo_txt = '%s.txt' % common_prefix
                self.commands.append('%s %s > %s' % (app,cmd_params,fo_txt)) 
                
                
            else: #do not merge, just annotate each VCF
                for fn_vcf in self.tg_vcfs:
                    app = 'snpEff.jar %s %s' % (self.genome,fn_vcf)
                    cmd_params = '-o txt'
                    user_params = self.params.get('-snpeff',None)
                    if user_params:
                        cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                    fo_txt = fn_vcf.rstrip('vcf')+'txt'
                    self.commands.append('%s %s > %s' % (app,cmd_params,fo_txt)) 

        else: #case-control
            #fn_vcf_intersection = 'case-control.vcf' 

            fn_vcf_control = 'c.snpeff.vcf'
             
            fn_vcf_treatment = 't.snpeff.vcf'
             
            fn_vcf_treatment_wo_control = 't-c.snpeff.vcf'
             
            self.commands.append('vcfop.py -i %s > %s' % (' '.join(self.bg_vcfs),fn_vcf_control))
            
            self.commands.append('vcfop.py -i %s > %s' % (' '.join(self.tg_vcfs),fn_vcf_treatment))
            
            self.commands.append('vcfop.py -c %s %s > %s' % (fn_vcf_treatment,fn_vcf_control,fn_vcf_treatment_wo_control))
            
            
            #merge_type = self.params.get('-m',None)
            #if merge_type:
            #    self.commands.append('merge_vcf.py -control %s -case %s -m %s > %s' % (' '.join(self.bg_vcfs),' '.join(self.tg_vcfs),merge_type[0],fn_vcf_intersection))
            #else:
            #    self.commands.append('merge_vcf.py -control %s -case %s -m ii > %s' % (' '.join(self.bg_vcfs),' '.join(self.tg_vcfs),fn_vcf_intersection))
        
            app = 'snpEff.jar %s %s' % (self.genome,fn_vcf_treatment_wo_control)
            cmd_params = '-o txt'
            user_params = self.params.get('-snpeff',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            fo_txt = fn_vcf_treatment_wo_control.rstrip('vcf')+'txt'
            self.commands.append('%s %s > %s' % (app,cmd_params,fo_txt)) 
        
        #remove the snpEff data dir
        cmd = 'rm -fr %s' % self.data_dir
        self.commands.append(cmd)
        
        #remove the config file
        cmd = 'rm -f %s' % self.config_file
        self.commands.append(cmd)
        
        
    def get_command(self):
        self.annotate()
        return self.commands
            
class ANNOVAR(Annotator):
    """
    tabix: wget http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.5.tar.bz2/download
    vcftools: wget http://sourceforge.net/projects/vcftools/files/vcftools_0.1.7.tar.gz/download
    annovar: wget http://www.openbioinformatics.org/annovar/download/annovar.latest.tar.gz.mirror
    
    tar -zxvf annovar.tar.gz
    cd annovar
    annotate_variation.pl -downdb -buildver hg19 refgene humandb
    annotate_variation.pl -downdb -buildver hg19 knownGene humandb               
    annotate_variation.pl -downdb -buildver hg19 ensGene humandb
    convert2annovar.pl x2.vcf -format vcf4 > x2.avr

    #NCBI RefSeq
    annotate_variation.pl --buildver hg19 -geneanno -dbtype refgene x2.avr human_db/
    
    #UCSC KnownGene
    annotate_variation.pl --buildver hg19 -geneanno -dbtype knownGene x2.avr human_db/
    
    #UCSC Gene
    annotate_variation.pl --buildver hg19 -geneanno -dbtype ensGene x2.avr human_db/
    
    
    """
    def __init__(self,reference_genome,background_variance_files,target_variance_files,params=None,reserve_intermediate_files=False):
        Annotator.__init__(self,reference_genome,background_variance_files,target_variance_files,params,reserve_intermediate_files)
        self.gene_model = 'knownGene'
        #self.gene_model = 'refgene'
        #self.gene_model = 'ensGene'

    def annotate(self):
        #from misc import CLUSTER_DATA_DIR
        fn_gene = ''
        fn_seq = ''
        fn_ref = ''
        path_libfiles='.'
        
        if self.gene_model=='refgene':
            fn_gene = 'hg19_refGene.txt'
            fn_seq = 'hg19_refGeneMrna.fa'
            fn_ref = 'refLink.txt.gz'
            
        elif self.gene_model=='knownGene':
            fn_gene = 'hg19_knownGene.txt'
            fn_seq = 'hg19_knownGeneMrna.fa'
            fn_ref = 'hg19_kgXref.txt'
            
        elif self.gene_model=='ensGene':
            fn_gene = 'hg19_ensGene.txt'
            fn_seq = 'hg19_ensGeneMrna.fa'
            
        else:
            raise Exception('Unsupported gene model: %s' % self.gene_model)
            
        cmd = 'cp %s .' % os.path.join(self.util.CLUSTER_DATA_DIR,fn_gene)
        self.commands.append(cmd)
        
        cmd = 'cp %s .' % os.path.join(self.util.CLUSTER_DATA_DIR,fn_seq)
        self.commands.append(cmd)
        
        if fn_ref:
            cmd = 'cp %s .' % os.path.join(self.util.CLUSTER_DATA_DIR,fn_ref)
            self.commands.append(cmd)
            
        #self.params
        reserved_params = []
        merge_type = self.params.get('-m',None)

        if not self.bg_vcfs: #no background, just annotate any input files use ANNOVAR
            if merge_type:
                common_prefix = os.path.commonprefix(self.tg_vcfs).strip()
                
                if not common_prefix:
                    common_prefix = 'merge'
                    
                if common_prefix.endswith(os.sep): #only the folder is same
                    common_prefix += 'merge'
                
                fvcf_merge = '%s.vcf' % common_prefix
                self.commands.append('vcfop.py -i %s > %s' % (' '.join(self.tg_vcfs),fvcf_merge))
                fn_annovar = fvcf_merge.rstrip('vcf')+'annovar'
                
                cmd = 'convert2annovar.pl %s -format vcf4 -includeinfo > %s' % (fvcf_merge,fn_annovar)
                self.commands.append(cmd)
                app = 'annotate_variation.pl --buildver %s -geneanno' % self.genome
                cmd_params = '-dbtype %s' % self.gene_model
                user_params = self.params.get('-annovar',None)
                if user_params:
                    cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                self.commands.append('%s %s %s %s' % (app,cmd_params,fn_annovar,path_libfiles)) 
                
                
            else: #do not merge, just annotate each VCF
                for fn_vcf in self.tg_vcfs:
                    fn_annovar = fn_vcf.rstrip('vcf')+'annovar'
                    
                    cmd = 'convert2annovar.pl %s -format vcf4 -includeinfo > %s' % (fn_vcf,fn_annovar)
                    self.commands.append(cmd)
                    app = 'annotate_variation.pl --buildver %s -geneanno' % self.genome
                    cmd_params = '-dbtype %s' % self.gene_model
                    user_params = self.params.get('-annovar',None)
                    if user_params:
                        cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
                    self.commands.append('%s %s %s %s' % (app,cmd_params,fn_annovar,path_libfiles)) 

        else: #case-control
            #fn_vcf_intersection = 'my.vcf' 
            fn_vcf_control = 'c.annovar.vcf' 
            fn_vcf_treatment = 't.annovar.vcf' 
            fn_vcf_treatment_wo_control = 't-c.annovar.vcf' 
            self.commands.append('vcfop.py -i %s > %s' % (' '.join(self.bg_vcfs),fn_vcf_control))
            self.commands.append('vcfop.py -i %s > %s' % (' '.join(self.tg_vcfs),fn_vcf_treatment))
            self.commands.append('vcfop.py -c %s %s > %s' % (fn_vcf_treatment,fn_vcf_control,fn_vcf_treatment_wo_control))
            
            #merge_type = self.params.get('-m',None)
            #if merge_type:
                #self.commands.append('merge_vcf.py -control %s -case %s -m %s > %s' % (' '.join(self.bg_vcfs),' '.join(self.tg_vcfs),merge_type[0],fn_vcf_intersection))
            #else:
                #self.commands.append('merge_vcf.py -control %s -case %s -m ii > %s' % (' '.join(self.bg_vcfs),' '.join(self.tg_vcfs),fn_vcf_intersection))
        
            fn_annovar = fn_vcf_treatment_wo_control.rstrip('vcf')+'annovar'
            cmd = 'convert2annovar.pl %s -format vcf4 > %s' % (fn_vcf_treatment_wo_control,fn_annovar)
            self.commands.append(cmd)
            app = 'annotate_variation.pl --buildver %s -geneanno' % self.genome
            cmd_params = '-dbtype %s' % self.gene_model
            user_params = self.params.get('-annovar',None)
            if user_params:
                cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
            self.commands.append('%s %s %s %s' % (app,cmd_params,fn_annovar,path_libfiles)) 

        cmd = 'rm -f %s' % fn_gene
        self.commands.append(cmd)
        
        cmd = 'rm -f %s' % fn_seq
        self.commands.append(cmd)
        
        if fn_ref:
            cmd = 'rm -f %s' % fn_ref
            self.commands.append(cmd)
    

    def get_command(self):
        self.annotate()
        return self.commands
    
class VAAST(Annotator):
    
    def __init__(self,reference_genome,background_variance_files,target_variance_files,params=None,reserve_intermediate_files=False):
        Annotator.__init__(self,reference_genome,background_variance_files,target_variance_files,params,reserve_intermediate_files)
        self.genome_fasta = '%s.fasta' % self.genome
        self.gene_reference = '%s.refgene.gff3' % self.genome
        
    def Vcf2Gvf(self,file_vcf):
        output = file_vcf.rstrip('vcf')+'gvf'
        cmd = 'vaast_converter %s -s -b %s > %s' % (file_vcf,self.genome,output)
        self.commands.append(cmd)
        return output
        
    def Vat(self,fn_gvf_input):
        fn_gvf_output = fn_gvf_input.rstrip('gvf')+'vat.gvf' 
        cmd = 'VAT -f %s -a %s  %s > %s' % (self.gene_reference,self.genome_fasta,fn_gvf_input,fn_gvf_output)
        self.commands.append(cmd)
        return fn_gvf_output
    
    def Vst(self,fn_gvf_inputs):                
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
        import time        
        num_gvfs = len(fn_gvf_inputs)        
        output = os.path.commonprefix(fn_gvf_inputs).rstrip('_').rstrip('.vat.gvf')
        if not output:
            output = 'tmp.'+str(time.time())
            time.sleep(1)            
        fn_cdr_output = output+'.cdr'
        if num_gvfs>1:
            cmd = "VST  -o 'I(0..%d)' -b %s %s > %s" % (num_gvfs-1,self.genome,' '.join(fn_gvf_inputs),fn_cdr_output)
        else: 
            cmd = "VST  -o 'I(0)' -b %s %s > %s" % (self.genome,' '.join(fn_gvf_inputs),fn_cdr_output)            
        self.commands.append(cmd)        
        return fn_cdr_output
                
    def Vaast(self,bg_cdr,tg_cdr):
        output = os.path.commonprefix((bg_cdr,tg_cdr)).rstrip('_').rstrip('.cdr')        
        output = 't-c'
        if not output:
            output = 'my'
        fn_vaast_output = output+'.vaast'
        app = 'VAAST'
        #cmd_params = "-iht r -lh n -fast_gp -d 1e4 -r 0.00035 -m lrt -k %s -o %s %s %s 2>%s.log" % (self.genome,fn_vaast_output,bg_cdr,tg_cdr,output)
        cmd_params = "-m lrt -o %s %s %s %s 2>vaast.log" % (output,self.gene_reference,bg_cdr,tg_cdr)
        reserved_params = {'-o'}
        user_params = self.params.get('-vaast',None)
        if user_params:
            cmd_params = self.util.update_cmd_params(cmd_params,reserved_params,user_params[0])
        self.commands.append('%s %s' % (app,cmd_params))
        return fn_vaast_output
        
    def get_command(self):
        bg_gvfs = []
        for b in self.bg_vcfs:
            m = self.Vat(self.Vcf2Gvf(b))
            bg_gvfs.append(m)
        
        x = self.Vst(bg_gvfs)
        tg_gvfs = []
        for t in self.tg_vcfs:
            m = self.Vat(self.Vcf2Gvf(t))
            tg_gvfs.append(m)
        y = self.Vst(tg_gvfs)
        self.Vaast(x,y)
        return self.commands
