# check for logic in osa! make sure there is none!

# make cool plots

# reduce logs to essential; marker stream on topic

# output effect

# evaluate time, estimate time to do

# gzip files always

# simlink import

# avoid optional caches

# deep cacges

# asynch logs
#
# there are two ways to define common origin:
#   define during class declaratopm
# dynamic redefinition of the classes, like ScWData, will invalidate (hopefully) the get call, so it should be reexecuted

# there is special treatement

# import as

# maks environment, osa, as input
# make aliases for evolutions

from dataanalysis import DataFile,shhash
import dataanalysis 


from pilton import heatool

import pprint
import os,shutil,re,time,glob
from astropy.io import fits as pyfits

def remove_repeating(inlist):
    if inlist==[]:
        return inlist

    outlist=[inlist[0]]
    for l in inlist[1:]:
        if l in outlist:
            print("repearing item:",l)
        else:
            outlist.append(l)
    return outlist

class MemCacheIntegralBaseOldPath: 
#class MemCacheIntegral(dataanalysis.MemCacheSqlite):
    def get_scw(self,hashe):                                                                                                                                       
        #if dataanalysis.printhook.global_log_enabled: print("search for scw in",hashe)
        if isinstance(hashe,tuple):                                                                                                                                
            if hashe[0]=="analysis": # more universaly                                                                                                             
                if hashe[2].startswith('ScWData'):
                    return hashe[1]
                return self.get_scw(hashe[1])
            if hashe[0]=="list": # more universaly                                                                                                                 
                for k in hashe[1:]:
                    r=self.get_scw(k)
                    if r is not None:
                        return r
                return None
            raise Exception("unknown tuple in the hash:"+str(hashe))                                                                                               
        if hashe is None:
            return None # 'Any'
        if isinstance(hashe,str):                                                                                                                                  
            return None                                                                                                                                           
        raise Exception("unknown class in the hash:"+str(hashe))                                                                                                   
    
    def get_rev(self,hashe):                                                                                                                                       
        #if dataanalysis.printhook.global_log_enabled: print("search for rev in",hashe)
        if isinstance(hashe,tuple):                                                                                                                                
            if hashe[0]=="analysis": # more universaly                                                                                                             
                if hashe[2].startswith('Revolution'):
                    return hashe[1]
                return self.get_rev(hashe[1])
            if hashe[0]=="list": # more universaly                                                                                                                 
                for k in hashe[1:]:
                    r=self.get_rev(k)
                    if r is not None:
                        return r
                return None
            raise Exception("unknown tuple in the hash:"+str(hashe))                                                                                               
        if hashe is None:
            return None
        if isinstance(hashe,str):                                                                                                                                  
            return None                                                                                                                                           
        raise Exception("unknown class in the hash:"+str(hashe))                                                                                                   
    
    def get_marked(self,hashe):                                                                                                                                       
        #if dataanalysis.printhook.global_log_enabled: print("search for marked in",hashe)
        if isinstance(hashe,tuple):                                                                                                                                
            if hashe[0]=="analysis": # more universaly                                                                                                             
                r=[]
                if hashe[2].endswith('..'):
                    r.append(hashe[2][:-2])
                r+=self.get_marked(hashe[1])
                return r
            if hashe[0]=="list": # more universaly                                                                                                                 
                r=[]
                for k in hashe[1:]:
                    r+=self.get_marked(k)
                return r
            raise Exception("unknown tuple in the hash:"+str(hashe))                                                                                               
        if hashe is None:
            return []
        if isinstance(hashe,str):                                                                                                                                  
            return []                                                                                                                                   
        raise Exception("unknown class in the hash:"+str(hashe))                                                                                                   

    def hashe2signature(self,hashe):                                                                                                           
        scw=self.get_scw(hashe)
        if scw is not None:
            if isinstance(hashe,tuple):          
                if hashe[0]=="analysis":                                                                                                           
                    return hashe[2]+":"+scw+":"+shhash(hashe)[:8]                                                                                          
            return shhash(hashe)[:8]                                                                                                               
        
        rev=self.get_rev(hashe)
        if rev is not None:
            if isinstance(hashe,tuple):          
                if hashe[0]=="analysis":                                                                                                           
                    return hashe[2]+":"+rev+":"+shhash(hashe)[:8]                                                                                          
            return shhash(hashe)[:8]                                                                                                               

        return hashe[2]+":"+shhash(hashe)[:16]                                               
        


    def construct_cached_file_path(self,hashe,obj=None):                                                                                                                        
        #print "will construct INTEGRAL cached file path for",hashe

        scw=self.get_scw(hashe)
        rev=self.get_rev(hashe)

        def hash_to_path2(hashe):                                                                                                                                      
            return dataanalysis.shhash(repr(hashe[1]))[:8]                                                                                                           
            
        marked=self.get_marked(hashe[1])
        marked=remove_repeating(marked)
        if dataanalysis.printhook.global_log_enabled: print("marked",marked)
        
        if not isinstance(scw,str):
            print("emergent scw:",scw)
            scw=None
        
        if scw=="Any":
            print("any scw:",scw,hashe)
            scw=None

        print "scw:",scw
        print "rev:",rev


        if scw is None:
            if dataanalysis.printhook.global_log_enabled: print("not scw-grouped cache")
            if rev is None:
                if dataanalysis.printhook.global_log_enabled: print("not rev-grouped cache")
                r=self.filecacheroot+"/global/"+hashe[2]+"/"+"/".join(marked)+"/"+hash_to_path2(hashe)+"/"
            else:
                if dataanalysis.printhook.global_log_enabled: print("cached rev:",rev)
                r=self.filecacheroot+"/byrev/"+rev+"/"+hashe[2]+"/"+"/".join(marked)+"/"+hash_to_path2(hashe)+"/" # choose to avoid overlapp    
        else:
            if dataanalysis.printhook.global_log_enabled: print("cached scw:",scw)
            print(scw,hashe[2],marked)
            r=self.filecacheroot+"/byscw/"+scw[:4]+"/"+scw+"/"+hashe[2]+"/"+"/".join(marked)+"/"+hash_to_path2(hashe)+"/" # choose to avoid overlapp    

        #if dataanalysis.printhook.global_log_enabled: print("cached path:",r)

        print self,"cached path:",r
                                                                                                                                                                       
        return r # choose to avoid overlapp    

class MemCacheIntegralBase: 
    def get_scw(self,hashe):                                                                                                                                       
        if isinstance(hashe,tuple):                                                                                                                                
            if hashe[0]=="analysis": # more universaly                                                                                                             
                if hashe[2].startswith('ScWData'):
                    return hashe[1]
                return self.get_scw(hashe[1])
            if hashe[0]=="list": # more universaly                                                                                                                 
                for k in hashe[1:]:
                    r=self.get_scw(k)
                    if r is not None:
                        return r
                return None
            raise Exception("unknown tuple in the hash:"+str(hashe))                                                                                               
        if hashe is None:
            return None #'Any'
        if isinstance(hashe,str):                                                                                                                                  
            return None                                                                                                                                           
        raise Exception("unknown class in the hash:"+str(hashe))                                                                                                   
    
    def get_rev(self,hashe):                                                                                                                                       
        #if dataanalysis.printhook.global_log_enabled: print("search for rev in",hashe)
        if isinstance(hashe,tuple):                                                                                                                                
            if hashe[0]=="analysis": # more universaly                                                                                                             
                if hashe[2].startswith('Revolution'):
                    return hashe[1]
                return self.get_rev(hashe[1])
            if hashe[0]=="list": # more universaly                                                                                                                 
                for k in hashe[1:]:
                    r=self.get_rev(k)
                    if r is not None:
                        return r
                return None
            raise Exception("unknown tuple in the hash:"+str(hashe))                                                                                               
        if hashe is None:
            return None
        if isinstance(hashe,str):                                                                                                                                  
            return None                                                                                                                                           
        raise Exception("unknown class in the hash:"+str(hashe))                                                                                                   
    
    def get_marked(self,hashe):                                                                                                                                       
        #if dataanalysis.printhook.global_log_enabled: print("search for marked in",hashe)
        if isinstance(hashe,tuple):                                                                                                                                
            if hashe[0]=="analysis": # more universaly                                                                                                             
                r=[]
                if hashe[2].endswith('..'):
                    r.append(hashe[2][:-2])
                r+=self.get_marked(hashe[1])
                return r
            if hashe[0]=="list": # more universaly                                                                                                                 
                r=[]
                for k in hashe[1:]:
                    r+=self.get_marked(k)
                return r
            raise Exception("unknown tuple in the hash:"+str(hashe))                                                                                               
        if hashe is None:
            return []
        if isinstance(hashe,str):                                                                                                                                  
            return []                                                                                                                                   
        raise Exception("unknown class in the hash:"+str(hashe))                                                                                                   

        
    def hashe2signature(self,hashe):                                                                                                           
        scw=self.get_scw(hashe)
        if scw is not None:
            if isinstance(hashe,tuple):          
                if hashe[0]=="analysis":                                                                                                           
                    return hashe[2]+":"+scw+":"+shhash(hashe)[:8]                                                                                          
            return shhash(hashe)[:8]                                                                                                               
        
        rev=self.get_rev(hashe)
        if rev is not None:
            if isinstance(hashe,tuple):          
                if hashe[0]=="analysis":                                                                                                           
                    return hashe[2]+":"+rev+":"+shhash(hashe)[:8]                                                                                          
            return shhash(hashe)[:8]                                                                                                               

        return hashe[2]+":"+shhash(hashe)[:16]                                               


    def construct_cached_file_path(self,hashe,obj=None):                                                                                                                        
        #print "will construct INTEGRAL cached file path for",hashe

        scw=self.get_scw(hashe)
        rev=self.get_rev(hashe)

        def hash_to_path2(hashe):                                                                                                                                      
            return dataanalysis.shhash(repr(hashe[1]))[:8]                                                                                                           
            
        marked=self.get_marked(hashe[1])
        marked=remove_repeating(marked)
        if dataanalysis.printhook.global_log_enabled: print("marked",marked)

        for mark in marked:
            hashe=dataanalysis.hashe_replace_object(hashe,mark+"..","any")

        if not isinstance(scw,str):
            scw=None
        
        if scw=="Any":
            scw=None

        if scw is None:
            if dataanalysis.printhook.global_log_enabled: print("not scw-grouped cache")
            if rev is None:
                if dataanalysis.printhook.global_log_enabled: print("not rev-grouped cache")
                r=self.filecacheroot+"/global/"+hashe[2]+"/"+"/".join(marked)+"/"+hash_to_path2(hashe)+"/"
            else:
                hashe=dataanalysis.hashe_replace_object(hashe,rev,"any")
                #print "reduced hashe",hashe
                if dataanalysis.printhook.global_log_enabled: print("cached rev:",rev)
                r=self.filecacheroot+"/byrev/"+rev+"/"+hashe[2]+"/"+"/".join(marked)+"/"+hash_to_path2(hashe)+"/" # choose to avoid overlapp    
        else:
            hashe=dataanalysis.hashe_replace_object(hashe,scw,"any")
            hashe=dataanalysis.hashe_replace_object(hashe,('analysis', scw[:4], 'Revolution.v0'),('analysis', 'any', 'Revolution.v0'))

            #str("reduced hashe:",hashe,hash_to_path2(hashe))
            #if dataanalysis.printhook.global_log_enabled: print("reduced hashe:",hashe,hash_to_path2(hashe))
            open("reduced_hashe.txt","w").write(hash_to_path2(hashe)+"\n\n"+pprint.pformat(hashe)+"\n")
            print(scw,hashe[2],marked)

            r=self.filecacheroot+"/byscw/"+scw[:4]+"/"+scw+"/"+hashe[2]+"/"+"/".join(marked)+"/"+hash_to_path2(hashe)+"/" # choose to avoid overlapp    

        #if dataanalysis.printhook.global_log_enabled: print("cached path:",r)

        print self,"cached path:",r
                                                                                                                                                                       
        return r # choose to avoid overlapp    


class MemCacheIntegral(MemCacheIntegralBase,dataanalysis.MemCacheMySQL):
    pass

class MemCacheIntegralLegacy(MemCacheIntegralBase,dataanalysis.MemCacheSqlite):
    pass

class MemCacheIntegralFallback(MemCacheIntegralBase,dataanalysis.MemCacheNoIndex):
    pass

class MemCacheIntegralFallbackOldPath(MemCacheIntegralBaseOldPath,dataanalysis.MemCacheNoIndex):
    readonly_cache=True

class MemCacheIntegralIRODS(MemCacheIntegralBase,dataanalysis.MemCacheIRODS):
    pass

#mc=dataanalysis.TransientCacheInstance
#mcg=MemCacheIntegral('/Integral/data/reduced/ddcache/')
#mc=mcg
#mc.parent=mcg
#mcgl=MemCacheIntegralLegacy('/Integral/data/reduced/ddcache/')
#mcg.parent=mcgl

try:
    IntegralCacheRoots=os.environ['INTEGRAL_DDCACHE_ROOT']
except:
    IntegralCacheRoots='/sps/integral/data/reduced/ddcache/'

CacheStack=[]

for IntegralCacheRoot in IntegralCacheRoots.split(":"):
    ro_flag=False
    if IntegralCacheRoot.startswith("ro="):
        ro_flag=True
        IntegralCacheRoot=IntegralCacheRoot.replace("ro=","")

    mcgfb=MemCacheIntegralFallback(IntegralCacheRoot)
    mcgfb.readonly_cache=ro_flag
    if CacheStack==[]:
        CacheStack=[mcgfb]
    else:
        CacheStack[-1].parent=mcgfb
        CacheStack.append(mcgfb)

    mcgfb_oldp=MemCacheIntegralFallbackOldPath(IntegralCacheRoot)
    mcgfb_oldp.readonly_cache=ro_flag
    mcgfb.parent=mcgfb_oldp
    CacheStack.append(mcgfb_oldp)

mcgirods=MemCacheIntegralIRODS('/tempZone/home/integral/data/reduced/ddcache/')
CacheStack[-1].parent=mcgirods
CacheStack.append(mcgirods)

mc=CacheStack[0]

print "cache stack:",CacheStack

class DataAnalysis(dataanalysis.DataAnalysis):
    cache=mc

    write_caches=[dataanalysis.TransientCache,MemCacheIntegralFallback]
    read_caches=[dataanalysis.TransientCache,MemCacheIntegralFallback,MemCacheIntegralFallbackOldPath]

    cached=False

    def get_scw(self):
        if self._da_locally_complete is not None:
            try:
                return "(completescw:%s)"%self.cache.get_scw(self._da_locally_complete)
            except:
                return "(complete)"

        for a in self.assumptions:
            if isinstance(a,ScWData):
                return "(assumescw:%s)"%str(a.input_scwid)

        return ""

    def __repr__(self):
        return "[%s%s%s%s%i]"%(self.get_version(),self.get_scw(),";Virtual" if self.virtual else "",";Complete" if self._da_locally_complete else "",id(self))

class NoScWData(dataanalysis.AnalysisException):
    pass
    

class ScWData(DataAnalysis):
    input_scwid=None

    cached=False # how do we implment that this can change?

    schema_hidden=True

    version="v1"


    def main(self):
        try:
            self.scwid=self.input_scwid.handle
        except:
            self.scwid=self.input_scwid
        self.scwver=self.scwid[-3:]
        self.revid=self.scwid[:4]
        
        try:
            print "searching in "+os.environ['REP_BASE_PROD']
            self.assume_rbp(os.environ['REP_BASE_PROD'])
        except dataanalysis.AnalysisException:
            if self.scwver=="000":
                print "searching in "+os.environ['REP_BASE_PROD']+"/nrt"
                self.assume_rbp(os.environ['REP_BASE_PROD']+"/nrt")
            else:
                raise

    def assume_rbp(self,rbp):
        self.scwpath=rbp+"/scw/"+self.revid+"/"+self.scwid #!!!!
        self.revdirpath=rbp+"/scw/"+self.revid+"/rev."+self.scwver # ver?
        self.auxadppath=rbp+"/aux/adp/"+self.revid+"."+self.scwver

        if not os.path.exists(self.scwpath+"/swg.fits"):
            if not os.path.exists(self.scwpath+"/swg.fits.gz"):
                print "failed searching for",self.scwpath+"/swg.fits"
                raise NoScWData("no scw data for: "+repr(self.scwid))
                #raise Exception("no scw data!")
            else:
                self.swgpath=self.scwpath+"/swg.fits.gz"
        else:
            self.swgpath=self.scwpath+"/swg.fits"
        print "swgpath:",self.swgpath

    def get_isgri_events(self):
        if hasattr(self,'isgrievents'):
            return self.isgrievents.get_path()
        return self.scwpath+"/isgri_events.fits.gz"

    def get_telapse(self):
        return pyfits.open(self.swgpath)[1].header['TELAPSE']
    
    def get_t(self):
        h=pyfits.open(self.swgpath)[1].header
        return (h['TSTOP']+h['TSTART'])/2.,(h['TSTOP']-h['TSTART'])/2.
    
    def get_t1_t2(self):
        h=pyfits.open(self.swgpath)[1].header
        return h['TSTART'],h['TSTOP']

    def __repr__(self):
        return "[%s:%s]"%(self.__class__.__name__,self.input_scwid)

class Revolution(DataAnalysis):
    input_revid=None

    scwver="001"

    def get_revid(self):
        return self.input_revid.handle

    def main(self):
        rbp=os.environ["REP_BASE_PROD"]
        self.revroot=rbp+"/scw/%s/"%self.get_revid()
        self.revdir=self.revroot+"/rev.001/"
        self.auxadppath=rbp+"/aux/adp/"+self.revid+"."+self.scwver

    def get_ijd(self):
        r1100=4306.5559396296
        r100=1315.4808007407

        r=int(self.get_revid())
        return r100+(r1100-r100)/1000*(r-100)

    
    def __repr__(self):
        return "[Revolution:%s]"%self.input_revid

class RevForScW(DataAnalysis):
    input_scw=ScWData    
    run_for_hashe=True
    allow_alias=False

    def __repr__(self):
        return "[RevForScW:for %s]"%repr(self.input_scw)

    def main(self):
        revid=self.input_scw.input_scwid.handle[:4]
        print "revolution id for scw:",revid
        return Revolution(input_revid=revid)

class Rev4ScW(Revolution):
    input_scw=ScWData    
    input_revid=dataanalysis.NoAnalysis

    def __repr__(self):
        return "[Rev4ScW:for %s]"%repr(self.input_scw)

    def get_revid(self):
        revid=self.input_scw.input_scwid.handle[:4]
        print "revolution id for scw:",revid
        return revid

class ICRoot(DataAnalysis):
    input="standard_IC"

    cached=False # level!

    schema_hidden=True
    version="v1"

    def main(self):
        self.icroot=os.environ['CURRENT_IC']
        self.icindex=self.icroot+"/idx/ic/ic_master_file.fits[1]"

        print('current IC:',self.icroot)


class IBIS_ICRoot(DataAnalysis):
    schema_hidden=True
    input="standard_IC"
    
    cached=False # level!

    def main(self):
        self.ibisicroot=os.environ['CURRENT_IC']+"/ic/ibis"
        print("current IBIS ic root is:"),self.ibisicroot

class GetLUT2(DataAnalysis):
    input="LUT2_standard"
    input_ibisic=IBIS_ICRoot

    cached=False

    def main(self):
        self.datafile=self.input_ibisic.ibisicroot+"/mod/isgr_3dl2_mod_0001.fits"
        #self.datafile=os.environ['REP_BASE_PROD']+"/ic/ibis/mod/isgr_3dl2_mod_0001.fits"
        #self.datafile=DataFile(os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_3dl2_mod_0001.fits")

class GetEcorrCalDB(DataAnalysis):
    input=["ecorr_standard_OSA10.2"]
    input_lut2=GetLUT2
    input_ibisic=IBIS_ICRoot
    input_scw=ScWData
    
    cached=False

    #ignore_input=["input"]


    def main(self):
        newest_ver="0001"
        newest_vstart=-1000

        t,dt=self.input_scw.get_t()

        for gmfile in sorted(glob.glob(self.input_ibisic.ibisicroot+"/mod/isgr_gain_mod_*.fits")): # by name?..
            ver=re.match(self.input_ibisic.ibisicroot+"/mod/isgr_gain_mod_(.*?).fits",gmfile).group(1)
            vstart=pyfits.open(gmfile)[1].header['VSTART']
        
            print vstart,ver
            if vstart>newest_vstart and t>vstart:
                newest_vstart=vstart
                newest_ver=ver

        self.godol=self.input_ibisic.ibisicroot+"/cal/ibis_isgr_gain_offset_0010.fits"
        self.supgdol=self.input_ibisic.ibisicroot+"/mod/isgr_gain_mod_"+newest_ver+".fits[ISGR-GAIN-MOD,1,BINTABLE]"
        self.supodol=self.input_ibisic.ibisicroot+"/mod/isgr_off2_mod_"+newest_ver+".fits[ISGR-OFF2-MOD,1,BINTABLE]"
        self.risedol=self.input_lut2.datafile


class ibis_isgr_energy_standard(DataAnalysis):
    cached=False
    
    input_scw=ScWData()
    #input_raw_events=None
    #input_ibis_hk=None
    # can also take only isgri events

    input_ecorrdata=GetEcorrCalDB

    version="v4_extras"
   
    def main(self):

        remove_withtemplate("isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)")
    
        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]" \
        ])

        bin="ibis_isgr_energy"
        ht=heatool(bin)
        ht['inGRP']="og.fits"
        ht['outCorEvts']="isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)"
        ht['useGTI']="n"
        ht['randSeed']=500
        ht['riseDOL']=self.input_ecorrdata.risedol
        ht['GODOL']=self.input_ecorrdata.godol
        ht['supGDOL']=self.input_ecorrdata.supgdol
        ht['supODOL']=self.input_ecorrdata.supodol
        ht['chatter']="4"
        ht.run()


        self.output_events=DataFile("isgri_events_corrected.fits")

class ibis_isgr_energy(DataAnalysis):
    #cached=True

    input_scw=ScWData()
    input_ecorrdata=GetEcorrCalDB

    version="v5_extras"

    binary=None

    def main(self):

        remove_withtemplate("isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)")

        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]", \
            self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]" \
        ])

        #bin=os.environ['COMMON_INTEGRAL_SOFTDIR']+"/spectral/ibis_isgr_energy/ibis_isgr_energy_102_pha2/ibis_isgr_energy"
        bin="ibis_isgr_energy"

        if self.binary is not None:
            bin=self.binary

        ht=heatool(bin)
        ht['inGRP']="og.fits"
        ht['outCorEvts']="isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)"
        ht['useGTI']="n"
        ht['randSeed']=500
        ht['riseDOL']=self.input_ecorrdata.risedol
        ht['GODOL']=self.input_ecorrdata.godol
        ht['supGDOL']=self.input_ecorrdata.supgdol
        ht['supODOL']=self.input_ecorrdata.supodol
        ht['chatter']="4"
        ht.run()

        self.output_events=DataFile("isgri_events_corrected.fits")
            

class ibis_isgr_evts_tag(DataAnalysis):
    cached=False
    
    input_events_corrected=ibis_isgr_energy
    input_scw=ScWData() # again, can get separate

    #input_ecorrdata=GetEcorrCalDB()

    version="v2"
   
    def main(self):

        #remove_withtemplate("isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)")

        cte="isgri_events_corrected_tagged.fits"
        
        shutil.copyfile(self.input_events_corrected.output_events.path,cte)
    
        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            cte \
        ])
            #self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]" \

        import_attr(self.input_scw.scwpath+"/swg.fits",['OBTSTART','OBTEND'])
        
        bin="ibis_isgr_evts_tag"
        ht=heatool(bin)
        ht['inGRP']="og.fits"
        ht['idxSwitch']=self.input_scw.revdirpath+"/idx/isgri_pxlswtch_index.fits[1]"
        ht['seleEXT']="ISGR-EVTS-COR"
        ht.run()

        self.output_events=DataFile(cte)

class UserGTI(DataAnalysis):
    pass

class gti_user(DataAnalysis):
    input_gti=UserGTI

    cached=True

    def main(self):
        fn="gti_user.fits"
        remove_withtemplate(fn)

        t1,t2=self.input_gti.gti

        bin="gti_user"
        ogc=heatool(bin)
        ogc['begin']=t1
        ogc['end']=t2
        ogc['gti']=fn
        ogc.run()

        self.gti=da.DataFile(fn)

# maybe split indeed,but try to show another case
class ibis_gti(DataAnalysis):
    input_scw=ScWData()
    input_ic=ICRoot()
    #input_gticreate=gti_create()

    cached=True
    
    version="v2"
    def main(self):
        # horrible horrible full OSA

        open("scw.list","w").write(self.input_scw.scwpath+"/swg.fits[1]")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))
        
        bin="og_create"
        ogc=heatool(bin)
        ogc['idxSwg']="scw.list"
        ogc['instrument']="IBIS"
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=os.getcwd().replace("[","_").replace("]","_") # dangerous
        ogc.run()
        
        scwroot="scw/"+self.input_scw.scwid

        bin="ibis_gti"
        ht=heatool(bin,wd="obs/"+ogc['ogid'].value) 
        ht['swgDOL']=scwroot+"/swg_ibis.fits"
        ht['GTI_Index']="ibis_gti.fits(IBIS-GNRL-GTI-IDX.tpl)"
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['disablePICsIT']="YES"
        ht['GTI_attTolerance_X']="0.05"
        ht['GTI_attTolerance_Z']="0.2"
        ht['GTI_BTI_Dol']=self.input_ic.icroot+"/ic/ibis/lim/isgr_gnrl_bti_0016.fits"
        ht['GTI_BTI_Names']="IBIS_CONFIGURATION IBIS_BOOT ISGRI_RISE_TIME VETO_PROBLEM SOLAR_FLARE BELT_CROSSING"

        if hasattr(self,'input_usergti'):
            path=self.input_usergti.gti.get_path()
            if os.path.abspath(path)==os.path.normpath(path):
                print "full path",path
            else:
                print "not a full path",path
                path="../../"+path
            ht['GTI_gtiUserI']=path
            ht['GTI_TimeFormat']='UTC'


        ht.run()

        shutil.copy(ht.cwd+"/ibis_gti.fits","./ibis_gti.fits")
        self.output_gti=DataFile("ibis_gti.fits")

        gti=pyfits.open("ibis_gti.fits")[-1].data
        print gti
        

# maybe split indeed,but try to show another way
class ibis_dead(DataAnalysis):
    input_scw=ScWData()
    input_ic=ICRoot()

    cached=True
    
    version="v2"
    def main(self):
        # horrible horrible full OSA

        #open("scw.list","w").write(self.input_scw.scwid)
        open("scw.list","w").write(self.input_scw.scwpath+"/swg.fits[1]")

        print("scw path:",self.input_scw.scwpath+"/swg.fits[1]")
        if not os.path.exists(self.input_scw.scwpath+"/swg.fits"):
            raise Exception("no scw? broken!")

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))
        
        wd=os.getcwd().replace("[","_").replace("]","_")
        bin="og_create"
        ogc=heatool(bin)
        ogc['idxSwg']="scw.list"
        ogc['instrument']="IBIS"
        ogc['ogid']="scw_"+self.input_scw.scwid
        ogc['baseDir']=wd # dangerous
        ogc.run()
        
        scwroot="scw/"+self.input_scw.scwid

        bin="ibis_dead"
        ht=heatool(bin,wd=wd+"/obs/"+ogc['ogid'].value) 
        ht['swgDOL']=scwroot+"/swg_ibis.fits"
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['isgroutDead']="isgri_dead.fits(ISGR-DEAD-SCP.tpl)"
        ht['picsoutDead']="picsit_dead.fits(PICS-DEAD-SCP.tpl)"
        ht['compoutDead']="compton_dead.fits(COMP-DEAD-SCP.tpl)"
        ht['disablePICsIT']="YES"
        ht['disableCompton']="YES"
        ht.run()

        shutil.copy(ht.cwd+"/isgri_dead.fits","./isgri_dead.fits")
        self.output_dead=DataFile("isgri_dead.fits")

cache_local=MemCacheIntegralFallback()
cache_local.filecacheroot=os.getcwd()
class ISGRIEvents(DataAnalysis):
    input_evttag=ibis_isgr_evts_tag
    
    cached=True
    cache=cache_local

    read_caches=[cache_local.__class__]
    write_caches=[cache_local.__class__]


    version="v3"
    def main(self):
        self.events=self.input_evttag.output_events

class ImageBins(DataAnalysis):
    input_binsname="g25-80"
    ebins=None

    rmfbins=False

    def main(self):
        if self.ebins is None:
            self.bins=[(25,80)]
        else:
            self.bins=self.ebins

class SpectraBins(DataAnalysis):
    input_binsname="spectral_bins_62"
    
    rmfbins=True

    version="v3"
    def main(self):
        self.binrmf=os.environ['INTEGRAL_DATA']+"/resources/rmf_62bands.fits" # noo!!!
        e=pyfits.open(self.binrmf)[3].data
        self.bins=zip(e['E_MIN'],e['E_MAX'])
        self.binrmfext=self.binrmf+'[3]'

    def get_binrmfext(self):
        return self.binrmfext

class ListBins(DataAnalysis):
    input_bins=ImageBins

    def main(self):
        open("f.txt","w").write(str(self.input_bins.bins))

class BinEventsVirtual(DataAnalysis):
    input_scw=ScWData

    input_events=ISGRIEvents
    input_gti=ibis_gti
    input_dead=ibis_dead

    target_level=None
    input_bins=None

    maxrisetime=116
    minrisetime=16

    version="v2"
    
    cached=True

    default_log_level="binevents"

    ii_shadow_build_binary="ii_shadow_build"


    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.maxrisetime!=116:
            v+=".lrt%i"%self.maxrisetime
        if self.minrisetime!=16:
            v+=".hrt%i"%self.minrisetime
        return v

    def main(self):
        if self.target_level is None or self.input_bins is None:
            raise Exception("VirtualAnalysis: please inherit!")
        
        self.pre_process()

        # ask stephane why need raw events
        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_events.events.get_path(), \
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]", \
            self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]" \
        ]) # get separately tc etc

        import_attr(self.input_scw.scwpath+"/swg.fits",['OBTSTART','OBTEND'])
        set_attr({'ISDCLEVL':self.target_level})

        det_fn="isgri_detector_shadowgram_%s.fits"%self.target_level
        det_tpl="(ISGR-DETE-SHD-IDX.tpl)"
        eff_fn="isgri_efficiency_shadowgram_%s.fits"%self.target_level
        eff_tpl="(ISGR-EFFI-SHD-IDX.tpl)"

        remove_withtemplate(det_fn+det_tpl)
        remove_withtemplate(eff_fn+eff_tpl)

        bin=self.ii_shadow_build_binary
        ht=heatool(bin)
        ht['outSWGGRP']="og.fits[GROUPING,1,BINTABLE]"
        ht['inDead']=self.input_dead.output_dead.get_path()
        ht['inGTI']=self.input_gti.output_gti.get_path()
        ht['gti_name'] = 'MERGED_ISGRI'
        ht['outputLevel'] = self.target_level

        print "target_level",self.target_level

        print "has rmfbins "+str(self.input_bins.rmfbins) if hasattr(self.input_bins,'rmfbins') else "no rmfbins"
        print "has binrmfext" if hasattr(self.input_bins,'binrmfext') else "no binrmfext"

        if ( self.target_level=="BIN_I" or not hasattr(self.input_bins,'rmfbins') or not self.input_bins.rmfbins or not hasattr(self.input_bins,'binrmfext') ) and not ( hasattr(self.input_bins,'rmfbins') and self.input_bins.rmfbins ): # fix!!
            ht['isgri_e_num'] = len(self.input_bins.bins)
            ht['isgri_e_min'] = " ".join([str(a[0]) for a in self.input_bins.bins])
            ht['isgri_e_max'] = " ".join([str(a[1]) for a in self.input_bins.bins])
        elif self.target_level=="BIN_S" or ( hasattr(self.input_bins,'rmfbins') and self.input_bins.rmfbins ) or hasattr(self.input_bins,'binrmfext'):
            ht['isgri_e_num'] = -1

            rmf=self.input_bins.binrmfext
            if isinstance(rmf,DataFile):
                ht['inEnergyValues'] = rmf.get_path()  #  +"[1]" will this work?
            else:
                ht['inEnergyValues'] = rmf #  +"[1]" will this work?
        else:
            raise Exception("neither bins given!")

        ht['isgri_min_rise'] = self.minrisetime
        ht['isgri_max_rise'] = self.maxrisetime
        ht['isgri_t_len'] = 10000000
        ht['idxLowThre']=self.input_scw.revdirpath+"/idx/isgri_context_index.fits[1]"
        ht['idxNoisy']=self.input_scw.revdirpath+"/idx/isgri_prp_noise_index.fits[1]"
        ht['outRawShadow']=det_fn+det_tpl
        ht['outEffShadow']=eff_fn+eff_tpl

        self.extra_pars(ht)

        ht.run()

        self.shadow_detector=DataFile(det_fn)
        self.shadow_efficiency=DataFile(eff_fn)
        
        self.post_process()

    def extra_pars(selt,ht):
        pass

    def post_process(self):
        pass
    
    def pre_process(self):
        pass


class BinEventsImage(BinEventsVirtual):
    target_level="BIN_I"
    input_bins=ImageBins

class BinEventsSpectra(BinEventsVirtual):
    target_level="BIN_S"
    input_bins=SpectraBins

class BinMapsVirtual(DataAnalysis):
    input_bins=ImageBins
    input_ic=IBIS_ICRoot
    # and dependency on bkg!
    
    target_level=None
    input_bins=None
    
    cached=True

    version="v2"
    def main(self):
        #construct_gnrl_scwg_grp(self.input_scw,[\
        #        ])
        if self.target_level is None or self.input_bins is None:
            raise Exception("VirtualAnalysis: please inherit!")

        construct_empty_shadidx(self.input_bins.bins,levl=self.target_level)


        maps={
                'back':('Bkg',self.input_ic.ibisicroot+"/"+'bkg/isgr_back_bkg_0007.fits[1]'),
                'corr':('Corr',self.input_ic.ibisicroot+"/"+'mod/isgr_effi_mod_0011.fits[1]'),
                'unif':('Uni',self.input_ic.ibisicroot+"/"+'bkg/isgr_unif_bkg_0002.fits[1]')
                }

        if hasattr(self,'input_unif'):
            maps['unif']=('Uni',self.input_unif.unif.get_path()+"[1]")
            print "will use uniformity:",maps['unif']
        
        if hasattr(self,'input_bkg'):
            maps['back']=('Bkg',self.input_bkg.bkg.get_path()+"[1]")
            print "will use background:",maps['back']

        level2key={
                'BIN_I':'ima',
                'BIN_S':'spe'
                }
        
        lk=level2key[self.target_level]
        
        for k,(k2,m) in maps.items():
            fn="rebinned_"+k+"_"+lk+".fits"

            remove_withtemplate(fn) # again

            bin="ii_map_rebin"
            ht=heatool(bin)
            ht['outSwg']='og.fits[GROUPING,1,BINTABLE]'
            ht['OutType']=self.target_level
            ht['slope']='-2'
            ht['arfDol']=self.input_ic.ibisicroot+'/mod/isgr_effi_mod_0011.fits[ISGR-ARF.-RSP]'
            ht['inp'+k2+'Dol']=m 
            ht['reb'+k2+'Dol']=fn
            ht.run()

            setattr(self,k,DataFile(fn))

class BinMapsImage(BinMapsVirtual):
    target_level="BIN_I"
    input_bins=ImageBins

class BinMapsSpectra(BinMapsVirtual):
    target_level="BIN_S"
    input_bins=SpectraBins


import subprocess,os
import dataanalysis as da
import ddosa
import ast
from astropy import wcs as pywcs

class GRcat(DataAnalysis):
    input="gnrl_ref_cat_40"
    input_ic=ICRoot
    
    suffix=None

    cached=False # again, this is transient-level cache

    def get_version(self):
        v=self.get_signature()+"."+self.version
        if self.suffix is not None:
            v=v+"."+self.suffix
        return v

    def main(self):
        if self.suffix is None:
            self.cat=os.environ['REP_BASE_PROD']+"/cat/hec/gnrl_refr_cat_0040.fits[1]"
        else:
            self.cat=os.environ['REP_BASE_PROD']+"/cat/hec/gnrl_refr_cat_0040_%s.fits[1]"%self.suffix


class BrightCat(DataAnalysis):
    input=GRcat
    #input_selection="flag5"
    input_selection="flag5andover100"

    cached=False

    def main(self):
        #self.cat=self.input.cat+"[ISGRI_FLAG2==5]"
        self.cat=self.input.cat+"[ISGRI_FLAG2==5&&ISGR_FLUX_1>100]"

class BrightPIFImage(DataAnalysis):
    input_scw=ScWData
    input_cat=BrightCat#CatExtract(assume=ddosa.ISGRIRefCat(input=GRcat38m))
    input_ic=IBIS_ICRoot
    input_bins=ImageBins
    input_gti=ibis_gti

    copy_cached_input=True


    def main(self):
        if isinstance(self.input_cat.cat,str):
            catfn=self.input_cat.cat
        else:
            catfn=self.input_cat.cat.get_path()

        construct_gnrl_scwg_grp(self.input_scw,[\
                    catfn,
                    self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]",
                    self.input_gti.output_gti.path
                ])

        import_attr(self.input_scw.scwpath+"/swg.fits",["OBTSTART","OBTEND","TSTART","TSTOP","SW_TYPE","TELAPSE"])
        set_attr({'ISDCLEVL':"BIN_I"})
        set_attr({'INSTRUME':"IBIS"},"og.fits")

        construct_gnrl_scwg_grp_idx([
                    "og.fits"]
                )
        set_attr({'ISDCLEVL':"BIN_I"},"og_idx.fits")

        construct_og([\
                    "og_idx.fits",
                ])
        set_attr({'ISDCLEVL':"BIN_I"},"ogg.fits")


        remove_withtemplate("isgri_model.fits(ISGR-PIF.-SHD.tpl)")



        #ht=heatool(os.environ['COMMON_INTEGRAL_SOFTDIR']+"/ii_pif/ii_pif_oof/ii_pif")
        ht=heatool("ii_pif")
        ht['inOG']=""
        ht['outOG']="ogg.fits[1]"
        ht['inCat']=catfn
        ht['mask']=self.input_ic.ibisicroot+"/mod/isgr_mask_mod_0003.fits[ISGR-MASK-MOD,1,IMAGE]"
#        ht['deco']=self.input_ic.ibisicroot+"/mod/isgr_deco_mod_0008.fits[ISGR-DECO-MOD,1,IMAGE]"
        ht['tungAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0010.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['aluAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0011.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['leadAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0012.fits[ISGR-ATTN-MOD,1,BINTABLE]"
 #       ht['covrMod']=self.input_ic.ibisicroot+"/mod/isgr_covr_mod_0002.fits[1]"

        ht['num_band'] = len(self.input_bins.bins)
        ht['E_band_min'] = " ".join([str(a[0]) for a in self.input_bins.bins])
        ht['E_band_max'] = " ".join([str(a[1]) for a in self.input_bins.bins])
        
        #ht['AllowOffEdge'] = 100
        ht.run()

        self.pifs=da.DataFile("isgri_model.fits")
       # self.skyima=DataFile("isgri_sky_ima.fits")
       # self.skyres=DataFile("isgri_sky_res.fits")

class BrightPIFSpectra(BrightPIFImage):
    input_scw=ScWData
    input_cat=BrightCat#CatExtract(assume=ddosa.ISGRIRefCat(input=GRcat38m))
    input_ic=IBIS_ICRoot
    input_bins=SpectraBins
    input_gti=ibis_gti


class ShadowUBCVirtual(DataAnalysis):
    input_scw=ScWData
    level=None
    input_shadows=None
    input_maps=None
    #input_brpif=None #BrightPIFImage

    version="v3"

    brPifThreshold=1e-4
    
    cached=False

    binary="ii_shadow_ubc"

    def get_version(self):
        return self.get_signature()+"."+self.version+(".brthr%.5lg"%self.brPifThreshold if self.brPifThreshold!=1e-4 else "")

    def main(self):
        construct_gnrl_scwg_grp(self.input_scw,[\
                self.input_shadows.shadow_detector.get_path(),
                self.input_shadows.shadow_efficiency.get_path(),
            ])
        
        keys=["SWID","SW_TYPE","REVOL","EXPID","OBTSTART","OBTEND","TSTART","TSTOP","SWBOUND"]
        import_attr(self.input_scw.scwpath+"/swg.fits",keys)

        fn,tpl="isgri_cor_shad_%s.fits"%self.level,"(ISGR-CEXP-SHD-IDX.tpl)"
        remove_withtemplate(fn+tpl)
        
        ht=heatool(self.binary)
        #ht=heatool("/workdir/lin/soft/ii_shadow_ubc_bestim/ii_shadow_ubc")
        ht['outSWGRP']="og.fits"
        ht['OutType']=self.level
        ht['outCorShadow']=fn+tpl
        ht['isgrUnifDol']=self.input_maps.unif.get_path()
        ht['isgrBkgDol']=self.input_maps.back.get_path()
        if hasattr(self,'input_brpif'):
            ht['brPif']=self.input_brpif.pifs.get_path()
            ht['brPifThreshold']=self.brPifThreshold
        ht['method_cor']=1 # keep separately background correction routines
        ht.run()

        self.corshad=DataFile(fn)

        self.post_process()
    
    def post_process(self):
        pass

class ShadowUBCImage(ShadowUBCVirtual):
    level="BIN_I"
    input_shadows=BinEventsImage
    input_maps=BinMapsImage
    input_brpif=BrightPIFImage

class ShadowUBCSpectra(ShadowUBCVirtual):
    level="BIN_S"
    input_shadows=BinEventsSpectra
    input_maps=BinMapsSpectra
    #input_brpif=BrightPIFSpectra


class GBcat(DataAnalysis):
    input=GRcat
    input_selection="flag5"

    cached=False

    def main(self):
        self.cat=self.input.cat+"[ISGRI_FLAG2==5]"

class ghost_bustersVirtual(DataAnalysis):
    input_shadow=None
    level=None

    input_ic=IBIS_ICRoot
    input_scw=ScWData()
    input_cat=GBcat

    version="v2"
    
    cached=True

    gb_binary=None

    def main(self):

        att=self.input_scw.auxadppath+"/attitude_historic.fits"
        if os.path.exists(att):
            att=self.input_scw.auxadppath+"/attitude_historic.fits[AUXL-ATTI-HIS,1,BINTABLE]"
            attp=att
        else:
            att=self.input_scw.auxadppath+"/attitude_snapshot.fits[AUXL-ATTI-SNA,1,BINTABLE]"
            attp_g=glob.glob(self.input_scw.auxadppath+"/attitude_predicted_*.fits*")
            attp=attp_g[0]+"[AUXL-ATTI-PRE,1,BINTABLE]"

        construct_gnrl_scwg_grp(self.input_scw,[\
                self.input_shadow.corshad.path,
                att, \
                attp \
            ])

        import_attr(self.input_scw.scwpath+"/swg.fits",["TSTART","TSTOP"])
        
        if self.gb_binary is not None:
            ht=heatool(self.gb_binary)
        else:
            ht=heatool("ghost_busters")
        ht['ogDOL']=""
        ht['sourcecat']=self.input_cat.cat
        ht['maskmod']=self.input_ic.ibisicroot+"/mod/isgr_ghos_mod_001.fits[ISGR-GHOS-MOD,1,IMAGE]"
        ht['inDOL']="og.fits"
        ht.run()

        r="isgri_cor_shad_%s_gb.fits"%self.level
        shutil.copyfile(self.input_shadow.corshad.path,r)

        self.corshad=DataFile(r)
 
class ghost_bustersImage(ghost_bustersVirtual):
    input_shadow=ShadowUBCImage
    level="BIN_I"

class ghost_bustersSpectra(ghost_bustersVirtual):
    input_shadow=ShadowUBCSpectra
    level="BIN_S"

class ISGRIRefCat(DataAnalysis):
    input=GRcat
    input_selection="onlyisgri33"

    cached=False # since the path is transient

    def main(self):
        self.cat=self.input.cat+"[ISGRI_FLAG==1 || ISGRI_FLAG==2]"
        import time
        time.sleep(3) # expensive to search for cat: test

class CatExtractEmpty(DataAnalysis):
    input_cat=GRcat
    input_scw=ScWData
    
    cached=True

    def main(self):
        construct_gnrl_scwg_grp(self.input_scw,[\
                ])
        import_attr(self.input_scw.scwpath+"/swg.fits",["RA_SCX","DEC_SCX"])

        remove_withtemplate("isgri_catalog.fits(ISGR-SRCL-CAT.tpl)")

        ht=heatool("cat_extract")
        ht['outGRP']="og.fits[1]"
        ht['outCat']="isgri_catalog.fits(ISGR-SRCL-CAT.tpl)"
        ht['outExt']='ISGR-SRCL-CAT'
        ht['instrument']='ISGRI'
        ht['clobber']='yes'
        ht['refCat']=self.input_cat.cat+"[ISGRI_FLAG==-1]"
        ht.run()

        self.cat=DataFile("isgri_catalog.fits")

class CatExtract(DataAnalysis):
    input_cat=ISGRIRefCat
    input_scw=ScWData
    
    cached=True

    def main(self):
        construct_gnrl_scwg_grp(self.input_scw,[\
                ])
        import_attr(self.input_scw.scwpath+"/swg.fits",["RA_SCX","DEC_SCX"])

        remove_withtemplate("isgri_catalog.fits(ISGR-SRCL-CAT.tpl)")

        ht=heatool("cat_extract")
        ht['outGRP']="og.fits[1]"
        ht['outCat']="isgri_catalog.fits(ISGR-SRCL-CAT.tpl)"
        ht['outExt']='ISGR-SRCL-CAT'
        ht['instrument']='ISGRI'
        ht['clobber']='yes'
        ht['refCat']=self.input_cat.cat
        ht.run()

        self.cat=DataFile("isgri_catalog.fits")

class ImagingConfig(DataAnalysis):
    input="onesource_negmod"

    SearchMode=3
    ToSearch=5
    CleanMode=1
    MinCatSouSnr=4
    MinNewSouSnr=5
    NegModels=1
    DoPart2=0
    SouFit=0

class ii_skyimage(DataAnalysis):
    input_gb=ghost_bustersImage
    input_maps=BinMapsImage
    input_bins=ImageBins
    input_cat=CatExtract
    input_ic=IBIS_ICRoot
    input_imgconfig=ImagingConfig
    input_scw=ScWData
    input_gti=ibis_gti

    cached=True

    ii_skyimage_binary=None

    save_image=True

    image_tag=None

    version="v2"

    outtype="BIN_I"

    def get_version(self):
        v=self.get_signature()+"."+self.version
        for k in ['SouFit','SearchMode','ToSearch','CleanMode','MinCatSouSnr','MinNewSouSnr','NegModels','DoPart2']: # dopart2 is flow control, separately
            if hasattr(self,'ii_'+k):
                v+="_"+k+"_"+str(getattr(self,'ii_'+k))
        return v

    def main(self):
        construct_gnrl_scwg_grp(self.input_scw,[\
                    self.input_gb.corshad.path,
                    self.input_cat.cat.path,
                    self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]",
                    self.input_gti.output_gti.path,
                ])
        
        import_attr(self.input_scw.scwpath+"/swg.fits",["OBTSTART","OBTEND","TSTART","TSTOP","SW_TYPE","TELAPSE","SWID","SWBOUND"])
        set_attr({'ISDCLEVL':self.outtype})
        set_attr({'INSTRUME':"IBIS"},"og.fits")

        construct_gnrl_scwg_grp_idx([\
                    "og.fits",
                ])
        set_attr({'ISDCLEVL':self.outtype},"og_idx.fits")
        
        construct_og([\
                    "og_idx.fits",
                ])
        set_attr({'ISDCLEVL':self.outtype},"ogg.fits")
        
        remove_withtemplate("isgri_srcl_res.fits(ISGR-SRCL-RES.tpl)")
        remove_withtemplate("isgri_mosa_ima.fits(ISGR-MOSA-IMA-IDX.tpl)")
        remove_withtemplate("isgri_mosa_res.fits(ISGR-MOSA-RES-IDX.tpl)")
        remove_withtemplate("isgri_sky_ima.fits")
        remove_withtemplate("isgri_sky_res.fits")

        if self.ii_skyimage_binary is None:
            ii_skyimage_binary="ii_skyimage"
        else:
            ii_skyimage_binary=self.ii_skyimage_binary

        ht=heatool(ii_skyimage_binary)
        ht['outOG']="ogg.fits[1]"
        ht['outCat']="isgri_srcl_res.fits(ISGR-SRCL-RES.tpl)"
        ht['mask']=self.input_ic.ibisicroot+"/mod/isgr_mask_mod_0003.fits[ISGR-MASK-MOD,1,IMAGE]"
        ht['deco']=self.input_ic.ibisicroot+"/mod/isgr_deco_mod_0008.fits[ISGR-DECO-MOD,1,IMAGE]"
        ht['tungAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0010.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['aluAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0011.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['leadAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0012.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['covrMod']=self.input_ic.ibisicroot+"/mod/isgr_covr_mod_0002.fits[1]"
        ht['outMosIma']="isgri_mosa_ima.fits(ISGR-MOSA-IMA-IDX.tpl)"
        ht['outMosRes']="isgri_mosa_res.fits(ISGR-MOSA-RES-IDX.tpl)"
        ht['ScwDir'] = './'
        ht['ScwType'] = 'pointing'
        ht['ExtenType'] = 2
        ht['OutType']=self.outtype
        ht['num_band'] = len(self.input_bins.bins)
        ht['E_band_min'] = " ".join([str(a[0]) for a in self.input_bins.bins])
        ht['E_band_max'] = " ".join([str(a[1]) for a in self.input_bins.bins])
        for k in ['SouFit','SearchMode','ToSearch','CleanMode','MinCatSouSnr','MinNewSouSnr','NegModels','DoPart2']: # dopart2 is flow control, separately
            ht[k]=getattr(self.input_imgconfig,k)
            if hasattr(self,'ii_'+k): ht[k]=getattr(self,'ii_'+k)
        ht['corrDol'] = self.input_maps.corr.path
        ht.run()


        if not os.path.exists("isgri_sky_ima.fits"):
            print "no image produced: since there was no exception in the binary, assuming empty results"
            self.empty_results=True
            return

        self.srclres=DataFile("isgri_srcl_res.fits")

        if self.save_image:
            self.skyima=DataFile("isgri_sky_ima.fits")
        self.skyres=DataFile("isgri_sky_res.fits")

        if self.image_tag is not None:      
            try:
                tag=str(self.image_tag)
                shutil.copyfile("isgri_sky_res.fits","isgri_sky_res_%s.fits"%tag)
                if self.save_image:
                    shutil.copyfile("isgri_sky_ima.fits","isgri_sky_ima_%s.fits"%tag)
            except:
                print "something went wrong"
    
        self.post_process()

    def post_process(self):
        pass
            

class CatForSpectraFromImaging(DataAnalysis):
    input_imaging=ii_skyimage

    minsig=None
    maxsources=None

    def get_version(self):
        return self.get_signature()+"."+self.version+("" if self.minsig is None else "minsig%.3lg"%self.minsig)

    def main(self):
        if hasattr(self.input_imaging,'empty_results'):
            print "no results here"
            self.empty_results=True
            return

        catfn="cat4spectra.fits"

        f=pyfits.open(self.input_imaging.srclres.path)

        print "image catalog contains",len(f[1].data)

        if self.minsig is not None:
            f[1].data=f[1].data[f[1].data['DETSIG']>self.minsig]
            print "selecting by significance",len(f[1].data)
        
        if self.maxsources is not None:
            raise Exception("not implemented")

        f.writeto(catfn,clobber=True)

        self.cat=DataFile(catfn)

#class CatForSpectra(da.DataAnalysis):
#    pass

class ISGRIResponse(DataAnalysis):
    input_ecorrdata=GetEcorrCalDB

    path="/sps/integral/analysis/savchenk/lut2tests/test_05351/single_scw/resources/rmf/rmf_62bands.fits"

class ii_spectra_extract(DataAnalysis):
    input_gb=ghost_bustersSpectra
    input_cat=CatForSpectraFromImaging # or virtual
    input_ic=IBIS_ICRoot
    input_response=ISGRIResponse
    input_scw=ScWData()
    input_maps=BinMapsSpectra

    input_gti=ibis_gti
    
    cached=True

    version="v3"

    #input_bins=SpectraBins
    #input_cat=CatExtract
    #input_imgconfig=ImagingConfig

    shdtype="BIN_S"
    binary="ii_spectra_extract"

    usebkg=True

    report_runtime_destination="mysql://pixels.runtime"

    def main(self):
        if hasattr(self.input_cat,'empty_results'):
            print "empty here"
            self.empty_results=True
            return

        att=self.input_scw.auxadppath+"/attitude_historic.fits"
        if os.path.exists(att):
            att=self.input_scw.auxadppath+"/attitude_historic.fits[AUXL-ATTI-HIS,1,BINTABLE]"
            attp=att
        else:
            att=self.input_scw.auxadppath+"/attitude_snapshot.fits[AUXL-ATTI-SNA,1,BINTABLE]"
            attp_fn=glob.glob(self.input_scw.auxadppath+"/attitude_predicted_*.fits*")[0]
            attp=attp_fn+"[AUXL-ATTI-PRE,1,BINTABLE]"
            

        construct_gnrl_scwg_grp(self.input_scw,[\
                    self.input_gb.corshad.path,
                    self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]",
                    att,
                    attp,
                    self.input_gti.output_gti.path
                ])
                    #self.input_cat.cat.path,
        
        import_attr(self.input_scw.scwpath+"/swg.fits",["OBTSTART","OBTEND","TSTART","TSTOP","SW_TYPE","TELAPSE","SWID"])
        set_attr({'ISDCLEVL':"BIN_S"})
        #set_attr({'INSTRUME':"IBIS"},"og.fits")

        #construct_gnrl_scwg_grp_idx(self.input_scw,[\
        #            "og.fits",
        #        ])
        #set_attr({'ISDCLEVL':"BIN_I"},"og_idx.fits")
        
      #  construct_og(self.input_scw,[\
      #              "og_idx.fits",
      #          ])
      #  set_attr({'ISDCLEVL':"BIN_I"},"ogg.fits")
        
        #remove_withtemplate("isgri_srcl_res.fits(ISGR-SRCL-RES.tpl)")

        pif_fn,pif_tpl="isgri_pif.fits","(ISGR-PIF.-SHD-IDX.tpl)"
        spec_fn,spec_tpl="isgri_spectrum.fits","(ISGR-EVTS-SPE-IDX.tpl)"

        remove_withtemplate(pif_fn+pif_tpl)
        remove_withtemplate(spec_fn+spec_tpl)

        ht=heatool(self.binary)
        ht['outSwg']="og.fits"
        ht['inCat']=self.input_cat.cat.get_path()
        ht['outPif']=pif_fn+pif_tpl
        ht['outSpec']=spec_fn+spec_tpl
        ht['mask']=self.input_ic.ibisicroot+"/mod/isgr_mask_mod_0003.fits[ISGR-MASK-MOD,1,IMAGE]"
        ht['tungAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0010.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['aluAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0011.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['leadAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0012.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['idx_isgrResp']=self.input_response.path
        ht['isgrUnifDol']=self.input_maps.unif.path
        if self.usebkg:
            ht['isgrBkgDol']=self.input_maps.back.path
        else:
            ht['isgrBkgDol']="-"
        ht['corrDol']=self.input_maps.corr.path
        ht['OutType']=self.shdtype
        ht['method_cor']=1
        
        if hasattr(self,'input_bins') and not self.input_bins.rmfbins:
            ebins=self.input_bins.bins
            ht['num_band']=len(ebins)
            ht['E_band_min'],ht['E_band_max']=[" ".join(["%.5lg"%b for b in a]) for a in zip(*ebins)]
        else:
            rmf=self.input_bins.get_binrmfext() if hasattr(self,'input_bins') else self.input_response.path
            ht['num_band'] = -1
            ht['idx_isgrResp'] = rmf #  +"[1]" will this work?


        #for k in ['SearchMode','ToSearch','CleanMode','MinCatSouSnr','MinNewSouSnr','NegModels','DoPart2']: # dopart2 is flow control, separately
        #    ht[k]=getattr(self.input_imgconfig,k)

        ht.run()

        self.spectrum=DataFile(spec_fn)
        self.pifs=DataFile(pif_fn)


class root(DataAnalysis):
    input=[ibis_gti()]

#
#
#### tools
#
#

def remove_withtemplate(fn):
    s=re.search("(.*?)\((.*?)\)",fn)
    if s is not None:
        try:
            os.remove(s.group(2))
        except OSError:
            pass
        fn=s.group(1)

    try:
        os.remove(fn)
    except OSError:
        pass
    
    try:
        os.remove(fn+".gz")
    except OSError:
        pass
    

def construct_gnrl_scwg_grp(scw,children=[],fn="og.fits"):
    try:
        os.remove(fn)
    except OSError:
        pass

    dc=heatool("dal_create")
    dc['obj_name']="!"+fn
    dc['template']="GNRL-SCWG-GRP.tpl"
    dc.run()
    
    da=heatool("dal_attr")
    da['indol']=fn
    da['keynam']="REVOL"
    da['action']="WRITE"
    da['type']="DAL_INT"
    da['value_i']=scw.revid
    da.run()
    
    if children!=[]:
        if len(children)<=4:
            dac=heatool("dal_attach")
            dac['Parent']=fn

            for i in range(4):
                dac['Child%i'%(i+1)]=""

            for i,c in enumerate(children):
                dac['Child%i'%(i+1)]=c
                if i>3:
                    raise Exception("can not attach more than 4 children to the group!")
            dac.run()
        else:
            for children_group_id in range(len(children)/4+1):
                children_group=children[children_group_id*4:children_group_id*4+4]
                dac=heatool("dal_attach")
                dac['Parent']=fn

                for i in range(4):
                    dac['Child%i'%(i+1)]=""

                for i,c in enumerate(children_group):
                    dac['Child%i'%(i+1)]=c
                    if i>3:
                        raise Exception("can not attach more than 4 children to the group!")
                dac.run()

def construct_gnrl_scwg_grp_idx(children=[],fn="og_idx.fits"):
    remove_withtemplate(fn+"(GNRL-SCWG-GRP-IDX.tpl)")

    open("swgs.txt","w").write("\n".join(children))
    dc=heatool("txt2idx")
    dc['element']="swgs.txt"
    dc['index']=fn
    dc['template']="GNRL-SCWG-GRP-IDX.tpl"
    dc.run()

def construct_og(children=[],fn="ogg.fits"):
    remove_withtemplate(fn+"(GNRL-OBSG-GRP.tpl)")

    dc=heatool("dal_create")
    dc['obj_name']="!"+fn
    dc['template']="GNRL-OBSG-GRP.tpl"
    dc.run()
    
    if children!=[]:
        dac=heatool("dal_attach")
        dac['Parent']=fn
        for i,c in enumerate(children):
            dac['Child%i'%(i+1)]=c
            if i>3:
                raise Exception("can not attach more than 4 children to the group!")
        dac.run()

def import_attr(obj,attr,target="og.fits"):
    da=heatool("dal_attr_copy")
    da['indol']=obj
    da['outdol']=target
    da['keylist']=",".join(attr)
    da.run()
    
def set_attr(attrs,fn="og.fits"):

    pt2dt={int:"DAL_INT",str:"DAL_CHAR"}
    pt2k={int:"i",str:"s"}

    for k,v in attrs.items():
        da=heatool("dal_attr")
        da['indol']=fn
        da['keynam']=k
        da['action']="WRITE"
        da['type']=pt2dt[type(v)]
        da['value_'+pt2k[type(v)]]=v
        da.run()
        
def construct_empty_shadidx_old(bins,fn="og.fits",levl="BIN_I"):
    remove_withtemplate(fn+"(ISGR-DETE-SHD-IDX.tpl)")

    ht=heatool("dal_create")
    ht['obj_name']=fn
    ht['template']="ISGR-DETE-SHD-IDX.tpl"
    ht.run()

    for e1,e2 in bins:
        tshad="shad_%.5lg_%.5lg.fits"%(e1,e2)
        remove_withtemplate(tshad)

        ht=heatool("dal_create")
        ht['obj_name']=tshad
        ht['template']="ISGR-DETE-SHD.tpl"
        ht.run()

        da=heatool("dal_attr")
        da['indol']=ht['obj_name'].value
        da['keynam']="E_MIN"
        da['action']="WRITE"
        da['type']="DAL_DOUBLE"
        da['value_r']=e1
        da.run()

        da=heatool("dal_attr")
        da['indol']=ht['obj_name'].value
        da['keynam']="E_MAX"
        da['action']="WRITE"
        da['type']="DAL_DOUBLE"
        da['value_r']=e2
        da.run()
        
        da=heatool("dal_attr")
        da['indol']=ht['obj_name'].value
        da['keynam']="ISDCLEVL"
        da['action']="WRITE"
        da['type']="DAL_CHAR"
        da['value_s']="BIN_I"
        da.run()

        da=heatool("dal_attach")
        da['Parent']=fn
        da['Child1']=ht['obj_name'].value
        da.run()

    # attaching does not create necessary fields: use txt2idx instead
    
    og=pyfits.open(fn) 
    for i,(e1,e2) in enumerate(bins):
        og[1].data[i]['E_MIN']=e1
        og[1].data[i]['E_MAX']=e2
        og[1].data[i]['ISDCLEVL']=levl
    og.writeto(fn,clobber=True)

def construct_empty_shadidx(bins,fn="og.fits",levl="BIN_I"):
    remove_withtemplate(fn+"(ISGR-DETE-SHD-IDX.tpl)")


    ht=heatool("dal_create")
    ht['obj_name']=fn
    ht['template']="ISGR-DETE-SHD-IDX.tpl"
    ht.run()

    for e1,e2 in bins:
        print "appending",e1,e2
        dc=heatool('dal_append')
        dc['grpDOL']=fn
        dc['element']='ISGR-DETE-SHD.tpl'
        dc.run()

    og=pyfits.open(fn) 
    for i,(e1,e2) in enumerate(bins):
        og[1].data[i]['E_MIN']=e1
        og[1].data[i]['E_MAX']=e2
        og[1].data[i]['ISDCLEVL']=levl
        og[2+i].header['E_MIN']=e1
        og[2+i].header['E_MAX']=e2
        og[2+i].header['ISDCLEVL']=levl
    og.writeto(fn,clobber=True)

class AnyScW(dataanalysis.AnyAnalysis):
    pass

class AnyRev(dataanalysis.AnyAnalysis):
    pass

class BinnedDataProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=BinEventsImage(assume=ScWData(input_scwid="any",use_abstract=True)) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
       # print "one scw hash:",ahash
        #ahash=dataanalysis.hashe_replace_object(ahash,'AnyScW','None')
        print "generalized hash:",ahash
        rh=dataanalysis.shhash(ahash)
        print "reduced hash",rh
        return [dataanalysis.DataHandle('processing_definition:'+rh[:8])]

class BasicEventProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=ISGRIEvents(assume=ScWData(input_scwid="any",use_abstract=True)) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
       # print "one scw hash:",ahash
        #ahash=dataanalysis.hashe_replace_object(ahash,'AnyScW','None')
        print "generalized hash:",ahash
        rh=dataanalysis.shhash(ahash)
        print "reduced hash",rh
        return [dataanalysis.DataHandle('processing_definition:'+rh[:8])]

class ImageProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=ii_skyimage(assume=ScWData(input_scwid="any",use_abstract=True)) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
        print "one scw hash:",ahash
        ahash=dataanalysis.hashe_replace_object(ahash,'AnyScW','None')
        print "generalized hash:",ahash
        rh=dataanalysis.shhash(ahash)
        print "reduced hash",rh
        d=dataanalysis.DataHandle('processing_definition:'+rh[:8])
        dataanalysis.AnalysisFactory.register_definition(d.handle,ahash)
        d.hash=ahash
        return [d]

class SpectraProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=ii_spectra_extract(assume=[ScWData(input_scwid="any",use_abstract=True),Revolution(input_revid="any",use_abstract=True)]) # arbitrary choice of scw, should be the same: assumption of course
        #mf=ii_spectra_extract(assume=[ScWData(input_scwid=AnyScW),Revolution(input_revid=AnyRevID)]) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
        #print "one scw hash:",ahash
        #ahash=dataanalysis.hashe_replacI#e_object(ahash,'AnyScW','None')
        #ahash=dataanalysis.hashe_replace_object(ahash,'AnyRevID','None')
        print "generalized hash:",ahash
        rh=dataanalysis.shhash(ahash)
        print "reduced hash",rh
        return [dataanalysis.DataHandle('processing_definition:'+rh[:8])]

class IDScWList(DataAnalysis):
    scwid_list=None
    allow_alias=True

    def main(self):
        self.scwlistdata=[ScWData(input_scwid=s.strip()) for s in self.scwid_list]


class RevScWList(DataAnalysis):
    input_rev=Revolution

    run_for_hashe=False
    allow_alias=True

    def main(self):
        import os

        event_files=glob.glob(self.input_rev.revroot+"/*/isgri_events.fits*")

        print "event files in",self.input_rev.revroot
        print "found event files",event_files

        scwids=sorted([fn.split("/")[-2] for fn in event_files])

        self.scwlistdata=[ScWData(input_scwid=s) for s in scwids if s[:-1].endswith("0010.00")]

    def __iter__(self):
        return iter(self.scwlistdata)

    def __repr__(self):
        return "[RevScWList:%s]"%repr(self.input_rev)

class ScWFilter(DataAnalysis):
    version="nofilter"
    input_list=RevScWList

    def main(self):
        return self.input_list

class BadListScWFilter(DataAnalysis):
    input="badlist"
    #input_lists

    def main(self):
        open(os.environ['REP_BASE_PROD_global']+"/")


class ScWListFiltered(DataAnalysis):
    input_filter=ScWFilter
    input_list=RevScWList
    #allow_alias=True

    def main(self):
        import os

        event_files=glob.glob(self.input_rev.revroot+"/*/isgri_events.fits*")

        scwids=[fn.split("/")[-2] for fn in event_files]

        self.scwlistdata=[ScWData(input_scwid=s) for s in scwids]

class PickFewScWList(DataAnalysis):
    nscw=10
    firstscws=True

    firstscw=0
    step=1

    input_list=RevScWList
    allow_alias=True

    version="v2.."

    def get_version(self):
        if self.firstscws:
            v=self.get_signature()+".firstscw.%i"%self.nscw+"."+self.version
        else:
            v=self.get_signature()+".lastscw.%i"%self.nscw+"."+self.version
        if self.firstscw>0:
            v=v+".from.%i"%self.firstscw
        if self.step>1:
            v=v+".step.%i"%self.step
        return v

    def main(self):
        thelist=sorted(self.input_list.scwlistdata,key=lambda x:x.input_scwid.handle)
        print "available list",len(thelist),thelist
        print "first scw:",self.firstscw
        if self.firstscws:
            self.scwlistdata=thelist[self.firstscw:self.nscw+self.firstscw]
        else:
            self.scwlistdata=thelist[-self.nscw:]
        if self.step!=1:
            self.scwlistdata=self.scwlistdata[::self.step]
        print "resulting list:",len(self.scwlistdata),self.scwlistdata

class FileScWList(DataAnalysis):
    input_fn=None
    allow_alias=True
    maxscw=None

    def main(self):
        if self.maxscw is not None:
            self.scwlistdata=[ScWData(input_scwid=s.strip()) for s in open(self.input_fn.handle.split(":",1)[0]).readlines()[:self.maxscw]] # omg!!
        else:
            self.scwlistdata=[ScWData(input_scwid=s.strip()) for s in open(self.input_fn.handle.split(":",1)[0])]


class ScWList(DataAnalysis):
    input_list=None
    allow_alias=True

    def main(self):
        self.scwlistdata=self.input_list.scwlistdata


def fromUTC(utc):                                                                                                                                                                                              
    r=subprocess.check_output(["converttime","UTC",utc,""])                                                                                                                                                    
    d={}                                                                                                                                                                                                       
    for l in r.split("\n"):                                                                                                                                                                                    
        t=re.search("Output Time\((.*?)\): (.*?)$",l,re.S)                                                                                                                                                     
        #print l,t                                                                                                                                                                                             
        if t:                                                                                                                                                                                                  
            g=t.groups()                                                                                                                                                                                       
            d[g[0]]=g[1]                                                                                                                                                                                       
            #print g                                                                                                                                                                                           
    return d    

