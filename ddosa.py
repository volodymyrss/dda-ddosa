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

from dataanalysis.core import DataFile
from dataanalysis.caches import cache_core
from dataanalysis import hashtools
from dataanalysis.hashtools import shhash
import dataanalysis.printhook
import dataanalysis.core as da

import pilton
from pilton import heatool 

import pprint
import os,shutil,re,time,glob
from os import access, R_OK
from os.path import isfile
from astropy.io import fits
from astropy import wcs
from astropy import wcs as pywcs
import subprocess,os
import ast
import copy
import re

try:
    import pandas as pd
except ImportError:
    print("no pandas!")

try:
    import yaml
except ImportError:
    print("no yaml!")


import numpy as np

if hasattr(da,'DataAnalysisPrototype'):
    DataAnalysisPrototype = da.DataAnalysisPrototype
else:
    DataAnalysisPrototype = da.DataAnalysis

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
        raise Exception("deprecated cache!")

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
        #print("will construct INTEGRAL cached file path for",hashe)

        scw=self.get_scw(hashe)
        rev=self.get_rev(hashe)

        def hash_to_path2(hashe):                                                                                                                                      
            return shhash(hashe[1])[:8]                                                                                                           
            
        marked=self.get_marked(hashe[1])
        marked=remove_repeating(marked)
        if dataanalysis.printhook.global_log_enabled: print("marked",marked)
        
        if not isinstance(scw,str):
            print("emergent scw:",scw)
            scw=None
        
        if scw=="Any":
            print("any scw:",scw,hashe)
            scw=None

        print("scw:",scw)
        print("rev:",rev)


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

        print(self,"cached path:",r)
                                                                                                                                                                       
        return r # choose to avoid overlapp    

class MemCacheIntegralBase: 
    def get_scw(self,hashe):                                                                                                                                       
        if isinstance(hashe,tuple):                                                                                                                                
            if hashe[0]=="analysis": # more universaly                                                                                                             
                if hashe[2].startswith('ScWData'):
                    if isinstance(hashe[1],tuple):
                        return hashe[1][-1]
                    else:
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
                if '...' in hashe[2]:
                    r+=hashe[2].split('...')[1:]
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
        #print("will construct INTEGRAL cached file path for",hashe)
        #hashe= hashtools.hashe_replace_object(hashe, None, "None")

        scw=self.get_scw(hashe)
        rev=self.get_rev(hashe)

        def hash_to_path2(hashe):                                                                                                                                      
            return shhash(hashe[1])[:8]                                                                                                           
            
        marked=self.get_marked(hashe[1])
        marked=remove_repeating(marked)
        if dataanalysis.printhook.global_log_enabled: print("marked",marked)

        for mark in marked:
            hashe=hashtools.hashe_replace_object(hashe,mark+"..","any")

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
                hashe=hashtools.hashe_replace_object(hashe,rev,"any")
                #print("reduced hashe",hashe)
                if dataanalysis.printhook.global_log_enabled: print("cached rev:",rev)
                r=self.filecacheroot+"/byrev/"+rev+"/"+hashe[2]+"/"+"/".join(marked)+"/"+hash_to_path2(hashe)+"/" # choose to avoid overlapp    
        else:
            hashe=hashtools.hashe_replace_object(hashe,scw,"any")
            hashe=hashtools.hashe_replace_object(hashe,('analysis', scw[:4], 'Revolution.v0'),('analysis', 'any', 'Revolution.v0'))

            #str("reduced hashe:",hashe,hash_to_path2(hashe))
            #if dataanalysis.printhook.global_log_enabled: print("reduced hashe:",hashe,hash_to_path2(hashe))
            #open("reduced_hashe.txt","w").write(hash_to_path2(hashe)+"\n\n"+pprint.pformat(hashe)+"\n")
            print(scw,hashe[2],marked)

            r=self.filecacheroot+"/byscw/"+scw[:4]+"/"+scw+"/"+hashe[2]+"/"+"/".join(marked)+"/"+hash_to_path2(hashe)+"/" # choose to avoid overlapp    

        #if dataanalysis.printhook.global_log_enabled: print("cached path:",r)

        print(self,"cached path:",r)
                                                                                                                                                                       
        return r # choose to avoid overlapp    



#class MemCacheIntegral(MemCacheIntegralBase,dataanalysis.MemCacheMySQL):
#    pass

#class MemCacheIntegralLegacy(MemCacheIntegralBase,dataanalysis.MemCacheSqlite):
#    pass


class MemCacheIntegralFallback(MemCacheIntegralBase,dataanalysis.caches.cache_core.CacheNoIndex):
    def store(self, hashe, obj):
        filepath=self.construct_cached_file_path(hashe,obj)

        return dataanalysis.caches.cache_core.CacheNoIndex.store(self,hashe,obj)

    def restore(self, hashe, obj, restore_config=None):
        filepath = self.construct_cached_file_path(hashe, obj)

        return dataanalysis.caches.cache_core.CacheNoIndex.restore(self, hashe, obj, restore_config)



        #class MemCacheIntegralFallbackOldPath(MemCacheIntegralBaseOldPath,dataanalysis.caches.core.CacheNoIndex):
    #readonly_cache=True

#class MemCacheIntegralIRODS(MemCacheIntegralBase,dataanalysis.MemCacheIRODS):
#    pass

#mc=dataanalysis.TransientCacheInstance
#mcg=MemCacheIntegral('/Integral/data/reduced/ddcache/')
#mc=mcg
#mc.parent=mcg
#mcgl=MemCacheIntegralLegacy('/Integral/data/reduced/ddcache/')
#mcg.parent=mcgl

IntegralCacheRoots=os.environ.get('INTEGRAL_DDCACHE_ROOT','/sps/integral/data/reduced/ddcache/')

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

    #mcgfb_oldp=MemCacheIntegralFallbackOldPath(IntegralCacheRoot)
    #mcgfb_oldp.readonly_cache=ro_flag
    #mcgfb.parent=mcgfb_oldp
    #CacheStack.append(mcgfb_oldp)

#mcgirods=MemCacheIntegralIRODS('/tempZone/home/integral/data/reduced/ddcache/')
#CacheStack[-1].parent=mcgirods
#CacheStack.append(mcgirods)

mc=CacheStack[0]

print("cache stack:",CacheStack)

class OSA_tool_kit_class(object):
    tool_versions=None

    def get_tool_version(self,name):
        if self.tool_versions is None:
            self.tool_versions={}

        if name not in self.tool_versions:
            cl=subprocess.Popen([name,"--version"],stdout=subprocess.PIPE)
            self.tool_versions[name]=cl.stdout.read().strip().split()[-1]

        return self.tool_versions[name]

OSA_tool_kit=OSA_tool_kit_class()

def get_OSA_tools(names=None):
    if names is None:
        return da.NoAnalysis

    if isinstance(names,str) or isinstance(names,tuple):
        names=[names]

    
    names=[ name if isinstance(name,tuple) else (name,None)
                for name in names ]
        

    class OSA_tools(DataAnalysisPrototype):
        osa_tools=names[:]
        
        def get_version(self):
            v=self.get_signature()+"."+self.version
            for osa_tool,default_version in self.osa_tools:
                current_version=OSA_tool_kit.get_tool_version(osa_tool)
                if default_version is None or current_version!=default_version:
                    v+="."+osa_tool+"_"+current_version
            return v

    return OSA_tools

class DataAnalysis(DataAnalysisPrototype):
    cache=mc

    write_caches=[cache_core.TransientCache,MemCacheIntegralFallback]
    read_caches=[cache_core.TransientCache,MemCacheIntegralFallback] #,MemCacheIntegralFallbackOldPath]

    input_osatools=get_OSA_tools()

    cached=False

    def get_scw(self):
        if self._da_locally_complete is not None:
            try:
                return "(completescw:%s)"%self.cache.get_scw(getattr(self,'_da_locally_complete',None))
            except:
                return "(complete)"

        for a in self.assumptions:
            if isinstance(a,ScWData):
                return "(assumescw:%s)"%str(a.input_scwid)

        return ""

    def __repr__(self):
        return "[%s%s%s%s%i]"%(self.get_version(),self.get_scw(),";Virtual" if self.virtual else "",";Complete" if self._da_locally_complete else "",id(self))


class NoScWData(da.AnalysisException):
    pass

class NoDeadData(da.AnalysisException):
    pass
    
class NoISGRIEvents(da.AnalysisException):
    pass

class EmptyScWList(da.AnalysisException):
    pass

class NoValidScW(da.AnalysisException):
    pass

class EmptyImageList(da.AnalysisException):
    pass

class ScWDataCorrupted(da.AnalysisException):
    pass

class FractionalEnergyBinsNotAllowed(da.AnalysisException):
    pass

class EfficiencyNotComputed(da.AnalysisException):
    pass

def good_file(fn):
    return os.path.exists(fn) and isfile(fn) and access(fn, R_OK)

class ScWData(DataAnalysis):
    input_scwid=None

    cached=False # how do we implment that this can change?

    schema_hidden=True

    version="v1"

    scwver="001"
    auxadpver="001"

    def main(self):
        try:
            self.scwid=self.input_scwid.handle
        except:
            self.scwid=self.input_scwid
        self.scwver=self.scwid[-3:]
        self.revid=self.scwid[:4]
        
        try:
            print("searching in "+detect_rbp(self.scwver))
            self.assume_rbp(detect_rbp(self.scwver))
        except da.AnalysisException:
            if self.scwver=="000":
                print("searching in "+detect_rbp(self.scwver)+"/nrt")
                self.assume_rbp(detect_rbp(self.scwver)+"/nrt")
            else:
                raise


    def test_scw(self):
        try:
            f=fits.open(self.scwpath+"/swg.fits")
            print("valid file:",f)
        except IOError as e:
            if e.message=="Header missing END card.":
                raise ScWDataCorrupted(self.scwpath,e.message)
            else:
                raise 

    def test_isgri_events(self):
        print("checking for readable events...")

        options=[self.scwpath+"/isgri_events.fits"]
        options.append(options[-1]+".gz")

        for option in options:
            if good_file(option):
                print("ok:",option)
                return

        raise NoISGRIEvents("no usable event data for: "+repr(self.scwid))


    def assume_rbp(self,rbp):
        self.revdirver = None
        tried = []
        for v in self.scwver, "002", "001", "000":
            p = rbp+"/scw/"+self.revid+"/rev."+v
            if os.path.exists(p):
                self.revdirver = v
                print("found revdir ver", v)
                break
            tried.append(p)

        if self.revdirver is None:
            raise Exception("no revdir available! tried: {}".format(tried))

        self.scwpath=rbp+"/scw/"+self.revid+"/"+self.scwid
        self.revdirpath=rbp+"/scw/"+self.revid+"/rev."+self.revdirver
        #self.auxadppath=rbp+"/aux/adp/"+self.revid+"."+self.auxadpver

        print("searching for auxadpver")
        for ver in "000", "001":
            self.auxadppath=rbp+"/aux/adp/"+self.revid+"."+ver
            self.auxadpver=ver
            if os.path.exists(self.auxadppath):
                print("auxadpver picked", self.auxadpver, self.auxadppath)
                break


        if not good_file(self.scwpath+"/swg.fits"):
            if not good_file(self.scwpath+"/swg.fits.gz"):
                print("failed searching for",self.scwpath+"/swg.fits")
                raise NoScWData("no scw data for: "+repr(self.scwid))
                #raise Exception("no scw data!")
            else:
                self.swgpath=self.scwpath+"/swg.fits.gz"
        else:
            self.swgpath=self.scwpath+"/swg.fits"
        print("swgpath:",self.swgpath)


    def get_isgri_events(self):
        if hasattr(self,'isgrievents'):
            return self.isgrievents.get_path()
        return self.scwpath+"/isgri_events.fits.gz"

    def get_telapse(self):
        return fits.open(self.swgpath)[1].header['TELAPSE']
    
    def get_t(self):
        h=fits.open(self.swgpath)[1].header
        return (h['TSTOP']+h['TSTART'])/2.,(h['TSTOP']-h['TSTART'])/2.
    
    def get_t1_t2(self):
        h=fits.open(self.swgpath)[1].header
        return h['TSTART'],h['TSTOP']

    def __repr__(self):
        return "[%s:%s]"%(self.__class__.__name__,self.input_scwid)

def detect_rbp(scwver="001"):
    if scwver=="001":
        if "REP_BASE_PROD_CONS" in os.environ:
            return os.environ["REP_BASE_PROD_CONS"]
    
    if scwver=="000":
        if "REP_BASE_PROD_NRT" in os.environ:
            return os.environ["REP_BASE_PROD_NRT"]

    return os.environ["REP_BASE_PROD"]

class ScWVersion(DataAnalysis):
    scwver="000"

    def get_version(self):
        return DataAnalysis.get_version(self)+".ver"+self.scwver

class Revolution(DataAnalysis):
    input_revid=None

    input_scwver=da.NoAnalysis

    scwver="001"
    auxadpver="001"

    def get_revid(self):
        return self.input_revid.handle

    def main(self):
        if not isinstance(self.input_scwver,da.NoAnalysis) and not self.input_scwver == da.NoAnalysis:
            self.scwver=self.input_scwver.scwver

        rbp=detect_rbp(scwver=self.scwver)

        self.revroot=rbp+"/scw/%s/"%self.get_revid()
        self.revdir=self.revroot+"/rev.%s/"%self.scwver
        self.auxadppath=rbp+"/aux/adp/"+self.get_revid()+"."+self.auxadpver

    def get_ijd(self):
        r1100=4306.5559396296
        r100=1315.4808007407

        r=int(self.get_revid())
        return r100+(r1100-r100)/1000*(r-100)

    def get_ijd_exact(self):
        ijd_b = map(float, converttime("REVNUM",self.get_revid(),"IJD"))
        return (ijd_b[1] + ijd_b[0]) / 2., (ijd_b[1] - ijd_b[0])
    
    def __repr__(self):
        return "[Revolution:%s:%s]"%(self.input_revid,self.scwver)

class RevForScW(DataAnalysis):
    input_scw=ScWData    
    run_for_hashe=True
    allow_alias=False

    def __repr__(self):
        return "[RevForScW:for %s]"%repr(self.input_scw)

    def main(self):
        revid=self.input_scw.input_scwid.handle[:4]
        scwver=self.input_scw.input_scwid.handle[-3:]
        print("revolution id for scw:",revid)
        return Revolution(input_revid=revid,use_scwver=scwver)
        #return Revolution(input_revid=revid,use_scwver=scwver)

class Rev4ScW(Revolution):
    input_scw=ScWData    
    input_revid=da.NoAnalysis

    def __repr__(self):
        return "[Rev4ScW:for %s]"%repr(self.input_scw)

    def get_revid(self):
        revid=self.input_scw.input_scwid.handle[:4]
        print("revolution id for scw:",revid)
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

class FailingMedia(DataAnalysis):
    def main(self):
        raise Exception("exampliary failure")



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
            vstart=fits.open(gmfile)[1].header['VSTART']
        
            print(vstart,ver)
            if vstart>newest_vstart and t>vstart:
                newest_vstart=vstart
                newest_ver=ver

        self.godol=self.input_ibisic.ibisicroot+"/cal/ibis_isgr_gain_offset_0010.fits"
        self.supgdol=self.input_ibisic.ibisicroot+"/mod/isgr_gain_mod_"+newest_ver+".fits[ISGR-GAIN-MOD,1,BINTABLE]"
        self.supodol=self.input_ibisic.ibisicroot+"/mod/isgr_off2_mod_"+newest_ver+".fits[ISGR-OFF2-MOD,1,BINTABLE]"
        self.risedol=self.input_lut2.datafile

class BinnedDataProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=BinEventsImage(assume=ScWData(input_scwid="any",use_abstract=True)) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
       # print("one scw hash:",ahash)
        #ahash=hashtools.hashe_replace_object(ahash,'AnyScW','None')
        print("generalized hash:",ahash)
        rh=shhash(ahash)
        print("reduced hash",rh)
        return [dataanalysis.DataHandle('processing_definition:'+rh[:8])]

class BasicEventProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=ISGRIEvents(assume=ScWData(input_scwid="any",use_abstract=True)) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
       # print("one scw hash:",ahash)
        #ahash=hashtools.hashe_replace_object(ahash,'AnyScW','None')
        print("generalized hash:",ahash)
        rh=shhash(ahash)
        print("reduced hash",rh)

        handle=dataanalysis.DataHandle('processing_definition:' + rh[:8])

        self.factory.note_factorization(dict(
            origin_object=self.__class__.__name__,
            origin_module=__name__,
            generalized_hash=ahash,
            reduced_hash=rh,
            handle=handle.handle,
        ))

        return [handle]

class ImageProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=ii_skyimage(assume=ScWData(input_scwid="any",use_abstract=True)) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
        print("one scw hash:",ahash)
        ahash=hashtools.hashe_replace_object(ahash,'AnyScW','None')
        print("generalized hash:",ahash)
        rh=shhash(ahash)
        print("reduced hash",rh)
        d=dataanalysis.DataHandle('processing_definition:'+rh[:8])
        dataanalysis.AnalysisFactory.register_definition(d.handle,ahash)

        self.factory.note_factorization(dict(
            origin_object=self.__class__.__name__,
            origin_module=__name__,
            generalized_hash=ahash,
            reduced_hash=rh,
            handle=d.handle,
        ))

        d.hash=ahash
        return [d]

class LCProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=ii_lc_extract(assume=ScWData(input_scwid="any",use_abstract=True)) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
        print("one scw hash:",ahash)
        ahash=hashtools.hashe_replace_object(ahash,'AnyScW','None')
        print("generalized hash:",ahash)
        rh=shhash(ahash)
        print("reduced hash",rh)

        handle = dataanalysis.DataHandle('processing_definition:' + rh[:8])
        self.factory.note_factorization(dict(
            origin_object=self.__class__.__name__,
            origin_module=__name__,
            generalized_hash=ahash,
            reduced_hash=rh,
            handle=handle,
        ))


        dataanalysis.AnalysisFactory.register_definition(handle.handle, ahash)
        handle.hash = ahash
        return [handle]

class SpectraProcessingSummary(DataAnalysis):
    run_for_hashe=True

    def main(self):
        mf=ii_spectra_extract(assume=[ScWData(input_scwid="any",use_abstract=True),Revolution(input_revid="any",use_abstract=True)]) # arbitrary choice of scw, should be the same: assumption of course
        #mf=ii_spectra_extract(assume=[ScWData(input_scwid=AnyScW),Revolution(input_revid=AnyRevID)]) # arbitrary choice of scw, should be the same: assumption of course
        ahash=mf.process(output_required=False,run_if_haveto=False)[0]
        #print("one scw hash:",ahash)
        #ahash=dataanalysis.hashe_replacI#e_object(ahash,'AnyScW','None')
        #ahash=hashtools.hashe_replace_object(ahash,'AnyRevID','None')
        print("generalized hash:",ahash)
        rh=shhash(ahash)
        print("reduced hash",rh)
        handle=dataanalysis.DataHandle('processing_definition:'+rh[:8])
        self.factory.note_factorization(dict(
            origin_object=self.__class__.__name__,
            origin_module=__name__,
            generalized_hash=ahash,
            reduced_hash=rh,
            handle=handle.handle,
        ))
        return [handle]

class ibis_isgr_energy_standard(DataAnalysis):
    cached=False
    
    input_scw=ScWData()
    #input_raw_events=None
    #input_ibis_hk=None
    # can also take only isgri events

    input_ecorrdata=GetEcorrCalDB

    version="v4_extras"

    osa_tools=["ibis_isgr_energy"]
   
    def main(self):
        self.input_scw.test_scw()
        self.input_scw.test_isgri_events()

        remove_withtemplate("isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)")
    
        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/isgri_events.fits[ISGR-EVTS-ALL]", \
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
        self.input_scw.test_isgri_events()

        remove_withtemplate("isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)")

        if not os.path.exists(self.input_scw.scwpath+"/isgri_events.fits") and not os.path.exists(self.input_scw.scwpath+"/isgri_events.fits.gz"):
            raise NoISGRIEvents()

        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/isgri_events.fits[ISGR-EVTS-ALL]", \
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
            self.input_scw.scwpath+"/isgri_events.fits[ISGR-EVTS-ALL]", \
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
    input_scw=ScWData
    input_ic=ICRoot
    #input_gticreate=gti_create()

    cached=True
    
    version="v2"
    def main(self):
        # horrible horrible full OSA

        self.input_scw.test_scw()
        self.input_scw.test_isgri_events()

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
                print("full path",path)
            else:
                print("not a full path",path)
                path="../../"+path
            ht['GTI_gtiUserI']=path
            ht['GTI_TimeFormat']='UTC'

        ht.run()

        shutil.copy(ht.cwd+"/ibis_gti.fits","./ibis_gti.fits")
        self.output_gti=DataFile("ibis_gti.fits")

        gti=fits.open("ibis_gti.fits")[-1].data
        print(gti)
        

# maybe split indeed,but try to show another way
class ibis_dead(DataAnalysis):
    input_scw=ScWData
    input_ic=ICRoot

    cached=True
    
    version="v2"
    def main(self):
        # horrible horrible full OSA
        self.input_scw.test_scw()
        self.input_scw.test_isgri_events()

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

        if not os.path.exists(ht.cwd+"/isgri_dead.fits"):
            print("not found dead!")
            raise NoDeadData()

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

    autoversion=False

    rmfbins=False

    def get_version(self):
        v=self.get_signature()+"."+self.version

        if self.autoversion:
            if self.ebins is None:
                v+=".std_one_25_80"
            else:
                if len(self.ebins)==1:
                    v+=".one_bin_%.15lg_%.15lg"%(self.ebins[0][0],self.ebins[0][1])
                else:
                    v+=".%i_bins"%len(self.ebins)
                    for ebin in self.ebins:
                        v+=".%.15lg_%.15lg"%(ebin[0],ebin[1])

        return v


    def main(self):
        if self.ebins is None:
            self.bins=[(25,80)]
        else:
            self.bins=self.ebins

        for e1,e2 in self.bins:
            if abs(round(e1*2)/2. - e1)>1e-5 or abs(round(e2*2)/2. - e2)>1e-5:
                raise FractionalEnergyBinsNotAllowed()

class SpectraBins(DataAnalysis):
    input_binsname="spectral_bins_62"

    rmfpath=None
    
    rmfbins=True
    rmfext=3

    version="v3"
    def main(self):
        if self.rmfpath is not None:
            self.binrmf=self.rmfpath
        else:
            self.binrmf=os.environ['INTEGRAL_DATA']+"/resources/rmf_62bands.fits"
        
        if not os.path.exists(self.binrmf):
            self.binrmf=os.environ.get("INTEGRAL_RESOURCES","/data/resources")+"/rmf_62bands.fits"

        if not os.path.exists(self.binrmf):
            self.binrmf="/unsaved_data/savchenk/rmf_62bands.fits"

        #self.binrmf=os.environ['CURRENT_IC']+"/ic/ibis/rsp/rmf_62bands.fits" # noo!!!
        #self.binrmf=os.environ['CURRENT_IC']+"/ic/ibis/rsp/isgr_ebds_mod_0001.fits" # noo!!!
        e=fits.open(self.binrmf)[self.rmfext].data
        self.bins=zip(e['E_MIN'],e['E_MAX'])
        self.binrmfext=self.binrmf+'[%i]'%self.rmfext

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
            self.input_scw.scwpath+"/isgri_events.fits[ISGR-EVTS-ALL]", \
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

        print("target_level",self.target_level)

        print("has rmfbins "+str(self.input_bins.rmfbins) if hasattr(self.input_bins,'rmfbins') else "no rmfbins")
        print("has binrmfext" if hasattr(self.input_bins,'binrmfext') else "no binrmfext")

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
        ht['isgri_t_len'] = 10000000 if not hasattr(self,'input_timebin') else self.input_timebin.time_bin_seconds
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

class LCEnergyBins(ImageBins):
    pass

class LCTimeBin(DataAnalysis):
    time_bin_seconds=100

    def get_version(self):
        return self.get_signature()+"."+self.version+".t%.5lg"%self.time_bin_seconds



class BinEventsLC(BinEventsVirtual):
    target_level="BIN_T"
    input_timebin=LCTimeBin
    input_bins = LCEnergyBins


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
            print("will use uniformity:",maps['unif'])
        
        if hasattr(self,'input_bkg'):
            maps['back']=('Bkg',self.input_bkg.bkg.get_path()+"[1]")
            print("will use background:",maps['back'])

        level2key={
                'BIN_I':'ima',
                'BIN_S':'spe',
                'BIN_T':'lcr',
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

class BinMapsLC(BinMapsVirtual):
    target_level="BIN_T"
    input_bins=LCEnergyBins

class GRcat(DataAnalysis):
    input_ic=ICRoot
    
    suffix=None

    cached=False # again, this is transient-level cache

    userefcatvar=False
    useresources=True

    refcat_version="42"

    def get_version(self):
        v=self.get_signature()+"."+self.version

        if self.useresources:
            self.cat=os.environ.get("INTEGRAL_RESOURCES","/data/resources")+"/gnrl_refr_cat_00%s.fits[1]"%self.refcat_version
            self.catname=self.cat.split("/")[-1]
            v+=".resources_"+self.catname
        else:
            if self.userefcatvar:
                self.cat=os.environ["ISDC_REF_CAT"]
                self.catname=self.cat.split("/")[-1]
                v+="var_"+self.catname
            else:
                if self.suffix is not None:
                    v=v+"."+self.suffix
        return v

    def main(self):
        if self.useresources:
            pass
        elif self.userefcatvar:
            pass
        else:
            if self.suffix is None:
                self.cat=detect_rbp()+"/cat/hec/gnrl_refr_cat_0043.fits[1]"
            else:
                self.cat=detect_rbp()+"/cat/hec/gnrl_refr_cat_0043_%s.fits[1]"%self.suffix


class BrightCat(DataAnalysis):
    input=GRcat
    #input_selection="flag5"
    input_selection="flag5andover100"

    cached=False

    def main(self):
        #self.cat=self.input.cat+"[ISGRI_FLAG2==5]"
        #self.cat=self.input.cat+"[ISGRI_FLAG2==5]"
        self.cat_path=self.input.cat+"[ISGRI_FLAG2==5&&ISGR_FLUX_1>100]"

        fn="very_bright_cat.fits"

        ht=heatool("fextract")
        ht['infile']=self.cat_path
        ht['outfile']=fn
        ht['clobber']='yes'
        ht.run()

        self.cat=da.DataFile(fn)


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
        
        att=self.input_scw.auxadppath+"/attitude_historic.fits"
        if os.path.exists(att):
            att=self.input_scw.auxadppath+"/attitude_historic.fits[AUXL-ATTI-HIS,1,BINTABLE]"
            attp=att
        else:
            att=self.input_scw.auxadppath+"/attitude_snapshot.fits[AUXL-ATTI-SNA,1,BINTABLE]"
            attp_g=glob.glob(self.input_scw.auxadppath+"/attitude_predicted_*.fits*")

            print("possible predicted attitude", attp_g)
            print("possible predicted attitude", self.input_scw.auxadppath+"/attitude_predicted_*.fits*")

            attp=attp_g[0]+"[AUXL-ATTI-PRE,1,BINTABLE]"


        construct_gnrl_scwg_grp(self.input_scw,[\
                    catfn,
                    self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]",
                    self.input_gti.output_gti.path,
#                    att,
#                    attp,
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
        print(self.input_shadows)


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

        try:
            ht.run()
        except Exception as e:
            print("problem:", repr(e))
            print("tool", ht.output)

            if 'Verif8Bins           Status :          0 Bands error 2 in TIME bins' in ht.output:
                raise EfficiencyNotComputed(ht.output)

            raise

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

class ShadowUBCLC(ShadowUBCVirtual):
    level="BIN_T"
    input_shadows=BinEventsLC
    input_maps=BinMapsLC
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
    
    cached=False

    gb_binary=None

    def main(self):
        self.input_scw.test_isgri_events()

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

    cached=True

class ghost_bustersSpectra(ghost_bustersVirtual):
    input_shadow=ShadowUBCSpectra
    level="BIN_S"


class ghost_bustersLC(ghost_bustersVirtual):
    input_shadow=ShadowUBCLC
    level="BIN_T"

class BinnedBackgroundSpectrumFromGB(DataAnalysis):
    input_gb=ghost_bustersSpectra

    copy_cached_input=False

    cached=True

    def main(self):
        corshad=fits.open(self.input_gb.corshad.get_path())

        corr_effi=[_e for _e in corshad[2:] if _e.header['SHD_TYPE']=='EFFICIENCY']
        corr_dete=[_e for _e in corshad[2:] if _e.header['SHD_TYPE']=='DETECTOR']
        corr_var=[_e for _e in corshad[2:] if _e.header['SHD_TYPE']=='VARIANCE']

        spectrum=[]
        for ef,dete,var in zip(corr_effi,corr_dete,corr_var):
            print("reading:",ef.header['E_MIN'])
            spectrum.append(dict(
                            e1=ef.header['E_MIN'],
                            e2=ef.header['E_MAX'],
                            effi_mean=ef.data.mean(),
                            effi_sum=np.nansum(ef.data),
                            counts_mean=dete.data.mean(),
                            counts_sum=np.nansum(dete.data),
                            var_mean=var.data.mean(),
                            var_sum=np.nansum(var.data),
                        ))

        spectrum=pd.DataFrame(spectrum)
        spectrum['ec']=(spectrum.e1+spectrum.e2)/2.

        fn="background_spectrum.csv"
        spectrum.to_csv(fn)
        self.spectrum=da.DataFile(fn)


class ISGRIRefCat(DataAnalysis):
    input=GRcat
    input_selection="onlyisgri33"

    cached=False # since the path is transient

    def main(self):
        self.cat=self.input.cat+"[ISGRI_FLAG==1 || ISGRI_FLAG==2]"

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


class SourceList(DataAnalysis):
    sources=[]

    def get_version(self):
        v=self.get_signature()+"."+self.version
        for source in self.sources:
            v+="%(name)s_%(ra).5lg_%(dec).5lg"%source
        return v

class CatExtract(DataAnalysis):
    input_cat=ISGRIRefCat
    input_scw=ScWData

    #input_extra_sources=SourceList
    
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

        fn="isgri_catalog.fits"
        if hasattr(self,'input_extra_sources'):
            f=fits.open("isgri_catalog.fits")
            t_orig=f[1]

            t_new=fits.BinTableHDU.from_columns(t_orig.columns,nrows=len(t_orig.data)+len(self.input_extra_sources.sources))
            t_new.data[:len(t_orig.data)]=t_orig.data[:]

            i_offset=len(t_orig.data)

            for i,source in enumerate(self.input_extra_sources.sources):
                print("adding",source)
                t_new.data[i_offset + i]['NAME'] = source['name']
                t_new.data[i_offset + i]['RA_OBJ'] = source['ra']
                t_new.data[i_offset + i]['DEC_OBJ'] = source['dec']
                t_new.data[i_offset + i]['ISGRI_FLAG'] = 1

            f[1].data=t_new.data

            fn="isgri_catalog_extra.fits"
            f.writeto(fn,clobber=True)

        print("storing cat as",fn)
        self.cat=DataFile(fn)


class ImagingConfig(DataAnalysis):
    input="onesource_negmod"

    SearchMode=3
    ToSearch=5
    CleanMode=1
    MinCatSouSnr=4
    MinNewSouSnr=5
    NegModels=1
    DoPart2=1
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
        
    empty_results=False

    def get_version(self):
        v=self.get_signature()+"."+self.version
        for k in ['SouFit','SearchMode','ToSearch','CleanMode','MinCatSouSnr','MinNewSouSnr','NegModels','DoPart2']: # dopart2 is flow control, separately
            if hasattr(self,'ii_'+k):
                v+="_"+k+"_"+str(getattr(self,'ii_'+k))
        return v

    def treat_input_analysis_exceptions(self,exceptions):
        for e in exceptions:
            print("ii_skyimage experienced",e)

            if isinstance(e[1],NoDeadData):
                self.empty_results=True
                continue

            if isinstance(e[1],NoISGRIEvents):
                self.empty_results=True
                continue

            return False

        return True


    def main(self):
        print("results marker",self.empty_results)
        print("")
        if self.empty_results:
            print("it's so empty... probably excepted")
            return

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
            print("no image produced: since there was no exception in the binary, assuming empty results")
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
                print("something went wrong")
    
        self.post_process()

    def post_process(self):
        pass

class ImageGroups(DataAnalysis):
    input_scwlist=None
    input_image_processing=ImageProcessingSummary

    allow_alias=True
    run_for_hashe=False

    copy_cached_input=False
    #copy_cached_input=True

    outtype="BIN_I"

    version="v1.1"

    def construct_og(self,og_fn):
        scw_og_fns = []

        if len(self.members)==0:
            raise EmptyImageList()

        total_extracted_cat=None
        total_skyres=None
        total_srclres=None

        for scw,image,gb,gti,cat in self.members:
            fn = "og_%s.fits" % scw.input_scwid.str()
            if not hasattr(image,'skyima'):
                print("skipping",scw)
                continue
            construct_gnrl_scwg_grp(scw, children=
                [
                    image.skyima.get_path(),
                    image.skyres.get_path(),
                    gti.output_gti.get_path(),
                    gb.corshad.get_path(),
     #               cat.cat._da_unique_local_path,
                    scw.auxadppath + "/time_correlation.fits[AUXL-TCOR-HIS]",
                ], fn=fn)

            import_attr(scw.scwpath + "/swg.fits",
                        ["OBTSTART", "OBTEND", "TSTART", "TSTOP", "SW_TYPE", "TELAPSE", "SWID", "SWBOUND"],fn)
            set_attr({'ISDCLEVL': self.outtype}, fn)
            set_attr({'INSTRUME': "IBIS"}, fn)

            scw_og_fns.append(fn)
    
            fe=fits.open(cat.cat.get_path())[1]
            #fe=fits.open(cat.cat._da_unique_local_path)[1]
            print(scw,"scw extracted cat has",len(fe.data),"sources",fe.data['NAME'][:10],"...")
            if total_extracted_cat is None:
                total_extracted_cat=fe
            else:
                total_extracted_cat.data=np.concatenate((total_extracted_cat.data,fe.data))
                print("total extracted cat has",len(total_extracted_cat.data),"sources before filtering")

                u,ui=np.unique(total_extracted_cat.data['NAME'],return_index=True)
                total_extracted_cat.data=total_extracted_cat.data[ui]
                print("total extracted cat has",len(total_extracted_cat.data),"sources after unique filtering")
            
            sfe=fits.open(image.skyres.get_path())[2] # one band
            if total_skyres is None:
                total_skyres=sfe
            else:
                total_skyres.data=np.concatenate((total_skyres.data,sfe.data))
            
            srfe=fits.open(image.srclres.get_path())[1] # one band
            if total_srclres is None:
                total_srclres=srfe
            else:
                total_srclres.data=np.concatenate((total_srclres.data,srfe.data))



        if total_extracted_cat is not None:
            total_extracted_cat.writeto("total_extracted_cat.fits",overwrite=True)
            self.total_extracted_cat=da.DataFile("total_extracted_cat.fits")

        if total_skyres is not None:
            total_skyres.writeto("total_skyres.fits",overwrite=True)
            self.total_skyres=da.DataFile("total_skyres.fits")
        
        if total_srclres is not None:
            total_srclres.writeto("total_srclres.fits",overwrite=True)
            self.total_srclres=da.DataFile("total_srclres.fits")
        
        construct_gnrl_scwg_grp_idx(scw_og_fns,fn="og_idx.fits")
        set_attr({'ISDCLEVL': self.outtype}, "og_idx.fits")


        construct_og(["og_idx.fits","total_skyres.fits","total_extracted_cat.fits"], fn=og_fn)
        #construct_og(["og_idx.fits","total_extracted_cat.fits"], fn=og_fn)

        set_attr({'ISDCLEVL': self.outtype}, og_fn)

    def main(self):
        self.members=[
            (
                scw,
                ii_skyimage(assume=[scw]),
                ghost_bustersImage(assume=[scw]),
                ibis_gti(assume=[scw]),
                CatExtract(assume=[scw])
            ) for scw in self.input_scwlist.scwlistdata
        ]

        if len(self.members)==0:
            raise EmptyScWList()

class LCGroups(DataAnalysis):
    input_scwlist=None
    input_lc_processing=LCProcessingSummary


    allow_alias=True
    run_for_hashe=True

    outtype="BIN_I"

    def construct_og(self,og_fn):
        scw_og_fns = []

        for scw,lc in self.members:
            if not hasattr(lc,'lightcurve'):
                print("lcgroups skipping this:",scw,lc)
                continue

            fn = "og_%s.fits" % scw.input_scwid.str()
            construct_gnrl_scwg_grp(scw, children=
                [
                    lc.lightcurve.get_path(),
                    scw.auxadppath + "/time_correlation.fits[AUXL-TCOR-HIS]",
                ], fn=fn)

            import_attr(scw.scwpath + "/swg.fits",
                        ["OBTSTART", "OBTEND", "TSTART", "TSTOP", "SW_TYPE", "TELAPSE", "SWID", "SWBOUND"],fn)
            set_attr({'ISDCLEVL': self.outtype}, fn)
            set_attr({'INSTRUME': "IBIS"}, fn)

            scw_og_fns.append(fn)

        construct_gnrl_scwg_grp_idx(scw_og_fns,fn="og_idx.fits")
        set_attr({'ISDCLEVL': self.outtype}, "og_idx.fits")

        construct_og(["og_idx.fits"], fn=og_fn)

        set_attr({'ISDCLEVL': self.outtype}, og_fn)

    def main(self):
        self.members=[
            (
                scw,
                ii_lc_extract(assume=[scw]),
            ) for scw in self.input_scwlist.scwlistdata
        ]

        if len(self.members)==0:
            raise EmptyScWList()

class lc_pick(DataAnalysis):
    input_lcgroups = LCGroups
    source_names=["Crab"]

    cached=True
    instrument="isgri"

    def get_version(self):
        try:
            return super(lc_pick, self).get_version()+"."+(".".join([m.replace(" ","_") for m in self.source_names]))
        except:
            return "lc_pick.UNDEFINED"

    def main(self):
        self.input_lcgroups.construct_og("ogg.fits")

        assert len(self.source_names)==1

        for source_name in self.source_names:
            fn = "lc_%s.fits" % source_name.replace(" ","_")
            remove_withtemplate(fn+"(ISGR-SRC.-LCR-IDX.tpl)")

            ht = heatool("lc_pick")
            ht['group'] = "ogg.fits[1]"
            ht['source']=source_name
            ht['instrument']=self.instrument
            ht['lc']=fn
            ht.run()
            setattr(self,'lightcurve',da.DataFile(fn))


    ## do this only isgri
#            d = fits.open(fn)[1]
#            t_lc = d.data['TIME']

#            timedel = data.header['TIMEDEL']
#            dt_lc = (timedel / 2) * np.ones(t_lc.shape)

#            for i in range(len(t_lc) - 1):
#                dt_lc[i + 1] = min(timedel / 2, t_lc[i + 1] - t_lc[i] - dt_lc[i])

#            d.data['XAX_E'] = dt_lc
#            d.writeto(fn, overwrite=True)
    ## do jemx dt_lc negative




class MosaicImagingConfig(DataAnalysis):
    input="500s2x88"

    SearchMode=2
    ToSearch=500
    CleanMode=1
    MinCatSouSnr=6
    MinNewSouSnr=6
    NegModels=0
    DoPart2=2
    SouFit=0



class mosaic_ii_skyimage(DataAnalysis):
    input_maps = BinMapsImage
    input_bins = ImageBins
    #input_cat = CatExtract
    input_ic = IBIS_ICRoot
    input_imgconfig = MosaicImagingConfig

    input_imagegroups=ImageGroups

    #input_gb = ghost_bustersImage
    #input_gti = ibis_gti

    cached = True
    copy_cached_input=False
    #copy_cached_input=True

    ii_skyimage_binary = None

    save_image = True

    image_tag = None

    version = "v2.2.4"

    outtype = "BIN_I"

    def get_version(self):
        v = self.get_signature() + "." + self.version
        for k in ['SouFit', 'SearchMode', 'ToSearch', 'CleanMode', 'MinCatSouSnr', 'MinNewSouSnr', 'NegModels']:
            if hasattr(self, 'ii_' + k):
                v += "_" + k + "_" + str(getattr(self, 'ii_' + k))
        return v

    def merge_cat(self):
        f_cat = fits.open(self.input_imagegroups.total_extracted_cat.get_path())
        f_sr = fits.open(self.input_imagegroups.total_skyres.get_path())

        m = f_sr[1].data['NEW_SOURCE'] != 0
        new_sources = f_sr[1].data[m]
        new_sources_srcl_cat = np.zeros(len(new_sources), f_cat[1].data.dtype)

        for k,v in [('RA_OBJ','RA_FIN'),('DEC_OBJ','DEC_FIN')]:
            new_sources_srcl_cat[k] = new_sources[v or k]

        f_cat[1].data = np.concatenate((f_cat[1].data, new_sources_srcl_cat))

        merged_cat_fn = "merged_srcl_cat.fits"
        f_cat.writeto(merged_cat_fn,clobber=True)

        return merged_cat_fn

    def main(self):

        #merged_cat_fn = self.merge_cat()
        #self.merged_cat=DataFile(merged_cat_fn)
        
        if self.ii_skyimage_binary is None:
            ii_skyimage_binary = "ii_skyimage"
        else:
            ii_skyimage_binary = self.ii_skyimage_binary


        def reset():
            self.input_imagegroups.construct_og("ogg.fits")
        
            self.total_extracted_cat=self.input_imagegroups.total_extracted_cat
            self.total_skyres=self.input_imagegroups.total_skyres

            remove_withtemplate("isgri_srcl_res.fits(ISGR-SRCL-RES.tpl)")
            remove_withtemplate("isgri_mosa_ima.fits(ISGR-MOSA-IMA-IDX.tpl)")
            remove_withtemplate("isgri_mosa_res.fits(ISGR-MOSA-RES-IDX.tpl)")
            remove_withtemplate("isgri_sky_ima.fits(ISGR-SKY-IMA-IDX.tpl)")
            remove_withtemplate("isgri_sky_res.fits(ISGR-SKY-RES-IDX.tpl)")


        warnings=[]


        reset()

        ht = heatool(ii_skyimage_binary)
        ht['outOG'] = "ogg.fits[1]"
        ht['outCat'] = "isgri_srcl_res.fits(ISGR-SRCL-RES.tpl)"
        ht['mask'] = self.input_ic.ibisicroot + "/mod/isgr_mask_mod_0003.fits[ISGR-MASK-MOD,1,IMAGE]"
        ht['deco'] = self.input_ic.ibisicroot + "/mod/isgr_deco_mod_0008.fits[ISGR-DECO-MOD,1,IMAGE]"
        ht['tungAtt'] = self.input_ic.ibisicroot + "/mod/isgr_attn_mod_0010.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['aluAtt'] = self.input_ic.ibisicroot + "/mod/isgr_attn_mod_0011.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['leadAtt'] = self.input_ic.ibisicroot + "/mod/isgr_attn_mod_0012.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['covrMod'] = self.input_ic.ibisicroot + "/mod/isgr_covr_mod_0002.fits[1]"
        ht['outMosIma'] = "isgri_mosa_ima.fits(ISGR-MOSA-IMA-IDX.tpl)"
        ht['outMosRes'] = "isgri_mosa_res.fits(ISGR-MOSA-RES-IDX.tpl)"
        ht['ScwDir'] = './'
        ht['ScwType'] = 'pointing'
        ht['ExtenType'] = 2
        ht['FastOpen'] = 1
        #ht['inCat'] = merged_cat_fn
        ht['MapSize'] = 80
        ht['OutType'] = self.outtype
        ht['num_band'] = len(self.input_bins.bins)
        ht['E_band_min'] = " ".join([str(a[0]) for a in self.input_bins.bins])
        ht['E_band_max'] = " ".join([str(a[1]) for a in self.input_bins.bins])
       # ht['DoPart2'] = 1
        for k in ['SouFit', 'SearchMode', 'ToSearch', 'CleanMode', 'MinCatSouSnr', 'MinNewSouSnr', 'NegModels', 'DoPart2']:
            ht[k] = getattr(self.input_imgconfig, k)
            if hasattr(self, 'ii_' + k): ht[k] = getattr(self, 'ii_' + k)
        ht['corrDol'] = self.input_maps.corr.path

        env = copy.deepcopy(os.environ)
            
        common_log_file = "common-log-file.txt"
        env['COMMONLOGFILE'] = "+"+os.getcwd()+"/"+common_log_file


        try:
            ht.run(env=env)
        except pilton.HEAToolException as e:
            if 'in MAIN__ at ii_skyimage_main.f90:94' in ht.output:
                new_catthr = int(float(ht['MinCatSouSnr'].value)*1.5)
                new_newthr = int(float(ht['MinNewSouSnr'].value)*1.5)
                warnings.append("""detected likely many-source segfault,
                                   increasing catalogue source significance threshold from %.5lg to %.5lg,
                                   increasing new source significance threshold from %.5lg to %.5lg"""%(
                                                ht['MinCatSouSnr'].value,
                                                new_catthr,
                                                ht['MinNewSouSnr'].value,
                                                new_newthr,
                                                ))
                print("WARNING:",warnings[-1])
                ht['MinCatSouSnr']=new_catthr
                ht['MinNewSouSnr']=new_newthr

                reset()

                ht.run(env=env)
            

        self.commonlog = DataFile(common_log_file)


        if not os.path.exists("isgri_mosa_ima.fits"):
            warnings.apped("no image produced: since there was no exception in the binary, assuming empty results")
            print("warnings", warnings[-1])
            self.empty_results = True
            return

        self.raw_srclres = DataFile("isgri_srcl_res.fits")

        if self.save_image:
            self.skyima = DataFile("isgri_mosa_ima.fits")

        self.skyres = DataFile("isgri_mosa_res.fits")

        self.srclres = self.merge_res()

        self.mosaic=self.skyima

        self.post_process()

        if len(warnings)>0:
            self.comment = "WARNING:" + ("; ".join(warnings))


    def merge_res(self):
        f_cat = fits.open(self.raw_srclres.get_path())
        f_sr = fits.open(self.input_imagegroups.total_srclres.get_path())

        new_sources_indices=[]

        for i in range(len(f_sr[1].data)):
            if f_sr[1].data[i]['NEW_SOURCE'] == 0:
                continue

            ra = f_sr[1].data[i]['RA_FIN']
            dec = f_sr[1].data[i]['DEC_FIN']


            if len(f_cat[1].data)>0 and min(abs(ra - f_cat[1].data['RA_OBJ']) + abs(f_cat[1].data['DEC_OBJ'])) < 15./60.:
                print(i,ra,dec,"already in the cat")
                continue

            if len(new_sources_indices)>0 and min([ (abs(ra - f_sr[1].data[o_i]['RA_FIN']) + abs(dec - f_sr[1].data[o_i]['DEC_FIN'])) for o_i in new_sources_indices]) < 15./60.:
                print(i,ra,dec,"already in the new")
                continue
                    

            new_sources_indices.append(i)


        new_sources_srcl_cat = np.zeros(len(new_sources_indices), f_cat[1].data.dtype)

        for j, i in enumerate(new_sources_indices):
            for k in f_cat[1].data.columns:
                new_sources_srcl_cat[j][k.name] = f_sr[1].data[i][k.name]

        print("new_sources_srcl_cat:",new_sources_srcl_cat)


        f_cat[1].data = np.concatenate((f_cat[1].data, new_sources_srcl_cat))

        merged_res_fn = "merged_res.fits"
        f_cat.writeto(merged_res_fn,clobber=True)

        return DataFile(merged_res_fn)

    def post_process(self):
        pass


class CatForSpectraFromImaging(DataAnalysis):
    input_imaging=ii_skyimage

    minsig=None
    maxsources=None

    def get_version(self):
        if isinstance(self.minsig,str):
            raise Exception("what is "+repr(self.minsig))
        v=self.get_signature()+"."+self.version+("" if self.minsig is None else "minsig%.3lg"%self.minsig)

        if hasattr(self,'input_extra_sources'):
            v+="_extrasv2"

        return v

    def main(self):
        if hasattr(self.input_imaging,'empty_results') and self.input_imaging.empty_results:
            print("no results here")
            self.empty_results=True
            return

        catfn="cat4spectra.fits"

        f=fits.open(self.input_imaging.srclres.path)

        print("image catalog contains",len(f[1].data))
        print("image catalog contains sig from",f[1].data['DETSIG'].min(),"to",f[1].data['DETSIG'].max())

        if self.minsig is not None:
            f[1].data=f[1].data[f[1].data['DETSIG']>self.minsig]
            print("selecting by significance",len(f[1].data))
        
        if self.maxsources is not None:
            raise Exception("not implemented")

        if hasattr(self, 'input_extra_sources'):
            t_orig = f[1]

            new_extra_sources=[]
            for source in self.input_extra_sources.sources:
                if source['name'] not in [str(s).strip() for s in t_orig.data['NAME']]:
                    new_extra_sources.append(source)

            print("new extra sources",new_extra_sources)

            t_new = fits.BinTableHDU.from_columns(t_orig.columns,
                                                    nrows=len(t_orig.data) + len(new_extra_sources))
            t_new.data[:len(t_orig.data)] = t_orig.data[:]

            i_offset = len(t_orig.data)

            for i, source in enumerate(new_extra_sources):
                print("adding", source)
                t_new.data[i_offset + i]['NAME'] = source['name']
                t_new.data[i_offset + i]['RA_OBJ'] = source['ra']
                t_new.data[i_offset + i]['DEC_OBJ'] = source['dec']
                t_new.data[i_offset + i]['ISGRI_FLAG'] = 1

            f[1].data = t_new.data

            catfn = "isgri_catalog_extra.fits"
            f.writeto(catfn, clobber=True)

        f.writeto(catfn,clobber=True)

        self.cat=DataFile(catfn)

#class CatForSpectra(DataAnalysisPrototype):
#    pass

class ISGRIResponse(DataAnalysis):
    input_ecorrdata=GetEcorrCalDB

    path=os.environ.get('ISGRI_RESPONSE',os.environ.get('INTEGRAL_DATA','')+"/resources/rmf_62bands.fits")

    
    def rmf_path(self):
        return self.path
        

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
        if hasattr(self.input_cat,'empty_results') and self.input_cat.empty_results:
            print("empty here")
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
        ht['idx_isgrResp']=self.input_response.rmf_path
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
            rmf=self.input_bins.get_binrmfext() if hasattr(self,'input_bins') else self.input_response.rmf_path
            ht['num_band'] = -1
            ht['idx_isgrResp'] = rmf #  +"[1]" will this work?


        #for k in ['SearchMode','ToSearch','CleanMode','MinCatSouSnr','MinNewSouSnr','NegModels','DoPart2']: # dopart2 is flow control, separately
        #    ht[k]=getattr(self.input_imgconfig,k)

        try:
            ht.run()
            self.spectrum=DataFile(spec_fn)
            self.pifs=DataFile(pif_fn)
        except pilton.HEAToolException as e:
            self.empty_results=True


class CatForLC(CatForSpectraFromImaging):
    pass


class ii_lc_extract(DataAnalysis):
    input_gb = ghost_bustersLC
    input_cat = CatForLC
    input_ic = IBIS_ICRoot
    input_scw = ScWData
    input_maps = BinMapsLC

    input_gti = ibis_gti

    cached = True

    version = "v1"

 #   input_bins=LCEnergyBins
    # input_cat=CatExtract
    # input_imgconfig=ImagingConfig

    shdtype = "BIN_T"
    binary = "ii_lc_extract"

    usebkg = True

    #report_runtime_destination = "mysql://pixels.runtime"

    def main(self):
        if hasattr(self.input_cat, 'empty_results') and self.input_cat.empty_results:
            print("empty here")
            self.empty_results = True
            return

        att = self.input_scw.auxadppath + "/attitude_historic.fits"
        if os.path.exists(att):
            att = self.input_scw.auxadppath + "/attitude_historic.fits[AUXL-ATTI-HIS,1,BINTABLE]"
            attp = att
        else:
            att = self.input_scw.auxadppath + "/attitude_snapshot.fits[AUXL-ATTI-SNA,1,BINTABLE]"
            attp_fn = glob.glob(self.input_scw.auxadppath + "/attitude_predicted_*.fits*")[0]
            attp = attp_fn + "[AUXL-ATTI-PRE,1,BINTABLE]"

        construct_gnrl_scwg_grp(self.input_scw, [ \
            self.input_gb.corshad.path,
            self.input_scw.auxadppath + "/time_correlation.fits[AUXL-TCOR-HIS]",
            att,
            attp,
            self.input_gti.output_gti.path
        ])
        # self.input_cat.cat.path,

        import_attr(self.input_scw.scwpath + "/swg.fits",
                    ["OBTSTART", "OBTEND", "TSTART", "TSTOP", "SW_TYPE", "TELAPSE", "SWID"])
        set_attr({'ISDCLEVL': self.shdtype})
        # set_attr({'INSTRUME':"IBIS"},"og.fits")

        # construct_gnrl_scwg_grp_idx(self.input_scw,[\
        #            "og.fits",
        #        ])
        # set_attr({'ISDCLEVL':"BIN_I"},"og_idx.fits")

        #  construct_og(self.input_scw,[\
        #              "og_idx.fits",
        #          ])
        #  set_attr({'ISDCLEVL':"BIN_I"},"ogg.fits")

        # remove_withtemplate("isgri_srcl_res.fits(ISGR-SRCL-RES.tpl)")

        #pif_fn, pif_tpl = "isgri_pif.fits", "(ISGR-PIF.-SHD-IDX.tpl)"
        lc_fn, lc_tpl = "isgri_lcr.fits", "(ISGR-SRC.-LCR-IDX.tpl)"

        remove_withtemplate(lc_fn + lc_tpl)

        ht = heatool(self.binary)
        ht['outSwg'] = "og.fits"
        ht['inCat'] = self.input_cat.cat.get_path()
        ht['outLC'] = lc_fn + lc_tpl
        ht['mask'] = self.input_ic.ibisicroot + "/mod/isgr_mask_mod_0003.fits[ISGR-MASK-MOD,1,IMAGE]"
        ht['tungAtt'] = self.input_ic.ibisicroot + "/mod/isgr_attn_mod_0010.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['aluAtt'] = self.input_ic.ibisicroot + "/mod/isgr_attn_mod_0011.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['leadAtt'] = self.input_ic.ibisicroot + "/mod/isgr_attn_mod_0012.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['isgrUnifDol'] = self.input_maps.unif.path
        if self.usebkg:
            ht['isgrBkgDol'] = self.input_maps.back.path
        else:
            ht['isgrBkgDol'] = "-"
        ht['corrDol'] = self.input_maps.corr.path
        ht['OutType'] = self.shdtype
        ht['method_cor'] = 1

        # for k in ['SearchMode','ToSearch','CleanMode','MinCatSouSnr','MinNewSouSnr','NegModels','DoPart2']: # dopart2 is flow control, separately
        #    ht[k]=getattr(self.input_imgconfig,k)

        ht.run()

        self.lightcurve = DataFile(lc_fn)



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

    CFITSIO_INCLUDE_FILES_shortened = re.sub("\:+",":", os.environ['CFITSIO_INCLUDE_FILES'])
    if os.environ['CFITSIO_INCLUDE_FILES'] != CFITSIO_INCLUDE_FILES_shortened:
        print("shortened CFITSIO_INCLUDE_FILES:", os.environ['CFITSIO_INCLUDE_FILES'] ,"to",CFITSIO_INCLUDE_FILES_shortened)

    os.environ['CFITSIO_INCLUDE_FILES'] = CFITSIO_INCLUDE_FILES_shortened

    dc=heatool("dal_create")
    dc['obj_name']="!"+fn
    dc['template']="GNRL-SCWG-GRP.tpl"

    try:
        dc.run()
    except Exception as e:
        print("unknown exception",e)
        print(subprocess.check_output("export",shell=True))
        print(subprocess.check_output("echo $CFITSIO_INCLUDE_FILES/GNRL-SCWG-GRP.tpl",shell=True))
        print(subprocess.check_output("ls -l $CFITSIO_INCLUDE_FILES/GNRL-SCWG-GRP.tpl",shell=True))
        print(subprocess.check_output("pwd",shell=True))
        print(subprocess.check_output("ls -lort",shell=True))
        raise
    
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
            for children_group_id in range(int(len(children)/4+1)):
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

    if len(children)==0:
        raise NoValidScW()

    open("swgs.txt","w").write("\n".join(children))
    dc=heatool("txt2idx")
    dc['element']="swgs.txt"
    dc['index']=fn
    dc['template']="GNRL-SCWG-GRP-IDX.tpl"
    dc.run()

def construct_gnrl_arbitrary_grp_idx(children=[],structure="GNRL-SCWG-GRP",fn="og_idx.fits"):
    remove_withtemplate(fn+"("+structure+"-IDX.tpl)")
    #remove_withtemplate(fn + "(" + structure + "-IDX.tpl)")

    open("swgs.txt","w").write("\n".join(children))
    dc=heatool("txt2idx")
    dc['element']="swgs.txt"
    dc['index']=fn
    dc['template']=structure+"-IDX.tpl"
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

    pt2dt={int:"DAL_INT",str:"DAL_CHAR", float:'DAL_DOUBLE'}
    pt2k={int:"i",str:"s", float:'r'}

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
    
    og=fits.open(fn) 
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
        print("appending",e1,e2)
        dc=heatool('dal_append')
        dc['grpDOL']=fn
        dc['element']='ISGR-DETE-SHD.tpl'
        dc.run()

    og=fits.open(fn) 
    for i,(e1,e2) in enumerate(bins):
        og[1].data[i]['E_MIN']=e1
        og[1].data[i]['E_MAX']=e2
        og[1].data[i]['ISDCLEVL']=levl
        og[2+i].header['E_MIN']=e1
        og[2+i].header['E_MAX']=e2
        og[2+i].header['ISDCLEVL']=levl
    og.writeto(fn,clobber=True)

#class AnyScW(da.AnyAnalysis):
#    pass

#class AnyRev(da.AnyAnalysis):
#    pass



class IDScWList(DataAnalysis):
    scwid_list=None
    allow_alias=True

    def get_version(self):
        v=self.get_signature()+"."+self.version

        if not isinstance(self.scwid_list,list):
            return v+repr(self.scwid_list)
            #raise Exception("scwid_list must be a list")

        if len(self.scwid_list)==1:
            v+="_one_"+self.scwid_list[0]
            return v
        
        v=("_n%i"%len(self.scwid_list))+"_"+shhash("_".join(self.scwid_list))[:4]

        revs=sorted(set([s[:4] for s in self.scwid_list]))
        if len(revs)==1:
            v+="_r..."+revs[0]
            return v

        v += "_r...from_" + revs[0] + "...to_" + revs[-1] 
        return v

    def main(self):
        self.scwlistdata=[ScWData(input_scwid=s.strip()) for s in self.scwid_list]


class RevScWList(DataAnalysis):
    input_rev=Revolution

    run_for_hashe=False
    allow_alias=True

    def main(self):
        import os

        event_files=[]
        for event_file in glob.glob(self.input_rev.revroot+"/*/isgri_events.fits*"):
            print(event_file)
            try:
                evts=fits.open(event_file)['ISGR-EVTS-ALL']
                print(evts)
                if evts.data.shape[0]<10000:
                    raise Exception("very few events!")
                event_files.append(event_file)
            except Exception as e:
                print(e)

        print("event files in",self.input_rev.revroot)
        print("found event files",event_files)

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
        print("available list",len(thelist),thelist)
        print("first scw:",self.firstscw)
        if self.firstscws:
            self.scwlistdata=thelist[self.firstscw:self.nscw+self.firstscw]
        else:
            self.scwlistdata=thelist[-self.nscw:]
        if self.step!=1:
            self.scwlistdata=self.scwlistdata[::self.step]
        print("resulting list:",len(self.scwlistdata),self.scwlistdata)

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

def converttime(informat,intime,outformat=""):
    r=subprocess.check_output(["converttime",informat,intime,outformat])
    d={}
    for l in r.split("\n"):
        t=re.search("Output Time\((.*?)\): (.*?)$",l,re.S)
        #print(l,t                                                                                                                                                                                             )
        if t:
            g=t.groups()
            if 'Boundary' in g[1]:
                d[g[0]]=g[1].split()[1:]
            else:
                d[g[0]]=g[1]

    if outformat == "":
        return d
    else:
        return d[outformat]

def fromUTC(utc):
    r=subprocess.check_output(["converttime","UTC",utc,""])
    d={}
    for l in r.split("\n"):
        t=re.search("Output Time\((.*?)\): (.*?)$",l,re.S)                                                                                                                                                     
        #print(l,t                                                                                                                                                                                             )
        if t:                                                                                                                                                                                                  
            g=t.groups()                                                                                                                                                                                       
            d[g[0]]=g[1]                                                                                                                                                                                       
            #print(g                                                                                                                                                                                           )
    return d    


import dataanalysis.callback

previously_accepted_classes=dataanalysis.callback.default_callback_filter.callback_accepted_classes

class CallbackRareDDOSAFilter(dataanalysis.callback.Callback):
    def extract_data(self,obj):
        data={'scwid':'inapplicable',}

        scw=obj.cache.get_scw(obj._da_locally_complete)
        
        expected_hashe=getattr(obj,'_da_expected_full_hashe',None)

        if expected_hashe is not None:
            data['node_id']=obj.cache.hashe2signature(expected_hashe)
        else:
            data['node_id']="undefined_expected_hashe_please_complain" # add sentry

        if scw is None:
            scw=obj.cache.get_scw(expected_hashe)

        if scw is not None:
            data.update({"scwid":scw})

        return data

dataanalysis.callback.default_callback_filter=CallbackRareDDOSAFilter

if previously_accepted_classes is not None:
    dataanalysis.callback.default_callback_filter.set_callback_accepted_classes(previously_accepted_classes)

dataanalysis.callback.default_callback_filter.set_callback_accepted_classes([mosaic_ii_skyimage, ii_skyimage, BinEventsImage, ibis_gti, ibis_dead, ISGRIEvents, ii_spectra_extract, BinEventsSpectra, ii_lc_extract, BinEventsLC, lc_pick])

