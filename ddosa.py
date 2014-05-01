# check for logic in osa! make sure there is none!

from dataanalysis import DataFile
#from dataanalysis import DataAnalysis,DataFile,Cache
import dataanalysis 

from pilton import heatool

import os,shutil,re,time

mcg=dataanalysis.MemCache('/sps/integral/analysis/cache/dda_global')

class DataAnalysis(dataanalysis.DataAnalysis):
    cache=mcg

#
# there are two ways to define common origin:
#   define during class declaratopm
# dynamic redefinition of the classes, like ScWData, will invalidate (hopefully) the get call, so it should be reexecuted

# there is special treatement
class ScWData(DataAnalysis):
    input_scwid=None

    def __init__(self,scwid=None):
        if scwid is not None:
            self.input_scwid=scwid

    def main(self):
        self.scwpath=os.environ['REP_BASE_PROD']+"/scw/"+self.input_scwid[:4]+"/"+self.input_scwid
        self.revdirpath=os.environ['REP_BASE_PROD']+"/scw/"+self.input_scwid[:4]+"/rev.001"

        if not os.path.exists(self.scwpath+"/swg.fits"):
            raise Exception("no scw data!")

class GetLUT2(DataAnalysis):
    input="LUT2_standard"

    def main(self):
        self.datafile=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_3dl2_mod_0001.fits"
        #self.datafile=DataFile(os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_3dl2_mod_0001.fits")

class GetEcorrCalDB(DataAnalysis):
    input=["ecorr_standard_OSA10"]
    input_lut2=GetLUT2()

    def main(self):
        self.godol=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/cal/ibis_isgr_gain_offset_0010.fits"
        self.supgdol=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_gain_mod_0001.fits[ISGR-GAIN-MOD,1,BINTABLE]"
        self.supodol=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_off2_mod_0001.fits[ISGR-OFF2-MOD,1,BINTABLE]"
        self.risedol=self.input_lut2.datafile

class ibis_isgr_energy(DataAnalysis):
    
    input_scw=ScWData()
    #input_raw_events=None
    #input_ibis_hk=None
    # can also take only isgri events

    input_ecorrdata=GetEcorrCalDB()

    version="v3"
   
    def main(self):

        remove_withtemplate("isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)")
    
        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]" \
        ])

        bin="/afs/in2p3.fr/throng/integral/software/spectral/ibis_isgr_energy/ibis_isgr_energy_pha2/ibis_isgr_energy"
        ht=heatool(bin)
        ht['inGRP']="og.fits"
        ht['outCorEvts']="isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)"
        ht['useGTI']="n"
        ht['riseDOL']=self.input_ecorrdata.risedol
        ht['GODOL']=self.input_ecorrdata.godol
        ht['supGDOL']=self.input_ecorrdata.supgdol
        ht['supODOL']=self.input_ecorrdata.supodol
        ht['chatter']="4"
        ht.run()


        self.output_events=DataFile("isgri_events_corrected.fits")
            

class ibis_isgr_evts_tag(DataAnalysis):
    
    input_events_corrected=ibis_isgr_energy()
    input_scw=ScWData() # again, can get separate

    #input_ecorrdata=GetEcorrCalDB()

    version="v0"
   
    def main(self):

        #remove_withtemplate("isgri_events_corrected.fits(ISGR-EVTS-COR.tpl)")
    
        construct_gnrl_scwg_grp(self.input_scw,[\
            self.input_scw.scwpath+"/isgri_events.fits[3]", \
            self.input_events_corrected.output_events.path \
        ])
            #self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]" \

        import_attr(self.input_scw.scwpath+"/swg.fits",['OBTSTART','OBTEND'])
        
        bin="ibis_isgr_evts_tag"
        ht=heatool(bin)
        ht['inGRP']="og.fits"
        ht['idxSwitch']=self.input_scw.revdirpath+"/idx/isgri_pxlswtch_index.fits[1]"
        ht['seleEXT']="ISGR-EVTS-COR"
        ht.run()

        self.output_events=DataFile("isgri_events_corrected.fits")

class ICRoot(DataAnalysis):
    input="standard_IC"

    version="v1"

    def main(self):
        self.icroot=os.environ['REP_BASE_PROD']+"/ic/ic_10/"
        self.icindex=os.environ['REP_BASE_PROD']+"/ic/ic_10/idx/ic/ic_master_file.fits[1]"


class IBIS_ICRoot(DataAnalysis):
    input="standard IC"

    def main(self):
        self.ibisicroot=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis"

# maybe split indeed,but try to show another case
class ibis_gti(DataAnalysis):
    input_scw=ScWData()
    input_ic=ICRoot()
    #input_gticreate=gti_create()
    
    version="v2"
    def main(self):
        # horrible horrible full OSA

        open("scw.list","w").write(self.input_scw.input_scwid)

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))
        
        bin="og_create"
        ogc=heatool(bin)
        ogc['idxSwg']="scw.list"
        ogc['instrument']="IBIS"
        ogc['ogid']="scw_"+self.input_scw.input_scwid[:12]
        ogc['baseDir']="./"
        ogc.run()
        
        scwroot="scw/"+self.input_scw.input_scwid

        bin="ibis_gti"
        ht=heatool(bin,wd="obs/"+ogc['ogid'].value) 
        ht['swgDOL']=scwroot+"/swg_ibis.fits"
        ht['GTI_Index']="ibis_gti.fits(IBIS-GNRL-GTI-IDX.tpl)"
        ht['IC_Group']=self.input_ic.icindex
        ht['IC_Alias']="OSA"
        ht['disablePICsIT']="YES"
        ht['GTI_attTolerance_X']="0.05"
        ht['GTI_attTolerance_Z']="0.2"
        ht['GTI_BTI_Dol']=self.input_ic.icroot+"/ic/ibis/lim/isgr_gnrl_bti_0012.fits"
        ht['GTI_BTI_Names']="IBIS_CONFIGURATION IBIS_BOOT ISGRI_RISE_TIME VETO_PROBLEM SOLAR_FLARE BELT_CROSSING"
        ht.run()

        shutil.copy(ht.cwd+"/ibis_gti.fits","./ibis_gti.fits")
        self.output_gti=DataFile("ibis_gti.fits")

# maybe split indeed,but try to show another way
class ibis_dead(DataAnalysis):
    input_scw=ScWData()
    input_ic=ICRoot()
    
    version="v2"
    def main(self):
        # horrible horrible full OSA

        open("scw.list","w").write(self.input_scw.input_scwid)

        if os.path.exists("obs"):
            os.rename("obs","obs."+str(time.time()))
        
        bin="og_create"
        ogc=heatool(bin)
        ogc['idxSwg']="scw.list"
        ogc['instrument']="IBIS"
        ogc['ogid']="scw_"+self.input_scw.input_scwid[:12]
        ogc['baseDir']="./"
        ogc.run()
        
        scwroot="scw/"+self.input_scw.input_scwid

        bin="ibis_dead"
        ht=heatool(bin,wd="obs/"+ogc['ogid'].value) 
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

class root(DataAnalysis):
    input=[ibis_gti()]

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
    

def construct_gnrl_scwg_grp(scw,children):
    dc=heatool("dal_create")
    dc['obj_name']="!og.fits"
    dc['template']="GNRL-SCWG-GRP.tpl"
    dc.run()
    
    da=heatool("dal_attr")
    da['indol']="og.fits"
    da['keynam']="REVOL"
    da['action']="WRITE"
    da['type']="DAL_INT"
    da['value_i']=scw.input_scwid[:4]
    da.run()
    
    if children!=[]:
        dac=heatool("dal_attach")
        dac['Parent']="og.fits"
        for i,c in enumerate(children):
            dac['Child%i'%(i+1)]=c
            if i>3:
                raise Exception("can not attach more than 4 children to the group!")
        dac.run()

def import_attr(obj,attr):
    da=heatool("dal_attr_copy")
    da['indol']=obj
    da['outdol']="og.fits"
    da['keylist']=",".join(attr)
    da.run()
