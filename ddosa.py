
from dataanalysis import DataAnalysis,DataFile,MemCache

from pilton import heatool

import os,shutil,shutil

mcg=MemCache('/sps/integral/analysis/cache/dda_global')

class ScWData(DataAnalysis):
    cache=mcg
    input_scwid=None

    def __init__(self,scwid=None):
        self.input_scwid=scwid

    def main(self):
        self.scwpath=os.environ['REP_BASE_PROD']+"/scw/"+self.input_scwid[:4]+"/"+self.input_scwid

class GetLUT2(DataAnalysis):
    cache=mcg
    input="LUT2_standard"

    def main(self):
        self.datafile=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_3dl2_mod_0001.fits"
        #self.datafile=DataFile(os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_3dl2_mod_0001.fits")

class GetEcorrCalDB(DataAnalysis):
    cache=mcg
    input=["standard_OSA10"]
    input_lut2=GetLUT2()
    def main(self):
        self.godol=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/cal/ibis_isgr_gain_offset_0010.fits"
        self.supgdol=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_gain_mod_0001.fits[ISGR-GAIN-MOD,1,BINTABLE]"
        self.supodol=os.environ['REP_BASE_PROD']+"/ic/ic_10/ic/ibis/mod/isgr_off2_mod_0001.fits[ISGR-OFF2-MOD,1,BINTABLE]"
        self.risedol=self.input_lut2.datafile

class ibis_isgr_energy(DataAnalysis):
    cache=mcg
    input_scw=None
    # can also take only isgri events

    input_ecorrdata=GetEcorrCalDB()

    version="v3"
   
    def main(self):
        dc=heatool("dal_create")
        dc['obj_name']="!og.fits"
        dc['template']="GNRL-SCWG-GRP.tpl"
        dc.run()
        
        da=heatool("dal_attr")
        da['indol']="og.fits"
        da['keynam']="REVOL"
        da['action']="WRITE"
        da['type']="DAL_INT"
        da['value_i']=self.input_scw.input_scwid[:4]
        da.run()
        
        dac=heatool("dal_attach")
        dac['Parent']="og.fits"
        dac['Child1']=self.input_scw.scwpath+"/isgri_events.fits[3]" # get properly
        dac['Child2']=self.input_scw.scwpath+"/ibis_hk.fits[IBIS-DPE.-CNV]" # get properly
        dac.run()

        
        # function
        try:
            os.remove("isgri_events_corrected.fits")
        except OSError:
            pass
        try:
            os.remove("ISGR-EVTS-COR.tpl")
        except OSError:
            pass

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


