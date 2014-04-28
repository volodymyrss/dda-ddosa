
from dataanalysis import DataAnalysis,DataFile,MemCacheLocal

from pilton import heatool

import os

class ScWData(DataAnalysis):
    input_scwid=None

    def __init__(self,scwid=None):
        self.input_scwid=scwid

    def main(self):
        self.scwpath=os.environ['REP_BASE_PROD']+"/scw/"+self.input_scwid[:4]+"/"+self.input_scwid

class ibis_isgr_energy(DataAnalysis):
    input_scw=None

    version="v2"
   
    def main(self):
        dc=heatool("dal_create")
        dc['obj_name']="og.fits"
        dc['template']="GNRL-SCWG-GRP.tpl"

        bin="/afs/in2p3.fr/throng/integral/software/spectral/ibis_isgr_energy/ibis_isgr_energy_pha2/ibis_isgr_energy"
        ht=heatool(bin)
        ht['inGRP']=""
        ht['OutCorEvts']="isgri_events_corrected.fits"
        #ht['inGRP']
        ht.run()


