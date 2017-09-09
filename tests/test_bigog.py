import ddosa
import os

def test_construct_many_images():
    ig=ddosa.ImageGroups(input_scwlist=ddosa.IDScWList(use_scwid_list=["066500230010.001","066500240010.001"]))
    ig=ig.get()

    print ig

    assert len(ig.members)==2
    assert len(ig.members[0]) == 5

    assert hasattr(ig.members[0][1],'skyima')
    assert os.path.exists(ig.members[0][1].skyima.get_path())

    ig.construct_og("og.fits")

def test_mosaic_ii_skyimage():
    #da.debug_output()

    ig = ddosa.ImageGroups(input_scwlist=ddosa.IDScWList(use_scwid_list=["066500230010.001", "066500240010.001"]))

    fa = ddosa.mosaic_ii_skyimage(input_imagegroups=ig)
    fa.read_caches = []

    fa.get()

#    assert os.path.exists(fa.lightcurve.get_path())


def test_lc_pick():
    #da.debug_output()
    from dataanalysis import analysisfactory

    analysisfactory.AnalysisFactory.WhatIfCopy("timebin",ddosa.LCTimeBin(use_time_bin_seconds=200))
    ig = ddosa.LCGroups(input_scwlist=ddosa.IDScWList(use_scwid_list=["066500230010.001", "066500240010.001"]))

    fa = ddosa.lc_pick(input_lcgroups=ig,use_source_names=["Crab"])
    fa.read_caches = []

    fa.get()

#    assert os.path.exi