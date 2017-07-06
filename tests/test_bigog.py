import ddosa
import os

def test_construct_many_images():
    ig=ddosa.ImageGroups(input_scwlist=ddosa.IDScWList(use_scwid_list=["066500230010.001","066500240010.001"]))
    ig=ig.get()

    print ig

    assert len(ig.members)==2
    assert len(ig.members[0]) == 4

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



