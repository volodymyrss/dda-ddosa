import glob

import astropy.io.fits as fits

import ddosa
import dataanalysis.core as da


def test_scw():
    da.debug_output()

    
    fa=ddosa.ScWData(input_scwid="066500230010.001")

    fa.get()
    print fa.scwpath


def test_gti():
    da.debug_output()

    return # disabled cache

    fa=ddosa.ibis_gti(assume=[
                            ddosa.ScWData(input_scwid="066500230010.001"),
                        ])
    fa.read_caches=[]
    

    fa.get()


def test_gti_cached():
    da.debug_output()

    fa = ddosa.ibis_gti(assume=[
        ddosa.ScWData(input_scwid="066500230010.001"),
    ])

    fa.get()


    # now get cached

    fb = ddosa.ibis_gti(assume=[
        ddosa.ScWData(input_scwid="066500230010.001"),
    ])

    fb.produce_disabled=True

    fb.get()


def test_one():
    da.debug_output()

    
    fa=ddosa.ii_skyimage(assume=[
                            ddosa.ScWData(input_scwid="066500230010.001"),
                        ])
    fa.read_caches=[]

    fa.get()


    print fa.skyima

