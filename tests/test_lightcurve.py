import glob
import os

import astropy.io.fits as fits

import ddosa
import dataanalysis.core as da


def test_lightcurves():
    #da.debug_output()

    fa = ddosa.ii_lc_extract(assume=[
        ddosa.ScWData(input_scwid="066500230010.001"),
    ])
    fa.read_caches = []

    fa.get()

    #assert os.path.exists(fa.spectrum.get_path())
