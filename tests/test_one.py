import glob

import astropy.io.fits as fits

import ddosa



def test_one():
    fa=ddosa.ii_skyimage(assume=[
                            ddosa.ScWData(input_scwid="066500220010.001"),
                        ])
    fa.get()
