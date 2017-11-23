def test_spectra():
    import dataanalysis.core as da
    #    da.debug_output()
    da.reset()
    import ddosa
    reload(ddosa)

    fa = ddosa.SpectraProcessingSummary()
    fa.read_caches = []

    fa.get()