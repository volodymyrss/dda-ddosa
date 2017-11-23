def test_spectra():
    import dataanalysis.core as da
    #    da.debug_output()
    da.reset()

    import os
    os.environ['INTEGRAL_DDCACHE_ROOT'] = "./"

    import ddosa
    reload(ddosa)


    class Analysis(ddosa.DataAnalysis):
        cached=True
        input_fact=ddosa.SpectraProcessingSummary

    fa = Analysis()
    fa.read_caches = []

    fa.get()

    print fa

    assert len(fa.factory.factorizations)==1