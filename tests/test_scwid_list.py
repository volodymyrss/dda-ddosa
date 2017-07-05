
def test_scwid_list():
    import ddosa
    import dataanalysis.core as da

    da.debug_output()

    fa=ddosa.IDScWList(use_scwid_list=["066500230010.001"])

    print(fa.get())
    print(fa.get_version())

    assert len(fa.scwlistdata) == 1

    fa=ddosa.IDScWList(use_scwid_list=["066500230010.001","066500250010.001"])

    print(fa.get())
    print(fa.get_version())

    assert len(fa.scwlistdata) == 2

    fa=ddosa.IDScWList(use_scwid_list=["066500230010.001","066600250010.001"])

    print(fa.get())
    print(fa.get_version())

    assert len(fa.scwlistdata) == 2
