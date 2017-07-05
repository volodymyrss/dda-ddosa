
def test_scwid_list():
    import ddosa

    fa=ddosa.IDScWList(use_scwid_list=["066500230010.001"])

    print(fa.get())
    print(fa.get_version())

    fa=ddosa.IDScWList(use_scwid_list=["066500230010.001","066500250010.001"])

    print(fa.get())
    print(fa.get_version())

    fa=ddosa.IDScWList(use_scwid_list=["066500230010.001","066600250010.001"])

    print(fa.get())
    print(fa.get_version())