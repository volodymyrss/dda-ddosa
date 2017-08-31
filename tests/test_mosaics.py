
def test_mosaic_ii_skyimage():
    import ddosa

    mosaic=ddosa.mosaic_ii_skyimage(
              assume=[
                  ddosa.ImageGroups(input_scwlist=ddosa.IDScWList(use_scwid_list=["066500330010.001","066500340010.001"])),
                  ddosa.ImageBins(use_ebins=[(25, 40)], use_version="onebin_25_40"),
                  ddosa.ImagingConfig(use_SouFit=0, use_version="soufit0")
              ]
            )

    mosaic.get()

