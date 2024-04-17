import webbpsf

filtername = 'F2550W'

mir = webbpsf.MIRI()
mir.filter = filtername

grid = mir.psf_grid(num_psfs=16, oversample=1, all_detectors=False, verbose=True, save=True)