from spectral_cube import SpectralCube
import numpy as np
from astropy.visualization import quantity_support
from astropy import units as u
from astropy import wcs
import matplotlib.pyplot as plt
from astropy.utils import data
from reproject import reproject_exact
from astropy.io import fits
import regions
import glob, os
import spectral_cube

vlow = -60*u.km/u.second
vhigh = -50*u.km/u.second

reg = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/filament_reg.reg')

fn_HNCO_mosaic = '/orange/adamginsburg/ACES/mosaics/cubes/HNCO_7m12mTP_CubeMosaic.fits'
cube_full = SpectralCube.read(filename_full).with_spectral_unit(u.km/u.s, velocity_convention='radio')

cube_full_filament = cube_full.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)
cube_full_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_HNCO_7m12mTP_filament.fits')

mom0_filament = cube_full_filament.moment0()
mom0_full_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_HNCO_7m12mTP_filament.fits')



