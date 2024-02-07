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

# Filename for large cube
filename_full = '/orange/adamginsburg/ACES/mosaics/cubes/HNCO_7m12mTP_CubeMosaic.fits'
line = "HNCO" ## CHANGE THIS FOR OTHER LINES

# All of the Cloud C cutouts
files = glob.glob('/orange/adamginsburg/jwst/cloudc/alma/cutouts/cloudc_*')

# Read the region files
print('Reading region files.')
reg = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/nircam_cloudc_fov_perp.reg')
reg_cloud = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/large_cloud.reg')
reg_sub = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/sub_cloud.reg')
reg_hmsfr = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/hmsfr.reg')
reg_filament = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/filament.reg')
reg_double = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/double_peak.reg')

reg_fov = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/fov.reg')

### Cubes
print('Opening Cubes')
cube_full = SpectralCube.read(filename_full).with_spectral_unit(u.km/u.s, velocity_convention='radio')

cube_full_reg = cube_full.subcube_from_regions(reg, minimize=False)
cube_full_cloud = cube_full.subcube_from_regions(reg_cloud, minimize=False).spectral_slab(8*u.km/u.s, 37*u.km/u.s)
cube_full_sub = cube_full.subcube_from_regions(reg_sub, minimize=False).spectral_slab(-27*u.km/u.s, 23*u.km/u.s)
cube_full_hmsfr = cube_full.subcube_from_regions(reg_hmsfr, minimize=False).spectral_slab(25*u.km/u.s, 50*u.km/u.s)
cube_full_filament = cube_full.subcube_from_regions(reg_filament, minimize=False).spectral_slab(-59*u.km/u.s, -52*u.km/u.s)
cube_full_double = cube_full.subcube_from_regions(reg_double, minimize=False).spectral_slab(0*u.km/u.s, 45*u.km/u.s)

### Moment Maps - put these in order of least and most likely to crash
print('Making moment map for double peaked cloud.')
mom0_full_double = cube_full_double.moment0()
mom0_full_double.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_{line}_double.fits')

print('Making moment map for smaller sub cloud.')
mom0_full_sub = cube_full_sub.moment0()
mom0_full_sub.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_{line}_sub.fits')

print('Making moment map for filament.')
mom0_full_filament = cube_full_filament.moment0()
mom0_full_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_{line}_filament.fits')

print('Making moment map for high mass star formation region.')
mom0_full_hmsfr = cube_full_hmsfr.moment0()
mom0_full_hmsfr.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_{line}_hmsfr.fits')

print('Making moment map for larger cloud.')
mom0_full_cloud = cube_full_cloud.moment0()
mom0_full_cloud.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_{line}_cloud.fits')

print('Making moment map for full fov of jwst data.')
mom0_full_reg = cube_full_reg.moment0()
mom0_full_reg.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_{line}_reg.fits')

### Attempting to write the cubes
cube_full_double.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_{line}_double.fits')
cube_full_sub.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_{line}_sub.fits')
cube_full_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_{line}_filament.fits')
cube_full_hmsfr.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_{line}_hmsfr.fits')
cube_full_cloud.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_{line}_cloud.fits')
cube_full_reg.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_{line}_reg.fits')
