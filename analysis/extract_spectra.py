import matplotlib.pyplot as plt
#import pyspeckit 

from astropy.io import fits
from astropy.wcs import WCS

from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
from regions import CircleSkyRegion

from spectral_cube import SpectralCube
from spectral_cube import Projection

def extract_spec(cube, pos, rad):
    circle_region = CircleSkyRegion(center=pos, radius=rad)
    subcube = cube.subcube_from_regions([circle_region], minimize=False)
    subcube.allow_huge_operations = True
    subcube = subcube.to(u.K)
    average_spectrum = subcube.mean(axis=(1, 2))
    return average_spectrum

pos_cloudc1 = SkyCoord('17:46:21.3048683891', '-28:35:33.1211282499', unit=(u.hourangle, u.deg))
pos_cloudc2 = SkyCoord('17:46:18.3316118680', '-28:34:48.4717811920', unit=(u.hourangle, u.deg))
pos_cloudd = SkyCoord('17:46:22.6563371259', '-28:33:27.5405071803', unit=(u.hourangle, u.deg))
pos_filament = SkyCoord('17:46:20.9063719501', '-28:37:51.6942550990', unit=(u.hourangle, u.deg))

rad_cloudc1 = 30*u.arcsec
rad_cloudc2 = 30*u.arcsec
rad_cloudd = 60*u.arcsec
rad_filament = 1.5*u.arcmin

cube = SpectralCube.read('/orange/adamginsburg/ACES/mosaics/cubes/HNCO_7m12mTP_CubeMosaic.fits')
#cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=87.925238*u.GHz)

# Cloud c1
average_spectrum_cloudc1 = extract_spec(cube, pos_cloudc1, rad_cloudc1)
average_spectrum_cloudc1.write('/orange/adamginsburg/jwst/cloudc/alma/spectra/spectrum_cloudc1_HNCO.fits')

average_spectrum_cloudc2 = extract_spec(cube, pos_cloudc2, rad_cloudc2)
average_spectrum_cloudc2.write('/orange/adamginsburg/jwst/cloudc/alma/spectra/spectrum_cloudc2_HNCO.fits')

average_spectrum_cloudd = extract_spec(cube, pos_cloudd, rad_cloudd)
average_spectrum_cloudd.write('/orange/adamginsburg/jwst/cloudc/alma/spectra/spectrum_cloudd_HNCO.fits')

average_spectrum_filament = extract_spec(cube, pos_filament, rad_filament)
average_spectrum_filament.write('/orange/adamginsburg/jwst/cloudc/alma/spectra/spectrum_filament_HNCO.fits')

