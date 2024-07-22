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

def extract_spectrum(obj, fn, line_name='HNCO', restfreq=87.925238*u.GHz):
    cube = SpectralCube.read(fn)
    cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=restfreq)
    pos = SkyCoord(obj['ra'], obj['dec'], unit=(u.hourangle, u.deg))
    rad = obj['rad']*u.arcsec
    average_spectrum = extract_spec(cube, pos, rad)
    average_spectrum.write(f'/orange/adamginsburg/jwst/cloudc/alma/spectra/spectrum_{obj["name"]}_{line_name}.fits')
    #return average_spectrum

cloudc1 = {'name': 'cloudc1', 'ra': '17:46:21.3048683891', 'dec': '-28:35:33.1211282499', 'rad': 30}
cloudc2 = {'name': 'cloudc2', 'ra': '17:46:18.3316118680', 'dec': '-28:34:48.4717811920', 'rad': 30}
cloudd = {'name': 'cloudd', 'ra': '17:46:22.6563371259', 'dec': '-28:33:27.5405071803', 'rad': 60}
filament = {'name': 'filament', 'ra': '17:46:20.9063719501', 'dec': '-28:37:51.6942550990', 'rad': 90}

fn_HNCO = '/orange/adamginsburg/ACES/mosaics/cubes/HNCO_7m12mTP_CubeMosaic.fits'
fn_12CO = '/orange/adamginsburg/cmz/nobeyama/12CO-2.BEARS.FITS'
#cube = cube.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=87.925238*u.GHz)

# Cloud c1
extract_spectrum(cloudc1, fn_HNCO, line_name='HNCO', restfreq=87.925238*u.GHz)
extract_spectrum(cloudc1, fn_12CO, line_name='12CO', restfreq=1.152712040000E+11*u.Hz)

# Cloud c2
extract_spectrum(cloudc2, fn_HNCO, line_name='HNCO', restfreq=87.925238*u.GHz)
extract_spectrum(cloudc2, fn_12CO, line_name='12CO', restfreq=1.152712040000E+11*u.Hz)

# Cloud d
extract_spectrum(cloudd, fn_HNCO, line_name='HNCO', restfreq=87.925238*u.GHz)
extract_spectrum(cloudd, fn_12CO, line_name='12CO', restfreq=1.152712040000E+11*u.Hz)

# Filament
extract_spectrum(filament, fn_HNCO, line_name='HNCO', restfreq=87.925238*u.GHz)
extract_spectrum(filament, fn_12CO, line_name='12CO', restfreq=1.152712040000E+11*u.Hz)


