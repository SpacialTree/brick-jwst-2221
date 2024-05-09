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

#reg_hmsfr = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions/hmsfr.reg')


