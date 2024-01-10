from scipy import ndimage
from astropy.io import fits
import scipy.ndimage
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'
pl.rcParams['image.origin'] = 'lower'
pl.rcParams['figure.figsize'] = (10,8)
from astropy.io import fits
import reproject
from astropy import convolution
from astropy.convolution import Gaussian2DKernel
from astropy import units as u
from astropy.table import Table
import pyavm
import regions
from astropy import coordinates
import PIL
#from spectral_cube import SpectralCube, Projection, Slice
from astropy import wcs
import matplotlib.pyplot as plt
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb
import numpy as np
import pylab as pl
from astropy.visualization import simple_norm
import glob
import warnings
import regions
from tqdm.notebook import tqdm
from astropy import log

basepath = '/orange/adamginsburg/jwst/cloudc/'

import sys
sys.path.append('/orange/adamginsburg/jwst/cloudc/reduction/brick-jwst-2221/reduction/')

from destreak import compute_zero_spacing_approximation

'''
for fn in glob.glob(f'{basepath}/images/jw02221-o002_t001_nircam_*nodestreak_i2d.fits'):
    print(fn)
    with warnings.catch_warnings():
        # specifically ignoring that we're using non-integer 
        warnings.simplefilter('ignore')
        dx = 128 if 'f4' in fn else 256
        sl = 'long' if 'f4' in fn else 'short'
        regs = regions.Regions.read(f'{basepath}/regions_/bright_stars_{sl}.reg')
        tp = compute_zero_spacing_approximation(fn, ext=1, dx=dx, regs=regs, percentile=10, progressbar=tqdm)
        #if np.nanpercentile(tp.data, 10) < 0:
        #    print(np.nanpercentile(tp.data, 10))
        #    tp.data -= np.nanpercentile(tp.data, 10)
        #tp.data[tp.data<0] = 0
        tp.writeto(fn.replace("_i2d.fits", f"_i2d_medfilt{dx}.fits"), overwrite=True)
        print(fn)
        print(fn.replace("_i2d.fits", f"_i2d_medfilt{dx}.fits"))
'''

log.info('Making Medium Filter Background without stars masked, percentile version, dx = 64 if "f4" in fn else 128')
for fn in glob.glob(f'{basepath}/images/jw02221-o002_t001_nircam_*nodestreak_i2d.fits'):
    log.info(fn)
    with warnings.catch_warnings():
        # specifically ignoring that we're using non-integer 
        warnings.simplefilter('ignore')
        dx = 64 if 'f4' in fn else 128 #64 # 128 if 'f4' in fn else 256
        sl = 'long' if 'f4' in fn else 'short'
        #regs = regions.Regions.read(f'{basepath}/regions_/bright_stars_{sl}.reg')
        log.info('Computing zero spacing approxmation.')
        tp = compute_zero_spacing_approximation(fn, ext=1, dx=dx, percentile=10, progressbar=tqdm)
        #if np.nanpercentile(tp.data, 10) < 0:
        #    print(np.nanpercentile(tp.data, 10))
        #    tp.data -= np.nanpercentile(tp.data, 10)
        #tp.data[tp.data<0] = 0
        log.info(f'Writing to {fn.replace("_i2d.fits", f"_i2d_nomask-medfilt{dx}.fits")}')
        tp.writeto(fn.replace("_i2d.fits", f"_i2d_nomask-medfilt{dx}.fits"), overwrite=True)
        log.info(fn)
        log.info(fn.replace("_i2d.fits", f"_i2d_nomask-medfilt{dx}.fits"))

'''
log.info('Making Medium Filter Background without stars masked, chunked version, dx = 128 if "f4" in fn else 256')
for fn in glob.glob(f'{basepath}/images/jw02221-o002_t001_nircam_*nodestreak_i2d.fits'):
    log.info(fn)
    with warnings.catch_warnings():
        # specifically ignoring that we're using non-integer 
        warnings.simplefilter('ignore')
        dx = 128 if 'f4' in fn else 256
        sl = 'long' if 'f4' in fn else 'short'
        regs = regions.Regions.read(f'{basepath}/regions_/bright_stars_{sl}.reg')
        log.info('Computing zero spacing approxmation.')
        tp = compute_zero_spacing_approximation(fn, ext=1, dx=dx, regs=regs, percentile=10, progressbar=tqdm, smooth=False)
        #if np.nanpercentile(tp.data, 10) < 0:
        #    print(np.nanpercentile(tp.data, 10))
        #    tp.data -= np.nanpercentile(tp.data, 10)
        #tp.data[tp.data<0] = 0
        log.info(f'Writing to {fn.replace("_i2d.fits", f"_i2d_chunk-medfilt{dx}.fits")}')
        tp.writeto(fn.replace("_i2d.fits", f"_i2d_chunk-medfilt{dx}.fits"), overwrite=True)
        log.info(fn)
        log.info(fn.replace("_i2d.fits", f"_i2d_chunk-medfilt{dx}.fits"))
'''