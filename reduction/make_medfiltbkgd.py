from scipy import ndimage
import scipy.ndimage
from astropy.io import fits
import pylab as pl
from astropy import convolution
from astropy.convolution import Gaussian2DKernel
from astropy import units as u
import pyavm
import regions
from astropy import wcs
import glob
import numpy as np
import warnings
import regions
from tqdm.notebook import tqdm
from astropy import log
import sys
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

def make_percentile_background(filtername, basepath='/orange/adamginsburg/jwst/cloudc/'):
    files = glob.glob(f'{basepath}/images/jw02221-o002_t001_nircam_*{filtername.lower()}*nodestreak_i2d.fits')
    for fn in files: 
        log.info(f'Making Percentile Filter Background for {fn}')
        with warnings.catch_warnings():
            # specifically ignoring that we're using non-integer 
            warnings.simplefilter('ignore')
            dx = 64 if 'f4' in fn else 128
            log.info('Computing zero spacing approxmation.')
            tp = compute_zero_spacing_approximation(fn, ext=1, dx=dx, percentile=10, progressbar=tqdm)
            log.info(f'Writing to {fn.replace("_i2d.fits", f"_i2d-perfilt{dx}.fits")}')
            tp.writeto(fn.replace("_i2d.fits", f"_i2d-perfilt{dx}.fits"), overwrite=True)
            log.info(fn)
            log.info(fn.replace("_i2d.fits", f"_i2d-perfilt{dx}.fits"))

def make_median_background(filtername, basepath='/orange/adamginsburg/jwst/cloudc/'):
    files = glob.glob(f'{basepath}/images/jw02221-o002_t001_nircam_*{filtername.lower()}*nodestreak_i2d.fits')
    for fn in files: 
        log.info(f'Making Median Filter Background for {fn}')
        with warnings.catch_warnings():
            # specifically ignoring that we're using non-integer 
            warnings.simplefilter('ignore')
            dx = 128 if 'f4' in fn else 256
            sl = 'long' if 'f4' in fn else 'short'
            regs = regions.Regions.read(f'{basepath}/regions_/bright_stars_{sl}.reg')
            log.info('Computing zero spacing approxmation.')
            tp = compute_zero_spacing_approximation(fn, ext=1, dx=dx, percentile=10, progressbar=tqdm, smooth=False)
            log.info(f'Writing to {fn.replace("_i2d.fits", f"_i2d-medfilt{dx}.fits")}')
            tp.writeto(fn.replace("_i2d.fits", f"_i2d-medfilt{dx}.fits"), overwrite=True)
            log.info(fn)
            log.info(fn.replace("_i2d.fits", f"_i2d-medfilt{dx}.fits"))


def main(filtername, basepath='/orange/adamginsburg/jwst/cloudc/', chunky=False):
    if not chunky: 
        make_percentile_background(filtername=filtername, basepath=basepath)
    else:
        make_median_background(filtername=filtername, basepath=basepath)
        return False

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--filternames", dest="filternames",
                      default='F466N,F405N,F410M,F212N,F182M,F187N',
                      help="filter name list", metavar="filternames")
    parser.add_option("-t", "--target", dest="target",
                    default='cloudc',
                    help="target field", metavar="target")
    parser.add_option("--chunky", dest="chunky",
                      default=False,
                      action='store_true',
                      help="Use alternative method?", metavar="chunky")
    (options, args) = parser.parse_args()

    filternames = options.filternames.split(",")
    target = options.target
    chunky = bool(options.chunky)
    
    print(filternames)
    print(target)
    print(chunky)

    if target == 'cloudc' or target == 'brick':
        basepath = f'/orange/adamginsburg/jwst/{target}/'
        for filt in filternames:
            main(filtername=filt, basepath=basepath, chunky=chunky)
    else: 
        print(f'Target {target} is not available.')
