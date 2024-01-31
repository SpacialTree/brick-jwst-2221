import sys
import regions
import pyavm
import numpy as np
import PIL

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import simple_norm

import reproject 
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
import matplotlib.pyplot as plt
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb

basepath = '/orange/adamginsburg/jwst/cloudc/'

files_with_stars = [
    fits.open(f'{basepath}/images/jw02221-o002_t001_nircam_clear-f405n-merged_i2d.fits'),
    fits.open(f'{basepath}/images/jw02221-o002_t001_nircam_clear-f466n-merged_i2d.fits'),
    fits.open(f'{basepath}/images/jw02221-o002_t001_nircam_clear-f410m-merged_i2d.fits'),
    fits.open(f'{basepath}/images/jw02221-o002_t001_nircam_clear-f212n-merged_i2d.fits'),
    fits.open(f'{basepath}/images/jw02221-o002_t001_nircam_clear-f182m-merged_i2d.fits'),
    fits.open(f'{basepath}/images/jw02221-o002_t001_nircam_clear-f187n-merged_i2d.fits'),
                    ]

target_wcs, target_shape = find_optimal_celestial_wcs(files_with_stars, hdu_in='SCI')

f405reprj_withstars, _ = reproject.reproject_interp(files_with_stars[0]['SCI'],
                           output_projection=target_wcs,
                           shape_out=target_shape)
(fits.PrimaryHDU(data=f405reprj_withstars,
                 header=target_wcs.to_header())
 .writeto(f'{basepath}/images/F405_reproj_merged-fortricolor.fits', overwrite=True)
)

f466reprj_withstars, _ = reproject.reproject_interp(files_with_stars[1]['SCI'],
                           output_projection=target_wcs,
                           shape_out=target_shape)
(fits.PrimaryHDU(data=f466reprj_withstars,
                 header=target_wcs.to_header())
 .writeto(f'{basepath}/images/F466_reproj_merged-fortricolor.fits', overwrite=True)
)

f410reprj_withstars, _ = reproject.reproject_interp(files_with_stars[2]['SCI'],
                           output_projection=target_wcs,
                           shape_out=target_shape)
(fits.PrimaryHDU(data=f410reprj_withstars,
                 header=target_wcs.to_header())
 .writeto(f'{basepath}/images/F410_reproj_merged-fortricolor.fits', overwrite=True)
)

f212reprj_withstars, _ = reproject.reproject_interp(files_with_stars[3]['SCI'],
                           output_projection=target_wcs,
                           shape_out=target_shape)
(fits.PrimaryHDU(data=f212reprj_withstars,
                 header=target_wcs.to_header())
 .writeto(f'{basepath}/images/F212_reproj_merged-fortricolor.fits', overwrite=True)
)

f182reprj_withstars, _ = reproject.reproject_interp(files_with_stars[4]['SCI'],
                           output_projection=target_wcs,
                           shape_out=target_shape)
(fits.PrimaryHDU(data=f182reprj_withstars,
                 header=target_wcs.to_header())
 .writeto(f'{basepath}/images/F182_reproj_merged-fortricolor.fits', overwrite=True)
)

f187reprj_withstars, _ = reproject.reproject_interp(files_with_stars[5]['SCI'],
                           output_projection=target_wcs,
                           shape_out=target_shape)
(fits.PrimaryHDU(data=f187reprj_withstars,
                 header=target_wcs.to_header())
 .writeto(f'{basepath}/images/F187_reproj_merged-fortricolor.fits', overwrite=True)
)