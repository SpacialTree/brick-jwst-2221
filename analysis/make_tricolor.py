import sys
import regions
import pyavm
import numpy as np
import PIL
import shutil

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D

import reproject 
from reproject.mosaicking import find_optimal_celestial_wcs, reproject_and_coadd
import matplotlib.pyplot as plt
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb

import cv2

available_filters = {
                     'brick': ['f410m', 'f212n', 'f466n', 'f405n', 'f187n', 'f182m', 'f444w', 'f356w', 'f200w', 'f115w'],
                     'cloudc': ['f410m', 'f212n', 'f466n', 'f405n', 'f187n', 'f182m'],
                    }

def save_rgb(img, filename, flip=-1, avm=None):
    #mask = (~np.isnan(img[:,:,0]))
    #mask_0 = np.logical_and(img[:,:,0]==0, img[:,:,1]==0, img[:,:,2]==0)

    #mask_nan = np.logical_and(np.isnan(img[:,:,0]), np.isnan(img[:,:,1]), np.isnan(img[:,:,2]))
    #mask = mask_nan #np.logical_or(mask_0, mask_nan)
    #kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3,3))
    #dilate = cv2.dilate(mask.astype('uint8'), kernel, iterations=1)
    #mask = np.array(dilate, dtype=bool)
    #
    #img = np.append(img, np.array([~mask]).swapaxes(0,2).swapaxes(0,1), axis=2) # 
    #img = (img*256)
    #img[img<0] = 0
    #img[img>255] = 255
    #img = img.astype('uint8')

    img = (img*256)
    img[img<0] = 0
    img[img>255] = 255

    mask_0 = np.logical_and(img[:,:,0]==0, img[:,:,1]==0, img[:,:,2]==0)
    mask_nan = np.logical_and(np.isnan(img[:,:,0]), np.isnan(img[:,:,1]), np.isnan(img[:,:,2]))
    mask = np.logical_or(mask_0, mask_nan)
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3,3))
    dilate = cv2.dilate(mask.astype('uint8'), kernel, iterations=1)
    mask = np.array(dilate, dtype=bool)

    img = np.append(img, np.array([~mask]).swapaxes(0,2).swapaxes(0,1)*255, axis=2)
    img = img.astype('uint8')

    #mask_arr = mask_arr.astype('uint8')
    #mask = PIL.Image.fromarray(mask_arr[::flip,:], mode='L')

    # saving
    img = PIL.Image.fromarray(img[::flip,:,:], mode='RGBA')
    img.save(filename)

    if avm is not None:
        avm.embed(filename, filename.replace('.png', '_avm.png'))
        shutil.move(filename.replace('.png', '_avm.png'), filename)
    
"""
def save_rgb(img, filename, flip=-1, avm=None):
    img = (img*256)
    img[img<0] = 0
    img[img>255] = 255
    img = img.astype('uint8')
    img = PIL.Image.fromarray(img[::-1,:,:])
    img.save(filename)
    if avm is not None:
        avm.embed(filename, 'avm_'+filename)
        shutil.move('avm_'+filename, filename)
"""

def two_filter_tricolor(filternames, target, module, basepath):
    rgb = [
        f'{basepath}/images/{filternames[0][:-1].upper()}_reproj_{module}-fortricolor.fits',
        f'{basepath}/images/{filternames[1][:-1].upper()}_reproj_{module}-fortricolor.fits',
        ]
    
    rgb_withstars = np.array(
      [
       fits.getdata(rgb[0]),
       fits.getdata(rgb[1]),
      ]
    ).swapaxes(0,2).swapaxes(0,1)

    rgb_withstars[rgb_withstars==0]=np.nan

    AVM = pyavm.AVM.from_header(fits.getheader(rgb[0]))

    narrowsum_withstars = rgb_withstars[:,:,0] + rgb_withstars[:,:,1]
    rgb_scaled = np.array([
                           simple_norm(rgb_withstars[:,:,0], stretch='asinh', min_cut=-1, max_cut=90)(rgb_withstars[:,:,0]),
                           simple_norm(narrowsum_withstars,  stretch='asinh', min_cut=-2, max_cut=210)(narrowsum_withstars),
                           simple_norm(rgb_withstars[:,:,1], stretch='asinh', min_cut=-1, max_cut=120)(rgb_withstars[:,:,1]),
                           #np.full(narrowsum_withstars.shape, 0), # transparency layer 
    ]).swapaxes(0,2).swapaxes(0,1)
    #hsv = rgb_to_hsv(rgb_scaled)
    #hsv[:,:,0] += -0.35  # 0.25 = 90/360
    #hsv[:,:,0] = hsv[:,:,0] % 1 
    #rgb_scaled = hsv_to_rgb(hsv)
    plt.figure(figsize=(24,10))
    plt.imshow(rgb_scaled, origin='lower')
    plt.xticks([]);
    plt.yticks([]);

    outfn = f"{basepath}/images/{target}JWST_{module}_R-{filternames[0]}_B-{filternames[1]}_rotated.png"
    save_rgb(rgb_scaled, outfn, flip=-1, avm=AVM)
    AVM.embed(outfn, outfn)
    

def three_filter_tricolor(filternames, target, module, basepath):
    rgb = [
        f'{basepath}/images/{filternames[0][:-1].upper()}_reproj_{module}-fortricolor.fits',
        f'{basepath}/images/{filternames[1][:-1].upper()}_reproj_{module}-fortricolor.fits',
        f'{basepath}/images/{filternames[2][:-1].upper()}_reproj_{module}-fortricolor.fits',
            ]
    rgb_withstars = np.array(
        [
            fits.getdata(rgb[0]),
            fits.getdata(rgb[1]),
            fits.getdata(rgb[2]),
        ]
    ).swapaxes(0,2).swapaxes(0,1)

    rgb_withstars[rgb_withstars==0]=np.nan

    AVM = pyavm.AVM.from_header(fits.getheader(rgb[0]))

    rgb_scaled = np.array([
                           simple_norm(rgb_withstars[:,:,0], stretch='asinh', min_cut=-1, max_cut=90)(rgb_withstars[:,:,0]),
                           simple_norm(rgb_withstars[:,:,1], stretch='asinh', min_cut=-1, max_cut=210)(rgb_withstars[:,:,1]),
                           simple_norm(rgb_withstars[:,:,2], stretch='asinh', min_cut=-1, max_cut=120)(rgb_withstars[:,:,2]),
                           #np.full(rgb_withstars[:,:,2].shape, 0), # transparency layer 
    ]).swapaxes(0,2).swapaxes(0,1)

    plt.figure(figsize=(24,10))
    plt.imshow(rgb_scaled, origin='lower')
    plt.xticks([]);
    plt.yticks([]);

    outfn = f"{basepath}/images/{target}JWST_{module}_R-{filternames[0]}_G-{filternames[1]}_B-{filternames[2]}_rotated.png"
    save_rgb(rgb_scaled, outfn, flip=-1, avm=AVM)
    AVM.embed(outfn, outfn)

def get_cutout(filename, position, l, w, format='fits'):
    if format == 'fits':
        try: 
            hdu = fits.open(filename, ext='SCI')[0]
        except: 
            hdu = fits.open(filename)[0]
    elif format == 'casa':
        hdu = SpectralCube.read(filename, format='casa').hdu
    data = np.squeeze(hdu.data)
    head = hdu.header

    #pixel_scale = head['PIXSCALE']*u.arcsec/u.pix
    ww = WCS(head).celestial
    size = (l, w)
    #((l/pixel_scale).to(u.pix), (w/pixel_scale).to(u.pix))
    cutout = Cutout2D(data, position=position, size=size, wcs=ww)
    return cutout

def spitzer_tricolor(position, l, w, target, basepath):
    fn_I2 = '/orange/adamginsburg/cmz/glimpse_data/GLM_00000+0000_mosaic_I2.fits'
    fn_I3 = '/orange/adamginsburg/cmz/glimpse_data/GLM_00000+0000_mosaic_I3.fits'
    fn_I4 = '/orange/adamginsburg/cmz/glimpse_data/GLM_00000+0000_mosaic_I4.fits'
    
    cutout_I2 = get_cutout(fn_I2, position, l, w)
    cutout_I3 = get_cutout(fn_I3, position, l, w)
    cutout_I4 = get_cutout(fn_I4, position, l, w)
        
    rgb = np.array([cutout_I2.data,cutout_I3.data,cutout_I4.data]).swapaxes(0,2).swapaxes(0,1)
    rgb[rgb==0]=np.nan
    
    AVM = pyavm.AVM.from_header(fits.getheader(fn_I2))
    
    rgb_scaled = np.array([
                           simple_norm(rgb[:,:,-1], stretch='linear', min_cut=10, max_cut=350)(rgb[:,:,-1]), # 350
                           simple_norm(rgb[:,:,-2], stretch='linear', min_cut=10, max_cut=200)(rgb[:,:,-2]), # 200
                           simple_norm(rgb[:,:,-3], stretch='linear', min_cut=10, max_cut=100)(rgb[:,:,-3]), # 100
                           #np.full(rgb[:,:,2].shape, 0), # transparency layer 
    ]).swapaxes(0,2).swapaxes(0,1)
    #return rgb_scaled, cutout_I4.wcs
    plt.figure(figsize=(24,10))
    plt.imshow(rgb_scaled, origin='lower')
    plt.xticks([]);
    plt.yticks([]);

    outfn = f"{basepath}/images/{target}Spitzer_R-I4_G-I3_B-I2.png"
    save_rgb(rgb_scaled, outfn, flip=-1, avm=AVM)
    AVM.embed(outfn, outfn)

# Run reproj_fortricolor first to have the files available
def main():
    #fn_182 = '/orange/adamginsburg/jwst/cloudc/images/F182_reproj_merged-fortricolor.fits'
    #fn_187 = '/orange/adamginsburg/jwst/cloudc/images/F187_reproj_merged-fortricolor.fits'
    #fn_212 = '/orange/adamginsburg/jwst/cloudc/images/F212_reproj_merged-fortricolor.fits'
    #fn_405 = '/orange/adamginsburg/jwst/cloudc/images/F405_reproj_merged-fortricolor.fits'
    #fn_410 = '/orange/adamginsburg/jwst/cloudc/images/F410_reproj_merged-fortricolor.fits'
    #fn_466 = '/orange/adamginsburg/jwst/cloudc/images/F466_reproj_merged-fortricolor.fits'

    target = 'cloudc'
    basepath = f'/orange/adamginsburg/jwst/{target}/'
    module = 'merged-reproject'
    
    filternames = ['F466N', 'F405N']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F212N', 'F187N']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F405N', 'F212N']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F410M', 'F182M']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F466N', 'F410M']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F212N', 'F182M']
    two_filter_tricolor(filternames, target, module, basepath)
    filternames = ['F405N', 'F182M']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F405N', 'F187N']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F466N', 'F182M']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F466N', 'F187N']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F410M', 'F405N']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F466N', 'F410M']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F182M', 'F187N']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F212N', 'F182M']
    two_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F466N', 'F405N', 'F187N']
    three_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F466N', 'F405N', 'F212N']
    three_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F405N', 'F212N', 'F187N']
    three_filter_tricolor(filternames, target, module, basepath)
    
    filternames = ['F466N', 'F212N', 'F187N']
    three_filter_tricolor(filternames, target, module, basepath)
#
    filternames = ['F466N', 'F410M', 'F405N']
    three_filter_tricolor(filternames, target, module, basepath)
#
    filternames = ['F212N', 'F182M', 'F187N']
    three_filter_tricolor(filternames, target, module, basepath)
#
    position = SkyCoord('17:46:21.4701708277', '-28:35:38.0673181068', unit=(u.hourangle, u.deg))
    l = 7*u.arcmin
    w = 10*u.arcmin
    spitzer_tricolor(position, l, w, target, basepath)

    '''from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-f", "--filternames", dest="filternames",
                      default='F466N,F405N',
                      help="filter name list", metavar="filternames")
    parser.add_option("-m", "--module", dest="module",
                    default='merged',
                    help="module", metavar="module")
    parser.add_option("-t", "--target", dest="target",
                    default='cloudc',
                    help="target field", metavar="target")
    (options, args) = parser.parse_args()

    filternames = options.filternames.split(",")
    module = options.module
    target = options.target

    basepath = f'/orange/adamginsburg/jwst/{target}/'
    
    for filt in filternames:
        if filt not in available_filters[target]:
            print(f"WARNING: filter {filt} is not in list of available filters {available_filters[target]}.")
            return False

    if len(filternames) > 3 or len(filternames) < 2:
        print(f'WARNING: Length of filternames list is {len(filternames)}, that is too few/many to make an rgb image out of!')
        return False

    elif len(filternames) == 2:
        two_filter_tricolor(filternames, target, module, basepath)

    elif len(filternames) == 3:
        three_filter_tricolor(filternames, target, module, basepath)
    '''
    

if __name__ == "__main__":
    main()
