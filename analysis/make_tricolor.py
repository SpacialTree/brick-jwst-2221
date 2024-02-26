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

available_filters = {
                     'brick': ['f410m', 'f212n', 'f466n', 'f405n', 'f187n', 'f182m', 'f444w', 'f356w', 'f200w', 'f115w'],
                     'cloudc': ['f410m', 'f212n', 'f466n', 'f405n', 'f187n', 'f182m'],
                    }

def save_rgb(img, filename, flip=-1):
    img = (img*256)
    img[img<0] = 0
    img[img>255] = 255
    img = img.astype('uint8')
    img = PIL.Image.fromarray(img[::flip,:,:])
    img.save(filename)

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

    AVM = pyavm.AVM.from_header(fits.getheader(rgb[0]))

    narrowsum_withstars = rgb_withstars[:,:,0] + rgb_withstars[:,:,1]
    rgb_scaled = np.array([
                           simple_norm(rgb_withstars[:,:,0], stretch='asinh', min_cut=-1, max_cut=90)(rgb_withstars[:,:,0]),
                           simple_norm(narrowsum_withstars,  stretch='asinh', min_cut=-2, max_cut=210)(narrowsum_withstars),
                           simple_norm(rgb_withstars[:,:,1], stretch='asinh', min_cut=-1, max_cut=120)(rgb_withstars[:,:,1]),
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
    save_rgb(rgb_scaled, outfn, flip=-1)
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

    AVM = pyavm.AVM.from_header(fits.getheader(rgb[0]))

    rgb_scaled = np.array([
                           simple_norm(rgb_withstars[:,:,0], stretch='asinh', min_cut=-1, max_cut=90)(rgb_withstars[:,:,0]),
                           simple_norm(rgb_withstars[:,:,1], stretch='asinh', min_cut=-1, max_cut=210)(rgb_withstars[:,:,1]),
                           simple_norm(rgb_withstars[:,:,2], stretch='asinh', min_cut=-1, max_cut=120)(rgb_withstars[:,:,2]),
    ]).swapaxes(0,2).swapaxes(0,1)

    plt.figure(figsize=(24,10))
    plt.imshow(rgb_scaled, origin='lower')
    plt.xticks([]);
    plt.yticks([]);

    outfn = f"{basepath}/images/{target}JWST_{module}_R-{filternames[0]}_G-{filternames[1]}_B-{filternames[2]}_rotated.png"
    save_rgb(rgb_scaled, outfn, flip=-1)
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
    module = 'merged'
    
    filternames = ['F466N', 'F405N']
    two_filter_tricolor(filternames, target, module, basepath)

    #filternames = ['F212N', 'F187N']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F405N', 'F212N']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F410M', 'F182M']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F466N', 'F410M']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F212N', 'F182M']
    #two_filter_tricolor(filternames, target, module, basepath)
    #filternames = ['F405N', 'F182M']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F405N', 'F187N']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F466N', 'F182M']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F466N', 'F187N']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F410M', 'F405N']
    #two_filter_tricolor(filternames, target, module, basepath)
    #
    #filternames = ['F466N', 'F410M']
    #two_filter_tricolor(filternames, target, module, basepath)
    #
    #filternames = ['F182M', 'F187N']
    #two_filter_tricolor(filternames, target, module, basepath)
#
    #filternames = ['F212N', 'F182M']
    #two_filter_tricolor(filternames, target, module, basepath)
    #
    #filternames = ['F466N', 'F405N', 'F187N']
    #three_filter_tricolor(filternames, target, module, basepath)
    #
    #filternames = ['F466N', 'F405N', 'F212N']
    #three_filter_tricolor(filternames, target, module, basepath)
    #    
    #filternames = ['F405N', 'F212N', 'F187N']
    #three_filter_tricolor(filternames, target, module, basepath)
    #
    #filternames = ['F466N', 'F212N', 'F187N']
    #three_filter_tricolor(filternames, target, module, basepath)


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
