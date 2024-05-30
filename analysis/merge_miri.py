import numpy as np
import time
import datetime
import os
import warnings
import regions
from astropy.io import fits
import glob
from photutils.background import MMMBackground, MADStdBackgroundRMS
from photutils.aperture import CircularAperture, CircularAnnulus
from photutils.detection import DAOStarFinder, IRAFStarFinder, find_peaks
from photutils.psf import (DAOGroup, IntegratedGaussianPRF, extract_stars,
                           IterativelySubtractedPSFPhotometry,
                           BasicPSFPhotometry, EPSFBuilder)
from astropy.modeling.fitting import LevMarLSQFitter
from astropy import stats
from astropy.table import Table, Column, MaskedColumn
from astropy.wcs import WCS
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import coordinates
from astropy.visualization import simple_norm
from astropy import wcs
from astropy import table
from astropy import units as u
from astroquery.svo_fps import SvoFps
import pylab as pl
pl.rcParams['figure.facecolor'] = 'w'
pl.rcParams['image.origin'] = 'lower'
pl.rcParams['figure.figsize'] = (10,8)
pl.rcParams['figure.dpi'] = 100


def getmtime(x):
    return datetime.datetime.fromtimestamp(os.path.getmtime(x)).strftime('%Y-%m-%d %H:%M:%S')

def merge_catalogs(basetable, basepath='/orange/adamginsburg/jwst/cloudc/'):
    print("Starting to merge catalogs")
    ref_filter = 'f2550w'
    max_offset = 0.15*u.arcsec

    fn_nir = f'{basepath}/catalogs/crowdsource_nsky1_merged_photometry_tables_merged.ecsv'
    nirtable = Table.read(fn_nir)

    jfilts = SvoFps.get_filter_list('JWST')
    jfilts.add_index('filterID')
    
    reffiltercol = [ref_filter] * len(basetable)
    
    nirtable_mask = np.logical_and(nirtable['qf_f405n']>0.6, nirtable['qf_f410m']>0.6)
    nirtable_use = tbl = nirtable[nirtable_mask]
    
    basecrds = basetable['skycoord']
    crds = tbl['skycoord_ref']
    
    meta = {}
    
    t0 = time.time()
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        wl = ref_filter
        matches, sep, _ = basecrds.match_to_catalog_sky(crds, nthneighbor=1)
        basetable.add_column(name=f"sep_{wl}", col=sep)
        basetable.add_column(name=f"id_{wl}", col=matches)
        matchtb = tbl[matches]
        badsep = sep > max_offset
        for cn in basetable.colnames:
            if isinstance(basetable[cn], SkyCoord):
                basetable.rename_column(cn, f"{cn}_{wl}")
                basetable[f'mask_{wl}'] = badsep
            else:
                basetable[f'{cn}_{wl}'] = MaskedColumn(data=basetable[cn], name=f'{cn}_{wl}')
                #basetable[f'{cn}_{wl}'].mask[badsep] = True
                if hasattr(basetable[cn], 'meta'):
                    basetable[f'{cn}_{wl}'].meta = basetable[cn].meta
                basetable.remove_column(cn)

    meta[f'{wl[1:-1]}pxdg'.upper()] = basetable.meta['pixelscale_deg2']
    meta[f'{wl[1:-1]}pxas'.upper()] = basetable.meta['pixelscale_arcsec']
    basetable = table.hstack([matchtb, basetable], join_type='exact')
    for key in basetable.meta:
        meta[f'{wl[1:-1]}{key[:4]}'.upper()] = basetable.meta[key]
    for key in nirtable.meta:
        meta[f'{wl[1:-1]}{key[:4]}'.upper()] = nirtable.meta[key]
    basetable.meta = meta
    #assert '212PXDG' in meta
    #assert '212PXDG' in basetable.meta
    #assert '2550PXDG' in meta
    #assert '2550PXDG' in basetable.meta
    
    tablename = f"{basepath}/catalogs/miri_combined_crowdsource_nsky1_photometry_tables_merged"
    t0 = time.time()
    print(f"Writing table {tablename}")
    # use caps b/c FITS will force it to caps anyway
    basetable.meta['VERSION'] = datetime.datetime.now().isoformat()
    # takes FOR-EV-ER
    basetable.write(f"{tablename}.ecsv", overwrite=True)
    print(f"Done writing table {tablename}.ecsv in {time.time()-t0:0.1f} seconds")
    t0 = time.time()
    # DO NOT USE FITS in production, it drops critical metadata
    # I wish I had noted *what* metadata it drops, though, since I still seem to be using
    # it in production code down the line...
    # OH, I think the FITS file turns "True" into "False"?
    # Yes, specifically: it DROPS masked data types, converting "masked" into "True"?
    for colname in basetable.colnames:    # Adding to try to fix FITS file turning "masked" into True
        try:  
            if basetable[colname].dtype.name == 'bool':
                basetable[colname][basetable[colname].mask] = False
        except:
            print('no dtype')
    basetable.write(f"{tablename}.fits", overwrite=True)
    print(f"Done writing table {tablename}.fits in {time.time()-t0:0.1f} seconds")
    

def merge_miri(basepath='/orange/adamginsburg/jwst/cloudc/'):
    print("Starting MIRI photometry calculations")

    fn_mir = f'{basepath}/F2550W/f2550w__crowdsource_nsky1_real.fits'

    img_fn = '/orange/adamginsburg/jwst/cloudc/images/jw02221-o001_t001_miri_f2550w_i2d.fits'
    img = fits.getdata(img_fn, ext=('SCI', 1))
    ww = wcs.WCS(fits.getheader(img_fn, ext=('SCI', 1)))

    jfilts = SvoFps.get_filter_list('JWST')
    jfilts.add_index('filterID')

    basetable = tbl = Table.read(fn_mir)

    tbl.meta['filename'] = fn_mir
    tbl.meta['filter'] = os.path.basename(fn_mir).split("_")[0]
    
    if 'skycoord' not in tbl.colnames:
        crds = ww.pixel_to_world(tbl['x'], tbl['y'])
        tbl.add_column(crds, name='skycoord')
    else:
        crds = tbl['skycoord']
    
    tbl.meta['pixelscale_deg2'] = ww.proj_plane_pixel_area()
    tbl.meta['pixelscale_arcsec'] = (ww.proj_plane_pixel_area()**0.5).to(u.arcsec)
    print('Calculating Flux [Jy]')
    flux_jy = (tbl['flux'] * u.MJy/u.sr * (2*np.pi / (8*np.log(2))) * tbl['fwhm']**2 * tbl.meta['pixelscale_deg2']).to(u.Jy)
    eflux_jy = (tbl['dflux'] * u.MJy/u.sr * (2*np.pi / (8*np.log(2))) * tbl['fwhm']**2 * tbl.meta['pixelscale_deg2']).to(u.Jy)
    with np.errstate(all='ignore'):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            filtername = tbl.meta["filter"]
            zeropoint = u.Quantity(jfilts.loc[f'JWST/MIRI.{filtername.upper()}']['ZeroPoint'], u.Jy)
            abmag = -2.5 * np.log10(flux_jy / zeropoint)
            abmag_err = 2.5 / np.log(10) * np.abs(eflux_jy / flux_jy)
            tbl.add_column(flux_jy, name='flux_jy')
            tbl.add_column(eflux_jy, name='eflux_jy')
            tbl.add_column(abmag, name='mag_ab')
            tbl.add_column(abmag_err, name='emag_ab')
    if hasattr(tbl['mag_ab'], 'mask'):
        print(f'ab mag tbl col has mask sum = {tbl["mag_ab"].mask.sum()} masked values')
    if hasattr(abmag, 'mask'):
        print(f'ab mag has mask sum = {abmag.mask.sum()} masked values')
    if hasattr(tbl['flux'], 'mask'):
        print(f'ab mag has mask sum = {tbl["flux"].mask.sum()} masked values')

    merge_catalogs(tbl, basepath)
    
def real_miri(basepath):
    # Make a mask that confirms all miri sources that are not ISM

    cat = Table.read(f'{basepath}/F2550W/f2550w__crowdsource_nsky1.fits')
    crds = cat['skycoord']
    manual_reg = regions.Regions.read(f'{basepath}/regions_/miri_catalog_new.reg', format='ds9')
    ra = []
    dec = []
    for reg in manual_reg: 
        ra.append(reg.center.ra)
        dec.append(reg.center.dec)
    manual_catalog = SkyCoord(ra=ra, dec=dec)
    
    matches, sep, _ = crds.match_to_catalog_sky(manual_catalog, nthneighbor=1)
    mask = sep < 0.33*u.arcsec
    masked_crds = crds[mask]
    
    cat['real_miri'] = mask
    
    cat.write(f'{basepath}/F2550W/f2550w__crowdsource_nsky1_real.fits', format='fits', overwrite=True)

def main():
    print('Starting main')
    
    target = 'cloudc'
    basepath = f'/orange/adamginsburg/jwst/{target}/'
    
    real_miri(basepath)
    merge_miri(basepath)
    
if __name__ == "__main__":
    main()
