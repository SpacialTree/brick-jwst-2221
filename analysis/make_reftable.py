import datetime

from astropy import units as u
from astropy.table import Table

from astropy.coordinates import SkyCoord
from astroquery.svo_fps import SvoFps
from astropy import wcs
from astropy.io import fits
import numpy as np

from astropy import stats


def main(basepath = '/orange/adamginsburg/jwst/cloudc/'):
    """
    June 28, 2023: decided to switch to F405N-only reference
    """
    #basepath = '/blue/adamginsburg/adamginsburg/jwst/brick/'

    tblfilename = (f'{basepath}/catalogs/f405n_merged_indivexp_merged_crowdsource_nsky0.fits')
    tbl = Table.read(tblfilename)

    try:
        sel = ((tbl['qf'] > 0.95) & (tbl['spread_model'] < 0.25) & (tbl['fracflux'] > 0.9) & (tbl['flux'] > 0))
    except KeyError:
        sel = ((tbl['qf'] > 0.95) & (tbl['fracflux'] > 0.9) & (tbl['flux'] > 0))

    print(f"QFs are good for {sel.sum()} out of {len(tbl)} catalog entries")
    print(f"Making the reference catalog from {sel.sum()} out of {len(tbl)} catalog entries")

    reftbl = tbl['skycoord', 'flux' ][sel]

    # Crossmatch to VVV and recenter
    vvvtb = Table.read(f'{basepath}/F405N/pipeline/jw02221-o001_t001_nircam_clear-f405n-merged_vvvcat.ecsv')
    vvvcrds = SkyCoord(vvvtb['RAJ2000'].quantity, vvvtb['DEJ2000'].quantity, frame='fk5')
    refcrds = reftbl['skycoord']

    jfilts = SvoFps.get_filter_list('JWST')
    jfilts.add_index('filterID')
    zeropoint = u.Quantity(jfilts.loc[f'JWST/NIRCam.F405N']['ZeroPoint'], u.Jy)
    #sqpixscale = wcs.WCS(fits.getheader(reftbl.meta['FILENAME'], ext=1)).proj_plane_pixel_area()
    flux_jy = (reftbl['flux'] * u.MJy/u.sr * (reftbl.meta['pixscale_as']*u.arcsec)**2).to(u.Jy)
    mag405 = -2.5 * np.log10(flux_jy / zeropoint) * u.mag

    total_dra, total_ddec = 0*u.arcsec, 0*u.arcsec

    for ii in range(8):
        keep = np.zeros(len(reftbl), dtype='bool')
        idx, sidx, sep, _ = refcrds.search_around_sky(vvvcrds, 0.2*u.arcsec)

        is_closest = np.array([(ii == idx).sum() == 1 or (sp == sep[idx == ii].min()) for ii, sp in zip(idx, sep)])
        idx = idx[is_closest]
        sidx = sidx[is_closest]

        # magnitude difference = ratio
        ratio = vvvtb['Ksmag3'][idx] - mag405[sidx]
        reject = np.zeros(ratio.size, dtype='bool')
        for jj in range(12):
            madstd = stats.mad_std(ratio[~reject])
            med = np.median(ratio[~reject])
            reject = (ratio < med - 3 * madstd) | (ratio > med + 3 * madstd) | reject
            ratio = 1 / ratio
            madstd = stats.mad_std(ratio[~reject])
            med = np.median(ratio[~reject])
            reject = (ratio < med - 3 * madstd) | (ratio > med + 3 * madstd) | reject
            ratio = 1 / ratio

        keep[sidx[~reject]] = True
        dra = (refcrds[sidx].ra - vvvcrds[idx].ra).to(u.arcsec)
        ddec = (refcrds[sidx].dec - vvvcrds[idx].dec).to(u.arcsec)
        dra_med, ddec_med = np.median(dra[~reject]), np.median(ddec[~reject])
        print(f"{len(idx)} VVV matches in f405n catalog.  ", end='')
        print(f'ratio = {med} +/- {madstd}   nkeep={(~reject).sum()}.  dra, ddec median = {dra_med}, {ddec_med}')
        total_dra += dra_med
        total_ddec += ddec_med

        refcrds_updated = SkyCoord(refcrds.ra - dra_med, refcrds.dec - ddec_med, frame=refcrds.frame)
        refcrds = refcrds_updated

        reftbl['VVV_matched'] = keep

        if (np.abs(dra_med) < 1*u.marcsec) & (np.abs(ddec_med) < 1*u.marcsec):
            break

    print(f"Shifted F405N coordinates by {total_dra}, {total_ddec} in {ii} iterations with stddev = {dra.std()}, {ddec.std()} ({(dra.var()+ddec.var())**0.5})")
    reftbl['skycoord'] = refcrds_updated

    # include two columns to make it a table, plus abmag for sorting
    reftbl['RA'] = reftbl['skycoord'].ra
    reftbl['DEC'] = reftbl['skycoord'].dec
    reftbl.sort('flux', reverse=True) # descending

    reftbl.meta['VERSION'] = datetime.datetime.now().isoformat()
    if 'VERSION' in tbl.meta:
        reftbl.meta['PARENT_VERSION'] = tbl.meta['VERSION']
        reftbl.meta['RAOFFSET'] = total_dra
        reftbl.meta['DECOFFSET'] = total_ddec

    reftbl.write(f'{basepath}/catalogs/crowdsource_based_nircam-f405n_reference_astrometric_catalog.ecsv', overwrite=True)
    reftbl.write(f'{basepath}/catalogs/crowdsource_based_nircam-f405n_reference_astrometric_catalog.fits', overwrite=True)

    return reftbl


def main_old():
    basepath = '/blue/adamginsburg/adamginsburg/jwst/cloudc/'
    long_filternames = ['f410m', 'f405n', 'f466n']

    # filtername = 'F410M'
    # module = 'merged'
    # tblfilename = f"{basepath}/{filtername}/{filtername.lower()}_{module}_crowdsource_nsky0.fits"

    # May 19, 2023: changed this to 'merged' b/c we can't keep going on with half a field; the workflow
    # relies on having a common catalog for both!
    # June 24, 2023: changed to merged-reproject, which worked, while merged did not.
    tblfilename = (f'{basepath}/catalogs/crowdsource_nsky0_merged-reproject_photometry_tables_merged.fits')
    tbl = Table.read(tblfilename)

    # reject sources with nondetections in F405N or F466N or bad matches
    sel = tbl['sep_f466n'].quantity < 0.1*u.arcsec
    sel &= tbl['sep_f405n'].quantity < 0.1*u.arcsec

    # reject sources with bad QFs
    goodqflong = ((tbl['qf_f410m'] > 0.90) |
                  (tbl['qf_f405n'] > 0.90) |
                  (tbl['qf_f466n'] > 0.90))
    goodspreadlong = ((tbl['spread_model_f410m'] < 0.25) |
                      (tbl['spread_model_f405n'] < 0.25) |
                      (tbl['spread_model_f466n'] < 0.25))
    goodfracfluxlong = ((tbl['fracflux_f410m'] > 0.8) |
                        (tbl['fracflux_f405n'] > 0.8) &
                        (tbl['fracflux_f466n'] > 0.8))

    any_saturated_ = [(tbl[f'near_saturated_{x}_{x}'] &
                      ~tbl[f'flux_{x}'].mask)
                      for x in long_filternames]
    any_saturated = any_saturated_[0]
    for col in any_saturated_[1:]:
        any_saturated = any_saturated | col

    any_replaced_saturated_ = [tbl[f'replaced_saturated_{x}'] &
                               ~tbl[f'flux_{x}'].mask for x in long_filternames]
    any_replaced_saturated = any_replaced_saturated_[0]
    for col in any_replaced_saturated_[1:]:
        any_replaced_saturated = any_replaced_saturated | col


    sel &= goodqflong & goodspreadlong & goodfracfluxlong
    print(f"QFs are good for {sel.sum()} out of {len(tbl)} catalog entries")
    print(f"Rejecting {(any_replaced_saturated & sel).sum()} replaced-saturated sources.")
    sel &= ~any_replaced_saturated
    print(f"Rejecting {(any_saturated & sel).sum()} near-saturated sources.")
    sel &= ~any_saturated
    print(f"Making the reference catalog from {sel.sum()} out of {len(tbl)} catalog entries")

    # include two columns to make it a table, plus abmag for sorting
    reftbl = tbl['skycoord_f410m', 'skycoord_f405n', 'mag_ab_f410m' ][sel]
    reftbl['RA'] = reftbl['skycoord_f410m'].ra
    reftbl['DEC'] = reftbl['skycoord_f410m'].dec
    reftbl.sort('mag_ab_f410m')

    reftbl.meta['VERSION'] = datetime.datetime.now().isoformat()
    reftbl.meta['PARENT_VERSION'] = tbl.meta['VERSION']

    reftbl.write(f'{basepath}/catalogs/crowdsource_based_nircam-long_reference_astrometric_catalog.ecsv', overwrite=True)
    reftbl.write(f'{basepath}/catalogs/crowdsource_based_nircam-long_reference_astrometric_catalog.fits', overwrite=True)

    return reftbl

if __name__ == "__main__":
    reftbl = main()
