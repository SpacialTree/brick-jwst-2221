from astropy.table import Table
import numpy as np

basepath = '/orange/adamginsburg/jwst/cloudc/'

basetable = Table.read(f'{basepath}/catalogs/crowdsource_nsky0_merged_photometry_tables_merged.ecsv')

print('Basetable', len(basetable))

short_near = np.logical_or(basetable['near_saturated_f182m_f182m'], basetable['near_saturated_f187n_f187n'], basetable['near_saturated_f212n_f212n'])
longg_near = np.logical_or(basetable['near_saturated_f410m_f410m'], basetable['near_saturated_f405n_f405n'], basetable['near_saturated_f466n_f466n'])
sat_neat = np.logical_or(short_near, longg_near)

print('Near Saturated Stars', np.sum(sat_neat))

print('F182M', np.sum(basetable['near_saturated_f182m_f182m']))
print('F187N', np.sum(basetable['near_saturated_f187n_f187n']))
print('F212N', np.sum(basetable['near_saturated_f212n_f212n']))
print('F410M', np.sum(basetable['near_saturated_f410m_f410m']))
print('F405N', np.sum(basetable['near_saturated_f405n_f405n']))
print('F466N', np.sum(basetable['near_saturated_f466n_f466n']))