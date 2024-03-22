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

print('Starting')

vlow = -65*u.km/u.second
vhigh = -45*u.km/u.second
print('Making cubes and moments from', vlow, 'to', vhigh)

reg_filament = regions.Regions.read('/orange/adamginsburg/jwst/cloudc/regions_/filament_reg.reg')

fn_HNCO_mosaic = '/orange/adamginsburg/ACES/mosaics/cubes/HNCO_7m12mTP_CubeMosaic.fits'

fn_s = ['/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a8/calibrated/working/uid___A001_X15a0_X1a8.s38_0.Sgr_A_star_sci.spw25.cube.I.iter1.image.pbcor.fits',
        '/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a8/calibrated/working/uid___A001_X15a0_X1a8.s38_0.Sgr_A_star_sci.spw29.cube.I.iter1.image.pbcor.fits',
        '/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a8/calibrated/working/uid___A001_X15a0_X1a8.s38_0.Sgr_A_star_sci.spw27.cube.I.iter1.image.pbcor.fits',
        '/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a8/calibrated/working/uid___A001_X15a0_X1a8.s38_0.Sgr_A_star_sci.spw35.cube.I.iter1.image.pbcor.fits',
        '/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a8/calibrated/working/uid___A001_X15a0_X1a8.s38_0.Sgr_A_star_sci.spw33.cube.I.iter1.image.pbcor.fits',
        '/orange/adamginsburg/ACES/rawdata/2021.1.00172.L/science_goal.uid___A001_X1590_X30a8/group.uid___A001_X1590_X30a9/member.uid___A001_X15a0_X1a8/calibrated/working/uid___A001_X15a0_X1a8.s38_0.Sgr_A_star_sci.spw31.cube.I.iter1.image.pbcor.fits']

# HNCO mosaic
print('HNCO')
cube_full = SpectralCube.read(fn_HNCO_mosaic).with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=87.925238*u.GHz)

cube_full_filament = cube_full.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_full_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_HNCO_mos_7m12mTP_filament.fits', overwrite=True)

mom0_HNCO_filament = cube_full_filament.moment0()
mom0_HNCO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_HNCO_mos_7m12mTP_filament.fits', overwrite=True)

# SiO 2-1 mosaic
print('SiO')
fn_SiO_mosaic = '/orange/adamginsburg/ACES/mosaics/cubes/SiO21_CubeMosaic.fits'
cube_SiO_mosaic = SpectralCube.read(fn_SiO_mosaic)

cube_SiO = cube_SiO_mosaic.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=86.84696*u.GHz)

cube_SiO_filament = cube_SiO.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_SiO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_SiO_mos_filament.fits', overwrite=True)

mom0_SiO_filament = cube_SiO_filament.moment0()
mom0_SiO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_SiO_mos_filament.fits', overwrite=True)

# CS 2-1 mosaic
print('CS')
fn_CS_mosaic = '/orange/adamginsburg/ACES/mosaics/cubes/CS21_CubeMosaic.fits'
cube_CS_mosaic = SpectralCube.read(fn_CS_mosaic)

cube_CS = cube_CS_mosaic.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=97.98095*u.GHz)

cube_CS_filament = cube_CS.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_CS_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_CS_mos_filament.fits', overwrite=True)

mom0_CS_filament = cube_CS_filament.moment0()
mom0_CS_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_CS_mos_filament.fits', overwrite=True)

# Other lines...

fn_spw25 = fn_s[0]
cube_25 = SpectralCube.read(fn_spw25)
cube_25.allow_huge_operations=True

#H13CN 1-0
print('H13CN')
cube_H13CN = cube_25.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=86.33992*u.GHz)

cube_H13CN_filament = cube_H13CN.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_H13CN_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_H13CN_filament.fits', overwrite=True)

mom0_H13CN_filament = cube_H13CN_filament.moment0()
mom0_H13CN_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_H13CN_filament.fits', overwrite=True)

fn_spw29 = fn_s[1]
cube_29 = SpectralCube.read(fn_spw29)
cube_29.allow_huge_operations=True

# HCO+ 1-0
print('HCO+')
cube_HCOp = cube_29.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=89.18852*u.GHz)

cube_HCOp_filament = cube_HCOp.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_HCOp_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_HCOp_filament.fits', overwrite=True)

mom0_HCOp_filament = cube_HCOp_filament.moment0()
mom0_HCOp_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_HCOp_filament.fits', overwrite=True)

fn_spw27 = fn_s[2]
cube_27 = SpectralCube.read(fn_spw27)
cube_27.allow_huge_operations=True

# SiO 2-1
print('SiO')
cube_SiO = cube_27.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=86.84696*u.GHz)

cube_SiO_filament = cube_SiO.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_SiO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_SiO_filament.fits', overwrite=True)

mom0_SiO_filament = cube_SiO_filament.moment0()
mom0_SiO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_SiO_filament.fits', overwrite=True)

# H13CO+ 1-0
print('H13CO+')
cube_H13COp = cube_27.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=86.7543*u.GHz)

cube_H13COp_filament = cube_H13COp.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_H13COp_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_H13COp_filament.fits', overwrite=True)

mom0_H13COp_filament = cube_H13COp_filament.moment0()
mom0_H13COp_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_H13COp_filament.fits', overwrite=True)

# HN13C 1-0
print('HN13C')
cube_HN13C = cube_27.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=87.09085*u.GHz)

cube_HN13C_filament = cube_HN13C.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_HN13C_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_HN13C_filament.fits', overwrite=True)

mom0_HN13C_filament = cube_HN13C_filament.moment0()
mom0_HN13C_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_HN13C_filament.fits', overwrite=True)
 
fn_spw35 = fn_s[3]
cube_35 = SpectralCube.read(fn_spw35)
cube_35.allow_huge_operations=True

# HC3N 11-10
print('HC3N')
cube_HC3N = cube_35.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=100.0763*u.GHz)

cube_HC3N_filament = cube_HC3N.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_HC3N_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_HC3N_filament.fits', overwrite=True)

mom0_HC3N_filament = cube_HC3N_filament.moment0()
mom0_HC3N_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_HC3N_filament.fits', overwrite=True)

fn_spw33 = fn_s[4]
cube_33 = SpectralCube.read(fn_spw33)
cube_33.allow_huge_operations=True

# CS 2-1
print('CS')
cube_CS = cube_33.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=97.98095*u.GHz)

cube_CS_filament = cube_CS.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_CS_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_CS_filament.fits', overwrite=True)

mom0_CS_filament = cube_CS_filament.moment0()
mom0_CS_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_CS_filament.fits', overwrite=True)

# SO 3(2)-2(1)
print('SO')
cube_SO = cube_33.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=99.29987*u.GHz)

cube_SO_filament = cube_SO.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_SO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_SO_filament.fits', overwrite=True)

mom0_SO_filament = cube_SO_filament.moment0()
mom0_SO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_SO_filament.fits', overwrite=True)

fn_spw31 = fn_s[5]
cube_31 = SpectralCube.read(fn_spw31)
cube_31.allow_huge_operations=True

# HNCO 4-3
print('HNCO')
cube_HNCO = cube_31.with_spectral_unit(u.km/u.s, velocity_convention='radio', rest_value=87.925238*u.GHz)

cube_HNCO_filament = cube_HNCO.subcube_from_regions(reg_filament, minimize=False).spectral_slab(vlow, vhigh)#.to(u.K)
cube_HNCO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/cube_HNCO_filament.fits', overwrite=True)

mom0_HNCO_filament = cube_HNCO_filament.moment0()
mom0_HNCO_filament.write(f'/orange/adamginsburg/jwst/cloudc/alma/moments/mom0_HNCO_filament.fits', overwrite=True)

print('Finished')