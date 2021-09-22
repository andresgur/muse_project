import sys
sys.path.append('/home/mparra/PARRA/Scripts/Python/HST')

from apt_phot import apt_phot
from astropy.wcs import WCS
from astropy.io import fits
from apt_phot import region_to_aperture as rta

import os
from mpdaf.obj import Image

os.chdir('/home/mparra/PARRA/Observ/NGC5408/HST/o3_calc')

file_502='ibde04010_drc_F502N.fits'
file_547='ibde53010_drc_F547N.fits'

bw_502=fits.open(file_502)[1].header['PHOTBW']
bw_547=fits.open(file_547)[1].header['PHOTBW']

wcs_502=WCS(fits.open(file_502)[1].header)
wcs_547=WCS(fits.open(file_547)[1].header)

#global nebula circular region
neb_reg='nebula.reg'
a_neb=rta(neb_reg,wcs_502).area
#ULX counterpart region, visible in both filters
xreg='x-1.reg'
a_xreg=rta(xreg,wcs_502).area

#star complex region near the ULX, visible in both filters
peskystar_reg='peskystar.reg'
a_pstar=rta(peskystar_reg,wcs_502).area

#annoying star only visible in the 547 filter (and thus removed from the continuum computation)
badstar_reg='badstar_547.reg'
a_badstar=rta(badstar_reg,wcs_547).area

#photometry of the nebula with the 0III5007 filter
flux_502_neb=apt_phot(file_502,neb_reg,background_reg=None,apt_corr=1)[6]
flux_502_outelem=apt_phot(file_502,xreg,background_reg=None,apt_corr=1)[6]\
                +apt_phot(file_502,peskystar_reg,background_reg=None,apt_corr=1)[6]
flux_502=flux_502_neb-flux_502_outelem
                
#photometry of the nebula with the continuum filter
flux_547_neb=apt_phot(file_547,neb_reg,background_reg=None,apt_corr=1)[6]
flux_547_outelem=apt_phot(file_547,xreg,background_reg=None,apt_corr=1)[6]\
                +apt_phot(file_547,peskystar_reg,background_reg=None,apt_corr=1)[6]
flux_547_badstar=apt_phot(file_547,badstar_reg,background_reg=None,apt_corr=1)[6]
flux_547=flux_547_neb-flux_547_outelem-flux_547_badstar

#since we take of a part of the area due to the annoying star, we readjust the area to compensate
a_corr=(a_neb-a_xreg-a_pstar)/(a_neb-a_xreg-a_pstar-a_badstar)
flux_547=flux_547*a_corr

#finally we can deduce both fluxes, with the bandwith as a conversion factor
flux_oIII=flux_502-flux_547*(bw_502/bw_547)

#areas in arcsec (manually taken from ds9 in the Analysis->statistics of the region)
aa_neb=5.86457
aa_xreg=0.0769175
aa_pstar=0.0973242
aa_source=aa_neb-aa_xreg-aa_pstar

#flux per arcsec squared : 
fpa=flux_oIII/aa_source

#the output can be used in https://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=MUSE 
#using a gaussian line profile as SED and the first observation setup to get time with the required SNR

#second method : renormalizing with a region with no diffuse emission but a lot of stars

renorm_reg='renorm.reg'

#photometry and ratio computation
flux_502_renorm=apt_phot(file_502,'renorm.reg',background_reg=None,apt_corr=1)[6]
flux_547_renorm=apt_phot(file_547,'renorm.reg',background_reg=None,apt_corr=1)[6]

renorm=flux_502_renorm/flux_547_renorm

flux_oIII_renorm=flux_502-flux_547*renorm
fpa_renorm=flux_oIII_renorm/aa_source

#visual check with MPDAF (which we already know will be imperfect because of the bad star but 
#the area correction should mostly cover that)

img_502=Image(file_502)
img_547=Image(file_547)

#adjusting rotation
img_547_derot=img_547.rotate(img_502.get_rot()-img_547.get_rot())

#resampling and adjusting eventual residual coordinate offsets
img_547_resampl=img_547_derot.align_with_image(img_502)

#renormalisation
img_547_renorm=img_547_resampl*renorm

#doesn't want to substract the images if I don't do this even tho I derotated and aligned IDK why
img_547_renorm.wcs=img_502.wcs

#substracting the continuum
img_502_sub=img_502-img_547_renorm

#hard to see anything in mpl
img_502_sub.plot(scale='linear',vmin=-0.5)

if os.path.exists('502_sub.fits'):
    os.system('rm 502_sub.fits')
img_502_sub.write('502_sub.fits')
