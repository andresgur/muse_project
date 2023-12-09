# @Author: Andrés Gúrpide <agurpide>
# @Date:   03-09-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 03-08-2021
# Script to create line map ratios out of two line map images fits file

# imports
import argparse
import os
import sys
import numpy as np
import warnings
import glob
from astropy.io import fits
from mpdaf.obj import Cube


def err_division(a, b, a_err, b_err):
    """Compute the variance on the division of two values or arrays (a/b)
    ----------
    Parameters
    a: value or np.array  of the numerator
    b: value or np.array of the denominator
    a_err: error on the numerator (or array)
    b_err: error on the denominator (or array)
    """

    return np.sqrt((a_err / b) ** 2 + (a * b_err / b ** 2) ** 2)


camel_root = "camel_*"
keywords = ["CRVAL3", "CRPIX3", "CD3_3", "CD3_1", "CD3_2", "CD1_3", "CD2_3", "CTYPE3", "CUNIT3", "CRDER3"]

# read arguments
ap = argparse.ArgumentParser(description='Computes the equivalent width (For a line AA = Line Flux/Continuum @ AA). It also propagates the uncertainties into the final image (variance stored as EXTENSION = STAT)')
ap.add_argument("lines", nargs='+', help="Lines(s) for which the equivalent widht is to be computed")
args = ap.parse_args()

continuum_cube = glob.glob("%s_contcube.fits" % camel_root)
if not len(continuum_cube):
    raise ValueError("Continuum cube not found! Check the root file is correct (%s) or your current directory" % camel_root)
else:
    cont_file = continuum_cube[0]
    continuum_cube = Cube(cont_file)
    econt_file = glob.glob("%s_econtcube.fits" % camel_root)[0]
    econtinuum_cube = Cube(econt_file)

    # take the average values and propagate uncertainties over the continuum
    cont_image = continuum_cube.mean(axis=0).data.data
    econt_image = np.sqrt((econtinuum_cube.data.data**2).sum(axis=0) / len(continuum_cube.wave.coord()))

for line in args.lines:
    flux_file = glob.glob("%s_flux_*%s.fits" % (camel_root, line))
    if not len(flux_file): # list is empty maps were not found
        print("Flux map for %s not found! Skipping" % line)
        continue
    eflux_file = glob.glob("%s_eflux_*%s.fits" % (camel_root, line))

    flux_data_hdu = fits.open(flux_file[0])
    flux_image = flux_data_hdu[0].data
    eflux_data_hdu = fits.open(eflux_file[0])
    eflux_image = eflux_data_hdu[0].data

    eqw_image = flux_image / cont_image
    eeqw_image = err_division(flux_image, cont_image, eflux_image, econt_image)

    # value map
    eqw_fits = fits.HDUList(fits.PrimaryHDU(header=flux_data_hdu[0].header))
    eqw_fits[0].header["EXTNAME"] = 'DATA'
    eqw_fits.append(fits.ImageHDU(data=eqw_image,
                    header=flux_data_hdu[0].header, name="DATA"))
    del eqw_fits[0].header["ERRDATA"]
    eqw_fits[1].header['COMMENT'] = "Ratio of %s/%s line maps" % (cont_file, flux_file[0])
    eqw_fits[1].header['WCSAXES'] = 2
    # error map
    eeqw_fits = fits.HDUList(fits.PrimaryHDU(header=flux_data_hdu[0].header))
    eeqw_fits[0].header["EXTNAME"] = 'DATA'
    eeqw_fits.append(fits.ImageHDU(eeqw_image, header=flux_data_hdu[0].header,name="DATA"))
    eeqw_fits[1].header['COMMENT'] = "Ratio of %s/%s line maps" % (econt_file, eflux_file[0])
    eeqw_fits[1].header['WCSAXES'] = 2
    del eeqw_fits[0].header["ERRDATA"]

    for hdu in eqw_fits:
        for key in keywords:
            if key in hdu.header:
                del hdu.header[key]

    for hdu in eeqw_fits:
        for key in keywords:
            if key in hdu.header:
                del hdu.header[key]

    outroot = cont_file.replace("_contcube.fits", "")
    outname = "%s_eqwidth_%s" % (outroot, line)

    eqw_fits.writeto("%s.fits" % (outname), overwrite=True)
    eeqw_fits.writeto("%s_eeqwidth_%s.fits" % (outroot, line), overwrite=True)
    print('Equivalent width %s/%s written to %s.fits' % (cont_file, flux_file[0], outname))
