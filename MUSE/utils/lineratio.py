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
from astropy.io import fits

def var_division(a, b, a_err, b_err):
    """Compute the variance on the division of two values or arrays (a/b)
    ----------
    Parameters
    a: value or np.array  of the numerator
    b: value or np.array of the denominator
    a_err: error on the numerator (or array)
    b_err: error on the denominator (or array)
    """

    return (a_err / b) ** 2 + (a * b_err / b ** 2) ** 2


def add_maps_quadrature(linemaps):
    """Add maps in quadrature given a list of files

    Parameters
    ----------
    linemaps: list
        List of strings (paths) of maps to be added
    """
    added_maps_log = ' '.join(linemaps)
    added_data = []
    # add line maps
    print('Adding %s maps together for numerator' % added_maps_log)
    for linemap in linemaps:
        if os.path.isfile(linemap):
            linemapfits = fits.open(linemap)
            if linemapfits[0].data is not None:
                extension = 0
            else:
                extension = 1
            if len(added_data) == 0:
                added_data = linemapfits[extension].data**2
            else:
                added_data += linemapfits[extension].data**2
        else:
            warnings.warn("Line map %s not found." % linemap)
            continue

    return np.sqrt(added_data), added_maps_log


def add_maps(linemaps):
    """Add maps given a list of files

    Parameters
    ----------
    linemaps: list
        List of strings of maps to be added

    """
    added_maps_log = ' '.join(linemaps)
    added_data = []
    # add line maps
    print('Adding %s maps together for numerator' % added_maps_log)
    for linemap in linemaps:
        if os.path.isfile(linemap):
            linemapfits = fits.open(linemap)
            if linemapfits[0].data is not None:
                extension = 0
            else:
                extension = 1
            if len(added_data) == 0:
                added_data = linemapfits[extension].data
            else:
                added_data += linemapfits[extension].data
        else:
            warnings.warn("Line map %s not found." % linemap)
            continue

    return added_data, added_maps_log

# read arguments
ap = argparse.ArgumentParser(description='Computes the line ratio between two maps. If associated error maps are given, it also propagates the uncertainties into the final image (variance stored as EXTENSION = STAT)')
ap.add_argument("-n", "--numerators", nargs='+', help="Line maps to be added in the numerator of the line ratio",
                required=True)
ap.add_argument("-en", "--enumerators", nargs='+', help="(Optional) Associated error line maps of the numerator(s) of the line ratio",
                required=False)
ap.add_argument("-d", "--denominator", nargs='+', help="The line maps to be added in the denominator of the line ratio",
                required=True)
ap.add_argument("-ed", "--edenominator", nargs='+', help="(Optional) Associated error line maps of the denominator(s) of the line ratio",
                required=False)
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='lineratios')
ap.add_argument("-f", "--file", nargs='?', help="Output file name (without fits ending)", default='lineratio.fits')
args = ap.parse_args()

if args.enumerators is not None:
    if len(args.numerators)!=len(args.enumerators):
        raise ValueError("The number of error maps (%d) and line maps (%s) has to be the same for the numerator(s)!" %(len(args.enumerators), len(args.numerators)))

if args.edenominator is not None:
    if len(args.denominator)!=len(args.edenominator):
        raise ValueError("The number of error maps (%d) and line maps (%s) has to be the same for the denominator(s)!"%(len(args.edenominator), len(args.denominator)))

linemaps_numerator = args.numerators
elinemaps_numerator = args.enumerators # associated error map(s)
linemaps_denominator = args.denominator
elinemaps_denominator = args.edenominator # associated error map(s)

outdir = args.outdir
outname = args.file

if not os.path.isdir(outdir):
    os.mkdir(outdir)


numerator_data, added_maps_numerator = add_maps(linemaps_numerator)

if not numerator_data: # list is empty maps were not found
    raise ValueError("Numerator map %s not found!" % added_maps_numerator)

denominator_data, added_maps_denominator = add_maps(linemaps_denominator)

if not denominator_data: # list is empty maps were not found
    raise ValueError("Numerator map %s not found!" % added_maps_numerator)

# get the header from one of the files (or the next one in case it doesn't exist)
if os.path.isfile(linemaps_denominator[0]):
    header = fits.open(linemaps_denominator[0])[0].header
else:
    header = fits.open(linemaps_denominator[1])[0].header

ratio_fits = fits.HDUList(fits.PrimaryHDU(header=header))
ratio_fits[0].header["EXTNAME"] = 'DATA'
ratio_fits[0].header["ERRDATA"] = 'STAT'

ratio_fits.append(fits.ImageHDU(data=numerator_data / denominator_data, header=header,
                         name="DATA"))
ratio_fits[1].header['COMMENT'] = "Ratio of %s/%s line maps" % (added_maps_numerator, added_maps_denominator)


if args.enumerators is not None and args.edenominator is not None:
    err_numerator, added_maps_numerator = add_maps_quadrature(elinemaps_numerator)
    err_denominator, added_maps_denominator = add_maps_quadrature(elinemaps_denominator)
    err_ratio = var_division(numerator_data, denominator_data, err_numerator, err_denominator)
    err_image = fits.ImageHDU(data=err_ratio, header=header, name="STAT")
    err_image.header['COMMENT'] = "Ratio of %s/%s line maps" % (added_maps_numerator, added_maps_denominator)
    ratio_fits.append(err_image)

#ratio_fits.data[np.where(denominator_data is None)] = None

if ratio_fits[0].header['WCSAXES'] == 3:
    for hdu in ratio_fits:
        hdu.header['WCSAXES'] = 2  # set number of axes 2 for the image

ratio_fits.writeto(outdir + "/" + outname + ".fits", overwrite=True)
print('Line ratio %s/%s written to %s/%s.fits' % (added_maps_numerator, added_maps_denominator, outdir, outname))
