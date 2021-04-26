# @Author: Andrés Gúrpide <agurpide>
# @Date:   03-09-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 26-04-2021
# Script to create line map ratios out of two line map images fits file

# imports
import argparse
import os
import sys
import numpy as np
import logging
from astropy.io import fits


def add_maps(linemaps):
    """Add maps given a list of files

    Parameters
    ----------
    linemaps: list
        List of strings of maps to be added

    """
    added_maps_log = ' '.join(linemaps_numerator)
    added_data = []
    # add line maps
    logging.info('Adding %s maps together for numerator' % added_maps_log)
    for linemap in linemaps:
        if os.path.isfile(linemap):
            linemapfits = fits.open(linemap)
            if len(added_data) == 0:
                added_data = linemapfits[0].data
            else:
                added_data += linemapfits[0].data
        else:
            logging.warning("Line map %s not found." % linemap)
            continue

    return added_data, added_maps_log


# read arguments
ap = argparse.ArgumentParser(description='Spectrum fits to be loaded')
ap.add_argument("-n", "--numerators", nargs='+', help="The line maps to be added in the numerator of the line ratio")
ap.add_argument("-d", "--denominator", nargs=1, help="The line maps to be added in the denominator of the line ratio")
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='lineratios')
ap.add_argument("-f", "--file", nargs='?', help="Output file name (without fits ending)", default='lineratio.fits')
args = ap.parse_args()

linemaps_numerator = args.numerators
linemaps_denominator = args.denominator
outdir = args.outdir
outname = args.file

if not os.path.isdir(outdir):
    os.mkdir(outdir)
numerator_data, added_maps_numerator = add_maps(linemaps_numerator)

if len(numerator_data) == 0:
    logging.error("No line maps were found for the numerator")
    sys.exit()

denominator_data, added_maps_denominator = add_maps(linemaps_denominator)

if len(denominator_data) == 0:
    logging.error("No line maps were found for the numerator")
    sys.exit()
# get the header from one of the files (or the next one in case it doesn't exist)
if os.path.isfile(linemaps_denominator[0]):
    header = fits.open(linemaps_denominator[0])[0].header
else:
    header = fits.open(linemaps_denominator[1])[0].header

ratio_fits = fits.PrimaryHDU(data=numerator_data / denominator_data, header=header)
ratio_fits.data[np.where(denominator_data is None)] = None

ratio_fits.header['COMMENT'] = "Ratio of %s/%s line maps" % (added_maps_numerator, added_maps_denominator)
if ratio_fits.header['WCSAXES'] == 3:
    ratio_fits.header['WCSAXES'] = 2  # set number of axes 2 for the image

ratio_fits.writeto(outdir + "/" + outname + ".fits", overwrite=True)
print('Line ratio %s/%s written to %s/%s.fits' % (added_maps_numerator, added_maps_denominator, outdir, outname))
