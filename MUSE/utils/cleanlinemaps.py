# @Author: Andrés Gúrpide <agurpide>
# @Date:   28-06-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 06-05-2021
# Script to clean the line maps based on a given SNR threshold

# imports
import argparse
import os
import sys
from astropy.io import fits
import numpy as np
import glob
from scipy.ndimage import gaussian_filter


def clean_images(images, snrmap, outdir="cleaned_images", smooth=False, sigma=1, **kargs):

    for image in images:
        if os.path.isfile(image):
            hdul = fits.open(image)
            if smooth:
                hdul[0].data = gaussian_filter(hdul[0].data, sigma=sigma, **kargs)

            if "vel" in image or 'z' in image:
                hdul[0].data[np.where(np.isnan(snr_map[0].data))] = np.nan
                hdul[0].header['COMMENT'] = "Used map %s to filter this image" % snrmap.filename
            else:
                hdul[0].data[np.where(hdul[0].data == 0)] = np.nan
                hdul[0].data[np.where(snr_map[0].data < sthreshold)] = np.nan
                hdul[0].header['COMMENT'] = "Used a signal to noise ratio of %.2f to filter this image" % sthreshold
            if hdul[0].header["WCSAXES"] == 3:
                hdul[0].header["WCSAXES"] = 2
            hdul.writeto("%s/clean%s" % (outdir, os.path.basename(image)), overwrite=True)
            print("Cleaned and stored %s image" % image)
            hdul.close()
        else:
            print("Image %s not found." % image)


# read arguments
ap = argparse.ArgumentParser(description='Clean line maps by thresholding them with the SNR maps.')
ap.add_argument("-r", "--rootname", nargs='?', help="Root name of the input files", type=str, default="")
ap.add_argument("-l", "--line", nargs='?', help='Line to clean maps', type=str)
ap.add_argument("-i", "--image", nargs='*', help='Image to clean', type=str)
ap.add_argument("--smooth", help='Whether to smooth the images prior to thresholding them', action='store_true')
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='cleaned_images')
ap.add_argument("-t", "--threshold", nargs='?', help="Threshold signal to noise ratio to clean the input images", default=5, type=float)
ap.add_argument("-s", "--signaltonoiseratiomap", nargs='?', help="Signal to noise ratio map to clean the input images", default="")
args = ap.parse_args()

rootname = args.rootname
line = args.line
sthreshold = args.threshold
outdir = args.outdir
file_map = args.signaltonoiseratiomap

if args.smooth:
    outdir += "_smoothed"

if not os.path.isdir(outdir):
    os.mkdir(outdir)

if file_map != "" and os.path.isfile(file_map):
    snr_map = fits.open(file_map)

elif rootname != "":
    file_map = glob.glob('./%s*snr*%s.fits' % (rootname, line))
    if len(file_map) == 0:
        print("SNR map with rootname %s for line %s not found!" % (rootname, line))
        sys.exit()
    elif os.path.isfile(file_map[0]):
        print("Found %s SNR map" % file_map[0])
        snr_map = fits.open(file_map[0])

    else:
        print("Signal to noise ratio map %s does not exist!" % file_map)
        sys.exit()
else:
    print("Unable to load signal noise ratio map %s. Provide a rootname or a file!" % file_map)
    sys.exit()

print('Loaded signal to noise ratio map %s successfully' % file_map)

if rootname != "":
    images_toclean = glob.glob('./%s*[!snr]*%s.fits' % (rootname, line))
    print(images_toclean)

    if len(images_toclean) == 0:
        print("No line maps found with root %s for line %s" % (rootname, line))
    else:
        clean_images(images_toclean, snr_map, outdir, smooth=args.smooth, truncate=3)

if args.image is not None:
    clean_images(args.image, snr_map, outdir, smooth=args.smooth, truncate=3)


snr_map.close()
