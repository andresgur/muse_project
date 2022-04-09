# @Author: Andrés Gúrpide <agurpide>
# @Date:   28-06-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 09-04-2022



# -*- coding: utf-8 -*
# Script to clean images using a signal-to-noise ratio
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 20-04-2019
    #!/usr/bin/env python3

# imports
import argparse
import os
from astropy.io import fits
import numpy as np
import glob
from scipy.ndimage import gaussian_filter


def clean_images(images, snrmap, outdir="cleaned_images", smooth=False, sigma=1, **kargs):

    for image in images:
        if os.path.isfile(image):
            hdul = fits.open(image)
            if hdul[0].data is not None:
                extension = 0
            else:
                extension = 1
            if snr_map[0].data is not None:
                snr_extension = 0
            else:
                snr_extension = 1
            if smooth:
                hdul[extension].data = gaussian_filter(hdul[extension].data, sigma=sigma, **kargs)

            if "vel" in image or 'z' in image:
                hdul[extension].data[np.where(np.isnan(snr_map[extension].data))] = np.nan
                hdul[extension].header['COMMENT'] = "Used map %s to filter this image" % snrmap.filename()
            else:
                hdul[extension].data[np.where(hdul[extension].data == 0)] = np.nan
                hdul[extension].data[np.where(snr_map[snr_extension].data < sthreshold)] = np.nan
                hdul[extension].header['COMMENT'] = "Used map %s to filter this image" % snrmap.filename()
                hdul[extension].header['COMMENT'] = "Used a signal to noise ratio of %.2f to filter this image" % sthreshold

            if "WCSAXES" in hdul[0].header:
                if hdul[0].header["WCSAXES"] == 3:
                    hdul[0].header["WCSAXES"] = 2
            # this keyword makes some functions crash beacuse the image is interpreted as a 3D cube
            if "CRDER3" in hdul[0].header:
                del hdul[0].header["CRDER3"]
            hdul.writeto("%s/clean_%s" % (outdir, os.path.basename(image)), overwrite=True)
            print("Cleaned and stored %s image" % image)
            hdul.close()
        else:
            print("Image %s not found." % image)


# read arguments
ap = argparse.ArgumentParser(description='Clean line maps by thresholding them with the SNR maps.')
ap.add_argument("-r", "--rootname", nargs='?', help="Root name of the input files", type=str, default="")
ap.add_argument("-l", "--line", nargs='?', help='Line to clean maps', type=str)
ap.add_argument("-i", "--image", nargs='*', help='Images to clean', type=str)
ap.add_argument("--smooth", help='Whether to smooth the images prior to thresholding them', action='store_true')
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='cleaned_images')
ap.add_argument("-t", "--threshold", nargs='?', help="Threshold signal to noise ratio to clean the input images", default=5, type=float)
ap.add_argument("-s", "--signaltonoiseratiomap", nargs='?', help="Signal to noise ratio map to clean the input images", default="")
args = ap.parse_args()

rootname = args.rootname
line = args.line
sthreshold = args.threshold
outdir = args.outdir
snr_input_map = args.signaltonoiseratiomap

if args.smooth:
    outdir += "_smoothed%d" % 2

if not os.path.isdir(outdir):
    os.mkdir(outdir)

snr_map = None
# By default try to locate the root snr file name
if rootname != "":
    file_map = glob.glob('./%s*snr*%s.fits' % (rootname, line))
    if len(file_map) == 0:
        raise ValueError("SNR map with rootname %s for line %s not found!" % (rootname, line))
    elif os.path.isfile(file_map[0]):
        print("Found %s SNR map" % file_map[0])
        snr_map = fits.open(file_map[0])

# If provided use the input one
if snr_input_map != "":
    if os.path.isfile(snr_input_map):
        print("Found new input %s SNR map. The map will be cleaned using this SNR map then." % snr_input_map)
        snr_map = fits.open(snr_input_map)
    else:
        raise ValueError("Unable to locate signal noise ratio map %s. Provide a rootname or a file!" % snr_input_map)

elif snr_map is None:
    raise ValueError('Unable to locate any SNR map. Provide a rootname or an input SNR map!')

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
