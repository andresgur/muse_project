# @Author: Andrés Gúrpide <agurpide>
# @Date:   07-07-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 11-03-2022



# -*- coding: utf-8 -*
# Script to clean a cube from sky residuqls
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 12-04-2019
# !/usr/bin/env python3

# imports
import numpy as np
from mpdaf.obj import Image, Cube
import zap
import argparse
import os
import logging
import muse_utils as mu
import multiprocessing


def create_mask_region(image, percentile):
    """Create mask region from white light image and percentile flux threshold.

    Creates mask region out of a white ligh image.

    Parameters
    ----------
    image : white light image with same size as the cube to be cleaned
    percentile : threshold percentile of the histogram flux to set mask regions
    """
    print("Creating masked region")

    threshold_flux = np.percentile(image.data, percentile)
    print("Masking every pixel above %.2f %s" % (threshold_flux, image.unit.to_string()))
    image.info()
    mask_image = image.clone(data_init=np.empty)

    mask_image[np.where(image.data >= threshold_flux)] = 1
    mask_image[np.where(image.data < threshold_flux)] = 0

    return mask_image


# read arguments
ap = argparse.ArgumentParser(description='Clean input cubes from sky residuals using the ZAP tool')
ap.add_argument("input_cube", nargs=1, help="Muse data cube to be cleaned from sky residuals")
ap.add_argument("-w", "--white_image", nargs='?', help="White light image from the cube to compute the zero order sky residuals spectrum. If not provided it is computed by the tool", type=str)
ap.add_argument("-p", "--percentatge", nargs='?', help="Percentatge threshold  of the white light image flux histogram to be masked \
 for the zero order sky calculation ", default="50", type=float)
ap.add_argument("-r", "--region", help='Region to create the mask for the cube', nargs='?', type=str)
ap.add_argument("-o", "--output", nargs='?', help="Name of the output cleaned cube. Default takes the input name as steem", default="")
ap.add_argument("-e", "--eigenvalues", nargs='?', help="Number of eigenvalues for ZAP. High values risk removing astrophysical signal, while low values will be less aggressive. Default will let zap work out the number needed.",
                default=-1, type=int)
ap.add_argument("-c", "--cores", nargs='?', help="Number of processing cores for the task.", default=-1, type=int)
ap.add_argument("--outdir", nargs='?', help="Name of the output directory", default="zap_cleaned", type=str)
ap.add_argument("--nanclean", help='Clean NaN values prior to cleaning the sky residuals, default false', action='store_true')
# parse args
args = ap.parse_args()

input_cube = args.input_cube[0]
percentile = args.percentatge
output_cube = args.output
mask_region = args.region
nanclean = args.nanclean
# use all cores - 4 by default
cores = multiprocessing.cpu_count() - 4 if args.cores ==-1 else args.cores

# logger
scriptname = os.path.basename(__file__)
# logger

logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# handler
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)
stream_handler.setFormatter(out_format)
logger.addHandler(stream_handler)

# create output dir
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

if args.white_image is None:
    white_light_image = Cube(input_cube).sum(axis=0)
    image_file = "w_light_image.fits"
else:
    white_light_image = Image(args.white_image)
    image_file = args.white_image

out_mask_file = "%s/%s" % (args.outdir, os.path.basename(image_file).replace(".fits", "masked.fits"))
# create mask region from white light image and percentile
if image_file is not None and mask_region is None:
    mask_image = create_mask_region(white_light_image, percentile)
    logger.debug("Saving masked image")
    mask_image.write(out_mask_file)
    # create mask image from region file
elif mask_region is not None:
    mask = mu.region_to_mask(white_light_image, mask_region)
    mask_image = white_light_image.clone(data_init=np.zeros)
    # we mask the values inside the circles!
    mask_image[np.where(mask == False)] = 1
    mask_image.write(out_mask_file)

else:
    mask_file = None
    print("Warning: no image was provided to compute the masked regions")

if output_cube == "":
    output_cube = "%s/%s" % (args.outdir, os.path.basename(input_cube).replace(".fits", "_zapcleaned.fits"))
    out_skycube = "%s/%s" % (args.outdir, os.path.basename(input_cube).replace(".fits", "_sky.fits"))
    out_varcurve = "%s/%s" % (args.outdir, os.path.basename(input_cube).replace(".fits", "_varcurve.fits"))
else:
    out_skycube = "%s/%s" % (args.outdir, output_cube.replace(".fits", "_sky.fits"))
    out_varcurve = "%s/%s" % (args.outdir, output_cube.replace(".fits", "_varcurve.fits"))

logger.debug("Name for the cleaned cube: %s" % output_cube)

# cfwidthSVD default is 300 for eigenvector calculation; good; SP continuum filter ~ 20 - 50 pixels

if args.eigenvalues==-1:
    zap.process(input_cube, skycubefits=out_skycube, mask=out_mask_file,
                outcubefits=output_cube, cfwidthSVD=300, cfwidthSP=300,
                varcurvefits=out_varcurve, overwrite=True, clean=nanclean, ncpu=cores)
else:
    logger.info("Using %d eigenvalues" % args.eigenvalues)
    zap.process(input_cube, skycubefits=out_skycube, mask=out_mask_file,
                outcubefits=output_cube, cfwidthSVD=300, cfwidthSP=300,
                varcurvefits=out_varcurve, overwrite=True, clean=nanclean, ncpu=cores, nevals=[args.eigenvalues])

'''
process.writecube(outcubefits='cube.fits', overwrite=True)
process.writeskycube(skycubefits='skycube.fits', overwrite=True)


process.plotvarcurve()
plt.figure()
plt.plot(process.cube[:, :, :].sum(axis=(1, 2)), 'b', alpha=0.5)
plt.plot(process.cleancube[:, :, :].sum(axis=(1, 2)), 'g')

plt.show()
'''
