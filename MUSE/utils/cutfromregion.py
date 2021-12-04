# @Author: Andrés Gúrpide <agurpide>
# @Date:   03-07-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 04-12-2021

# imports
import argparse
import os
from mpdaf.obj import Image, Cube
import logging
import muse_utils as mu
import numpy as np

ap = argparse.ArgumentParser(description='Cuts Cube or Image from a region file')
ap.add_argument("input_file", nargs=1, help="Image or cube to be cut")
ap.add_argument("-o", "--outname", nargs='?', help="Output file name, autonaming otherwise (cutinput_file.fits)", default=None)
ap.add_argument("-r", "--region", nargs=1, help='Region file to cut the cube')

args = ap.parse_args()

# logger
scriptname = os.path.basename(__file__)
logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# handler
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(out_format)
logger.addHandler(stream_handler)

region_file = args.region[0]
input_file = args.input_file[0]

try:
    # is the input a cube?
    input_type = 'cube'
    logger.info("Trying to load input %s as a %s..." % (input_file, input_type))
    input_mpdaf = Cube(input_file)
    img = input_mpdaf[0, :, :]
    logger.info("Masking pixels from region %s ..." % region_file)
    mask = mu.region_to_mask(img, region_file)
    img = None  # free memory

    # mask every pixel of the cube with the mask
    input_mpdaf.mask[:, ] = mask
    #mask = None  # free memory
    # get rid of the masked values
    logger.info("Cropping cube...")
    input_mpdaf.crop()
    #input_mpdaf.unmask()

except ValueError:

    input_type = 'image'
    # it is an image
    logger.info("Trying to load input %s as an %s..." % (input_file, input_type))
    input_mpdaf = Image(input_file, ext=1)
    logger.info("Cropping %s with region %s" % (input_file, region_file))
    mask = mu.region_to_mask(input_mpdaf, region_file)
    #input_mpdaf.data[mask] = np.ma.masked
    #input_mpdaf.crop()
    input_mpdaf[np.where(mask)] = True
    input_mpdaf.crop()
    input_mpdaf.unmask()

input_mpdaf.primary_header["HISTORY"] = "Cut with region %s" % region_file
input_mpdaf.write("cut%s" % os.path.basename(input_file))
logger.info("Written cut fits to cut%s " % os.path.basename(input_file))
