# @Author: Andrés Gúrpide <agurpide>
# @Date:   10-06-2019
# @Email:  agurpidelash@soton.ac.uk
# @Last modified by:   agurpide
# @Last modified time: 10-12-2021
# Script to adjust the coordinates of a cube from an input image by cross-correlating it( Preferabley HST)


# imports
import argparse
from mpdaf.obj import Cube, Image
import os
import logging
from astropy.io import fits
import gc
# read arguments
ap = argparse.ArgumentParser(description='Adjust coordinates from HST image')
ap.add_argument("input_cube", nargs=1, help="Cube to be corrected", type=str)
ap.add_argument("--hst", nargs=1, help="HST image", type=str)
ap.add_argument("-o", "--outdir", nargs='?', type=str, help="Output directory", default="coordadjusted")

# parse args
args = ap.parse_args()
nsigma = 3

# create output dir
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

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


input_cube = args.input_cube[0]
hst_image = args.hst[0]
outcube = '%s/corr%s' % (args.outdir, input_cube)

if os.path.isfile(input_cube) and os.path.isfile(hst_image):
    logger.info("Loading cube %s" % input_cube)
    logger.info("Loading succesful")
    hdul = fits.open(hst_image)
    if hdul[0].header["TELESCOP"] == "HST":
        # http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=HST/ACS_WFC.F555W&&mode=browse&gname=HST&gname2=ACS_WFC#filter
        instrument = hdul[0].header["INSTRUME"]
        if "FILTER" in hdul[0].header:
            hst_filter = hdul[0].header["FILTER"]
        elif "FILTNAM1" in hdul[0].header:
            hst_filter = hdul[0].header["FILTNAM1"]
        else:
            hst_filter = hdul[0].header["FILTER1"]
        filter_key = "%s_%s" % (instrument, hst_filter)
        logger.info("Constructing image in the %s filter" % filter_key)
        # free memory
        del hdul        
        gc.collect()
        cube = Cube(input_cube)

        try:
            cube_image = cube.get_band_image(filter_key)
        except ValueError:
            print(f"Filter {filter_key} not found, attempting to use ACS filter version instead")
            filter_key = filter_key.replace("WFPC2", "ACS")
            try:
                cube_image = cube.get_band_image(filter_key)
                print(f"Successfully created {filter_key} image")
            except ValueError:
                print(f"Filter {filter_key} not found, attempting to use WFC3 filter version instead")
                filter_key = filter_key.replace("ACS", "WFC3")
            try:
                cube_image = cube.get_band_image(filter_key)
            except ValueError:
                print(f"Filter {filter_key} not found, reverting to white light image")
                cube_image = cube.sum(axis=0)
                filter_key = "white_light"
        print(f"Successfully created {filter_key} image")
            
    else:
        cube = Cube(input_cube)
        filter_key = "white_light"
        print("Image format not HST, using cube white light image for the coordinate correction")
        cube_image = cube.sum(axis=0)
    del cube
    wlname = "%s/%s" % (args.outdir, input_cube.replace(".fits", "_imgcorr_%s.fits" % filter_key))
    # free memory

    hst_im = Image(hst_image)
    logger.info("Estimating coordinate offset...")
    dy, dx = cube_image.estimate_coordinate_offset(hst_im, nsigma=nsigma)
    cube_image.adjust_coordinates(hst_im, nsigma=nsigma, inplace=True)
    del hst_im
    cube_image.write(wlname)
    # update cube WCS header
    cube = Cube(input_cube)
    cube.wcs.set_crpix1(cube_image.wcs.get_crpix1())
    cube.wcs.set_crpix2(cube_image.wcs.get_crpix2())
    cube.wcs.set_cd(cube_image.wcs.get_cd())
    pixel_scale = cube.primary_header['HIERARCH ESO OCS IPS PIXSCALE']
    cube.primary_header["ASTRO_CORR_PX"] = ("%.2f, %.2f" % (dy, dx), 'Coordinate shift dy, dx in pixels')
    print(dy * pixel_scale)
    cube.primary_header["ASTRO_CORR_AS"] = ("%.2f, %.2f" % ((-dy * pixel_scale), (-dx * pixel_scale)), 'Coordinate shift dy, dx in arcseconds')
    cube.primary_header['COMMENT'] = "Coordinates corrected using %s image" % hst_image
    cube.write(outcube)
    logger.info("Output successfully stored to %s" % args.outdir)
else:
    logger.error("Input cube %s or HST image %s not found" % (input_cube, hst_image))
