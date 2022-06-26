
# -*- coding: utf-8 -*-
# Script to extract a spectrum from a certain region (circle) around a center of the cube. \
# The final spectrum is the mean spectrum of the region
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 25-03-2019
# !/usr/bin/env python3
# imports
from mpdaf.obj import Cube, Image
from mpdaf.sdetect import compute_optimal_spectrum
import argparse
import os
import sys
import matplotlib.pyplot as plt
import muse_utils as mu
import logging
import mpdaf.MUSE.PSF as psf
from regions import read_ds9
import numpy as np

# read arguments
ap = argparse.ArgumentParser(description='Parameters for spectrum extraction of a subcube region of the input cube')
ap.add_argument("input_files", nargs='+', help="List of fits muse data cubes to be loaded")
ap.add_argument("-r", "--region", nargs=1, help="Region file with the regions to be extracted. Label regions for autonaming", type=str)
ap.add_argument("-o", "--outdir", nargs='?', help="Name of the output directory", default="subcubes", type=str)
ap.add_argument("--lmin", help='Minimum wavelength to cut the cube spectrally', type=float, nargs='?', default=0)
ap.add_argument("--lmax", help='Maximum wavelength to cut the cube spectrally', type=float, nargs='?', default=100000)
ap.add_argument("--psf", help='Fits file containing the FWHM as a function of wavelength', type=str, nargs='?', default=None)
ap.add_argument("--e_fraction", help='Energy fraction to be enclosed in the extraction region', type=float, nargs='?', default=0.5)
ap.add_argument("-mode", choices=["sum", "mean", "psf"], default="sum")
ap.add_argument("--ext", help='Extension to read from the PSF fits file', type=int, nargs='?', default=2)

args = ap.parse_args()

muse_cubes = args.input_files
region_file = args.region[0]
outdir = args.outdir
extraction_mode = args.mode
extension = args.ext

e_fraction = args.e_fraction

scriptname = os.path.basename(__file__)
# logger

logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# handler
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(out_format)
logger.addHandler(stream_handler)

figure_comparison = None
white_light_image = None

if os.path.isfile(region_file):
    regs = read_ds9(region_file)
    logger.info("Loaded %i region(s) from %s" % (len(regs), region_file))
    logger.info(regs)
else:
    logger.error("Region file %s not found" % region_file)
    sys.exit()

if not os.path.isdir(outdir):
    os.mkdir(outdir)


for cubefile in muse_cubes:

    if os.path.isfile(cubefile):
        logger.debug('Loading cube %s ...' % cubefile)
        cube = Cube(cubefile)
        logger.debug('Loading successful')

    else:
        logger.warning("Cube %s not found. Skipping..." % cubefile)
        continue

    for reg_index, reg in enumerate(regs):

        centerra = reg.center.fk5.ra
        centerdec = reg.center.fk5.dec

        region_type = type(reg).__name__

        if region_type == 'CircleSkyRegion':
            radius = reg.radius
            try:
                subcube = cube.subcube_circle_aperture((centerdec.value, centerra.value), radius.value,
                                                       lbda=(args.lmin, args.lmax), unit_center=centerra.unit,
                                                       unit_radius=radius.unit, unit_wave=cube.wave.unit).copy()
            except Exception:
                logger.error("Exception occurred", exc_info=True)
                continue

        elif region_type == 'EllipseSkyRegion':
            width = reg.width
            height = reg.height
            angle = reg.angle
            subcube = cube.copy()

            subcube.mask_ellipse((centerdec.value, centerra.value), (width.value, height.value), angle.value, inside=False,
                                 lmin=args.lmin, lmax=args.lmax, unit_center=centerra.unit, unit_radius=height.unit)

        else:
            logger.error("Region %s not implemented yet" % region_type)
            continue

        logger.info("Extracting spectrum number %i from cube %s around position (%s) for region %s" % (reg_index, cubefile,
                    (centerra, centerdec), region_type))

        subcube.crop()
        subcube.info()

        if extraction_mode == 'sum':
            logger.info("Summing output spectra")
            spe = subcube.sum(axis=(1, 2))
        elif extraction_mode == 'mean':
            logger.info("Averaging output spectra")
            spe = subcube.mean(axis=(1, 2))
        elif extraction_mode == 'psf':
            if args.psf is None:
                logger.error("PSF mode requires a fits file PSF")
                sys.exit()
            else:

                white_light_image = subcube.sum(axis=0)

                psf_table, header, headername = mu.read_psf_fits(args.psf, extension)

                if "Beta" in header:
                    avg_beta = float(header['Beta'])
                else:
                    avg_beta = None

                logger.info("Creating PSF cube...")
                psf_cube = psf.create_psf_cube(subcube.shape, psf_table["FWHM"], avg_beta, subcube.wcs)
                psf_moffat = Image(wcs=white_light_image.wcs, data=np.sum(psf_cube, axis=0), cmap='inferno', mask=white_light_image.mask)
                peak_center = psf_moffat.peak()
                energy_square_size_x, energy_square_size_y = psf_moffat.ee_size(center=(peak_center['y'], peak_center['x']), cont=0, frac=e_fraction)
                center_pix = psf_moffat.wcs.sky2pix((peak_center['y'], peak_center['x']))[0]
                radius = energy_square_size_x / 2
                logger.info("Extraction radius determining energy fraction of %.1f is %.1f arcsec" % (e_fraction, radius))
                figure_comparison, axes = plt.subplots(nrows=1, ncols=2, figsize=(16.0, 10.0), subplot_kw={'projection': white_light_image.wcs.wcs},
                                                       gridspec_kw={'wspace': 0, 'hspace': 0}, sharex=True)
                figure_comparison.suptitle(t="Extraction region for region %i" % reg_index)
                moffat_ax = axes[1]

                # we need one region for each plot
                # X is the size of a square so its half of the radius
                ex_region_moffat = plt.Circle(center_pix, radius / white_light_image.get_step('arcsec')[0], color='b', fill=False)
                moffat_ax.add_artist(ex_region_moffat)
                moffat_ax.set_yticklabels([])
                moffat_ax.set_title("PSF")
                psf_moffat.plot(ax=moffat_ax)
                moffat_ax.set_yticklabels([])

                image_ax = axes[0]
                peak_wimage = white_light_image.peak()
                image_ax.set_title("Extraction region")
                image_ax.errorbar(peak_wimage["p"], peak_wimage["q"], lw=0, fmt='x', color='red')
                ex_region_image = plt.Circle(center_pix, radius / white_light_image.get_step('arcsec')[0], color='b', fill=False)
                image_ax.add_artist(ex_region_image)
                white_light_image.plot(ax=image_ax)
                plt.show()
                logger.info("Cutting extraction region encircling %.1f %% of the energy..." % e_fraction)
                subcube = subcube.subcube_circle_aperture((peak_wimage['y'], peak_wimage['x']), radius, unit_center='deg', unit_radius='arcsec',
                                                          unit_wave=cube.wave.unit)
                subcube.crop()
                # cut PSF again with new cube (maybe unnecessary but just so the mask and the cube sizes match)
                psf_cube = psf.create_psf_cube(subcube.shape, psf_table["FWHM"], avg_beta, subcube.wcs)

                logger.info("Extracting spectrum weighted by PSF")
                dummy_mask = np.ones((subcube.shape[1], subcube.shape[2]))  # dummy mask required by the method
                spe = compute_optimal_spectrum(subcube, dummy_mask, psf_cube)
        # create white light image
        if white_light_image is None:
            white_light_image = subcube.sum(axis=0)

        plt.figure()
        spe.info()
        spe.plot(title="Extracted spectrum for region %i" % reg_index)
        plt.show()

        try:
            outname = reg.meta['label']
        # region without text label
        except KeyError:
            logger.warning("Region without name. Set text to the region to autonamed it.")
            outname = 'unnamed_region'

        if outname == "":
            outname = 'subcube%i' % reg_index

        outcubename = os.path.basename(cubefile).replace(".fits", "")
        logger.debug('Writing outputs...')
        spe.write("%s/%s%s_spec.fits" % (outdir, outcubename, outname))
        white_light_image.write("%s/%s%s_img.fits" % (outdir, outcubename, outname))
        subcube.write("%s/%s%s_subcube.fits" % (outdir, outcubename, outname))
        if figure_comparison is not None:
            figure_comparison.savefig("%s/%s%sextraction_region.pdf" % (outdir, outcubename, outname))

        subcue = None # free memory
    cube = None
