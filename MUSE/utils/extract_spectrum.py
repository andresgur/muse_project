# @Author: Andrés Gúrpide <agurpide>
# @Date:   02-09-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 12-07-2021
# Script to extract a spectrum from a certain region around a center of the cube. Background subtraction is also possible taking into acocunt the scaling of the areas

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
import numpy as np
from regions import read_ds9, DS9RegionParserError


def extract_spectrum(cube, reg_file, mode="sum"):
    """Extract spectrum from a cube given a region (from read_ds9) and a mode.

    Parameters
    ----------
    cube: mpdaf.obj.Cube
        Cube to extract the spectrum from
    reg: str
        Region file
    mode: str
        Mode can be sum for integrated spectrum, average or psf
    """
    try:
        reg = read_ds9(reg_file)[0]

        centerra = reg.center.fk5.ra
        centerdec = reg.center.fk5.dec

        region_type = type(reg).__name__

        subcube = cube.copy()

        if region_type == 'CircleSkyRegion':
            radius = reg.radius
            try:
                subcube = subcube.subcube_circle_aperture((centerdec.value, centerra.value), radius.value,
                                                          lbda=None, unit_center=centerra.unit,
                                                          unit_radius=radius.unit, unit_wave=cube.wave.unit)
                logger.info("Extracting spectrum from cube %s around position (%s) for region %s" % (cubefile,
                            (centerra, centerdec), region_type))
            except Exception:
                logger.error("Exception occurred", exc_info=True)
                return None, None

        elif region_type == 'EllipseSkyRegion':
            width = reg.width
            height = reg.height
            angle = reg.angle

            subcube.mask_ellipse((centerdec.value, centerra.value), (width.value, height.value), angle.value, inside=False,
                                 lmin=None, lmax=None, unit_center=centerra.unit, unit_radius=height.unit)
            logger.info("Extracting spectrum from cube %s around position (%s) for region %s" % (cubefile,
                        (centerra, centerdec), region_type))
        elif region_type == 'CircleAnnulusSkyRegion':
            inner_radius = reg.inner_radius
            outer_radius = reg.outer_radius
            subcube.mask_region((centerdec.value, centerra.value), inner_radius.value, inside=True,
                                 lmin=None, lmax=None, unit_center=centerra.unit, unit_radius=inner_radius.unit)
            subcube.mask_region((centerdec.value, centerra.value), outer_radius.value, inside=False,
                                 lmin=None, lmax=None, unit_center=centerra.unit, unit_radius=inner_radius.unit)
            logger.info("Extracting spectrum from cube %s around position (%s) for region %s" % (cubefile,
                        (centerra, centerdec), region_type))

    except DS9RegionParserError as e:
        print(e)
        logger.warning("Region %s not implemented yet. Will use a mask to approximate the region" % region_type)
        mask = mu.region_to_mask(cube[0, :, :], region_file)
        # mask every pixel of the cube with the mask
        subcube.mask[:, :, :] = mask

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
    return spe, subcube


def make_extraction_image(image, source_region, background_region=None, outdir="."):

    img_figure, ax = plt.subplots(1, subplot_kw={'projection': image.wcs.wcs})

    image.plot(ax=ax, scale='linear', colorbar="v", show_xlabel=True, show_ylabel=True, zscale=True)
    if source_region is not None:
        mu.plot_regions(source_region, ax, image.data_header)
    if background_region is not None:
        mu.plot_regions(background_region, ax, image.data_header)
    im = ax.images
    cb = im[-1].colorbar
    cb.ax.set_ylabel("Flux (%s)" % image.unit)

    ax.set_xlabel('Ra')
    ax.set_ylabel('Dec')
    img_figure.savefig("%s/extraction_regions.png" % outdir)


# read arguments
ap = argparse.ArgumentParser(description='Parameters for spectrum extraction of a subcube region of the input cube')
ap.add_argument("input_files", nargs='+', help="List of fits muse data cubes to be loaded")
ap.add_argument("-r", "--region", nargs=1, help="Region file with the region to be extracted (only one). Output names are taken from the region", type=str)
ap.add_argument("-b", "--background", nargs="?", help="Background region", type=str, default=None)
ap.add_argument("-o", "--outdir", nargs='?', help="Name of the output directory", default="subcubes", type=str)
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

if not os.path.isfile(region_file):
    logger.error("Region file %s not found" % region_file)
    sys.exit()

if not os.path.isdir(outdir):
    os.mkdir(outdir)

for cubefile in muse_cubes:

    outname = 'subcube%s' % (region_file.replace(".reg", ""))

    outcubename = os.path.basename(cubefile).replace(".fits", "")
    if os.path.isfile(cubefile):
        logger.debug('Loading cube %s ...' % cubefile)
        cube = Cube(cubefile)
        logger.debug('Loading successful')
        print(cube)

    else:
        logger.error("Cube %s not found. Skipping..." % cubefile)
        continue

    source_spe, source_subcube = extract_spectrum(cube, region_file, args.mode)
    # create white light image
    if white_light_image is None:
        white_light_image = source_subcube.sum(axis=0)
        white_light_image.write("%s/%s%s_img.fits" % (outdir, outcubename, outname))
        white_light_image = None # free memory
    source_subcube.write("%s/%s%s_subcube.fits" % (outdir, outcubename, outname))

    # background region is present, subtract it from the source spectrum
    if args.background is not None:
        back_spe, bkg_subcube = extract_spectrum(cube, args.background, args.mode)
        bkg_subcube_wl = bkg_subcube.sum(axis=0)
        bkg_subcube_wl.write("%s/%s%s_bkgimg.fits" % (outdir, outcubename, outname))
        # free memory
        bkg_subcube = None
        bkg_subcube_wl = None

        logger.info("Subtracting scaled background...")
        bkg_reg = read_ds9(args.background)[0]
        bkg_area = mu.region_to_aperture(bkg_reg, cube.wcs.wcs).area
        source_reg = read_ds9(region_file)[0]
        source_area = mu.region_to_aperture(source_reg, cube.wcs.wcs).area
        corrected_spe = source_spe - back_spe * source_area / bkg_area
        back_spe.write("%s/%s%s_bkgspec.fits" % (outdir, outcubename, outname))
        corrected_spe.write("%s/%s%s_source_subspec.fits" % (outdir, outcubename, outname))
        plt.figure()
        corrected_spe.info()
        corrected_spe.plot(title="Extracted spectrum for region %s, subtracted background from %s" % (region_file, args.background))
        plt.savefig("%s/%s%s_source_sub_spe.png" % ((outdir, outcubename, outname)))
    if source_spe is None and source_subcube is None:
        logger.error("Region not implemented yet, skipping cube %s" % cubefile)
        continue

    source_subcube = None # free memory
    make_extraction_image(cube.sum(axis=0), region_file, args.background, outdir)
    cube = None # free memory

    plt.figure()
    source_spe.info()
    source_spe.plot(title="Extracted spectrum for region %s" % region_file)
    plt.savefig("%s/%s%s_source_spe.png" % ((outdir, outcubename, outname)))

    logger.debug('Writing outputs to %s...' % outdir)
    source_spe.write("%s/%s%s_sourcespec.fits" % (outdir, outcubename, outname))
    mu.plot_image("%s/%s%s_whiteimage" % (outdir, outcubename, outname), "%s/%s%s_img.fits" % (outdir, outcubename, outname))
    mu.plot_image("%s/%s%s_bkgimg" % (outdir, outcubename, outname), "%s/%s%s_bkgimg.fits" % (outdir, outcubename, outname))

    if figure_comparison is not None:
        figure_comparison.savefig("%s/%s%sextraction_region.pdf" % (outdir, outcubename, outname))
