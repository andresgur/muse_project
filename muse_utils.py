# @Author: Andrés Gúrpide <agurpide>
# @Date:   29-08-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 22-06-2021



# -*- coding: utf-8 -*-
# Methods to read config line files
# Author: Andres Gurpide Lasheras  andres.gurpide@gmail.com
# 04-04-2019
# !/usr/bin/env python3

# imports
import matplotlib.pyplot as plt
import numpy as np
import os
from regions import read_ds9
from mpdaf.obj import Image
from photutils.aperture import CircularAperture, CircularAnnulus, EllipticalAperture, SkyCircularAperture, SkyCircularAnnulus, SkyEllipticalAperture
import pyregion
import astropy.units as u
from collections import OrderedDict
import sys
from astropy.io import fits


def region_to_aperture(region, wcs=None):
    """Convert region object to photutils.aperture.aperture_photometry object. The wcs object is needed only if the input regions are in sky coordinates.

    Parameters
    ----------
    region: regions.Region
        Output of read_ds9 method
    wcs: astropy.wcs.WCS
        A world coordinate system if the region in sky coordinates IS needed to convert it to pixels.
    """

    region_type = type(region).__name__
    if "Pixel" in region_type:
        source_center = (region.center.x, region.center.y)
        if region_type == 'CirclePixelRegion':
            return CircularAperture(source_center, r=region.radius)
        elif region_type == "CircleAnnulusPixelRegion":
            return CircularAnnulus(source_center, r_in=region.inner_radius, r_out=region.outer_radius)
        elif region_type == "EllipsePixelRegion":
            # to be tested
            return EllipticalAperture(source_center, a=region.width, b=region.height, angle=region.angle)
    elif "Sky" in region_type:
        if wcs is None:
            print("Error, cannot obtain aperture without a wcs.")
            return None
        center = region.center.fk5
        if region_type == "CircleSkyRegion":
            return SkyCircularAperture(center, r=region.radius).to_sky(wcs)
        elif region_type == "CircleAnnulusSkyRegion":
            print("Region %s not implemented")
        elif region_type == "EllipseSkyRegion":
            return SkyEllipticalAperture(center, a=region.width, b=region.height, angle=region.angle).to_sky(wcs)
        elif region_type == "CircleAnnulusSkyRegion":
            return SkyCircularAnnulus(center, r_in=region.inner_radius, r_out=region.outer_radius).to_sky(wcs)
    else:
        print("Error region not implemented")
        return None



def read_psf_fits(fits_file, extension=1):
    """Read PSF file produced by determine_psf script
    Returns astropy table"""
    if not os.path.isfile(fits_file):
        sys.exit()

    hdul = fits.open(fits_file)
    print("Reading data from %s" % hdul[extension].name)
    return hdul[extension].data, hdul[extension].header, hdul[extension].name


def get_atm_seeing(input_mpdaf):
    try:
        return input_mpdaf.primary_header["HIERARCH ESO TEL IA FWHMLIN"] * u.arcsec
    except KeyError:
        return None


def get_sky_res(input_mpdaf):
    try:
        return input_mpdaf.primary_header["SKY_RES"], input_mpdaf.primary_header["SKY_RERR"]
    except KeyError:
        return None, None


def region_to_mask(image, region):
    """Creates a mask out of a region file. Masks everything outside the region.
    Parameters:
    fits : astropy.io.fits
        The fits file you want to create the mask for (Cube or Image)
    region : str
        The region file with which you want to create the mask (path)
    returns a ndarray with 0 for masked values and 1 for non-masked values
    """
    mask_region = pyregion.open(region)
    print(mask_region)
    mask = mask_region.get_mask(hdu=image.get_data_hdu())
    mask = ~mask
    return mask


def load_map_file(file):
    """Load config file for line map with line names and min and max wavelenghts to be masked.
    Parameters
    ----------
    file : the config file containing the lines to be processed."""
    with open(file) as line_file:
        info_lines = np.loadtxt(line_file, dtype='str', delimiter='\t\t')
        line_names = info_lines[:, 0]
        min_wa = np.array(info_lines[:, 1])
        max_wa = np.array(info_lines[:, 2])
        print("Lines found in %s \n" % file)
        print(line_names)
        return line_names, min_wa, max_wa


def load_spectrum_file(spectrum_file):
    """Load config file for spectrum fitting.
    Parameters
    ----------
    spectrum_file : the config file containing the lines to be fitted."""
    with open(spectrum_file) as config_file:
        info_lines = np.loadtxt(config_file, dtype='str', delimiter='\t\t')
        line_names = info_lines[:, 0]
        min_wa = np.array(info_lines[:, 1])
        max_wa = np.array(info_lines[:, 2])
        line_type = info_lines[:, 3]
        fwhms = info_lines[:, 5]
        return line_names, min_wa, max_wa, line_type, fwhms


def add_zlabel(image_file, units):
    '''Decide the z label for the input image based on the file name or the units.'''

    if units != "None":
        zlabel = units

    if 'flux' in image_file:
        zlabel = "1e-20 erg/cm$^2$/s"
    elif 'wave' in image_file:
        if 'disp' in image_file:
            zlabel = '$\Delta\lambda (\AA)$'
        else:
            zlabel = '$\lambda (\AA)$'

    elif 'disp' in image_file:
        zlabel = '$\Delta$v (km/s)'

    elif 'vel' in image_file:
        zlabel = 'v (km/s)'
    elif 'snr' in image_file:
        zlabel = '$\sigma$'
    if 'ratio' in image_file:
        zlabel = ""

    return zlabel


def create_grid(outname, *image_files, region=None, cutting_region=None, vmin=None, nrows=1, ncols=1, colorbar='v'):

    img_figure, axes = plt.subplots(nrows, ncols, figsize=(16, 10), sharey=True,
                                    subplot_kw={'projection': Image(image_files[0]).wcs.wcs})
    # add a big axes, hide frame
    plt.rcParams.update({'font.size': 20})
    img_figure.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.xlabel("Ra")
    plt.ylabel("Dec")
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    plt.grid(False)

    for image_file, ax in zip(image_files, axes):

        image = Image(image_file)

        if cutting_region is not None:
            mask = region_to_mask(image, cutting_region)
            image[np.where(mask)] = True
            image.crop()

        image.plot(ax=ax, scale='linear', colorbar=colorbar, show_xlabel=False, show_ylabel=False,
                   zscale=False, vmin=vmin)

        if region is not None:
            plot_regions(region, ax)

        if colorbar is not None:
            im = ax.images
            cb = im[-1].colorbar
            cb.ax.set_ylabel(add_zlabel(image_file, image.unit), fontsize=14)
        format_axis(ax, 15)

    img_figure.subplots_adjust(hspace=0)
    plt.show()

    save_plot(img_figure, outname)
    plt.close(img_figure)


def plot_image(outname, image_file, region=None, cutting_region=None, vmin=None, vmax=None, colorbar='v'):
    """Generate plot from input image file

    Parameters:
    -----------
    image_files :str
        The path to the fits images file to be plotted together
    outputname : str
        The output name to be given to the produced image
    region : str
     region file indicating the regions to be overlaid on the image
    cutting_region : str
        region file with region to make a cut on the image
    vmin :float
        the minimum value for the colorbar scale

    """

    image = Image(image_file)
    img_figure, ax = plt.subplots(1, subplot_kw={'projection': image.wcs.wcs})

    if cutting_region is not None:
        mask = region_to_mask(image, cutting_region)
        image[np.where(mask)] = True
        image.crop()

    image.plot(ax=ax, scale='linear', colorbar=colorbar, show_xlabel=True, show_ylabel=True, zscale=True, vmin=vmin, vmax=vmax)

    ax.set_xlabel('Ra')
    ax.set_ylabel('Dec')

    if region is not None:
        plot_regions(region, ax)

    if colorbar is not None:
        im = ax.images
        cb = im[-1].colorbar
        cb.ax.set_ylabel(add_zlabel(image_file, image.unit))
    save_plot(img_figure, outname)
    plt.close(img_figure)


def tokovinin_model(seeing_ref, wavelengths, wave_front=22 * u.m):
    """Compute FWHM as a function of wavelength from Tokovinin, A. 2002, PASP, 114, 1156. Returns the of the FWHM in arcseconds.
    Parameters
    seeing_ref : the atmospheric seeing in arcseconds for lambda = 5000 A in angle units (i.e. FWHM)
    wave_front : the wavefront outer scale lenght. Default to the Paranal value Conan, R., Ziad, A., Borgnino, J., Martin, F., & Tokovinin, A. A. 2000,
    in Interferometry in Optical Astronomy, eds. P. Léna, & A. Quirrenbach,
    Proc. SPIE, 4006, 963
    wavelengths : in angstroms
    """
    prop_constant = seeing_ref * (5000 * u.angstrom) ** (0.2)
    epsilon_0 = prop_constant * wavelengths.to(u.angstrom) ** (-0.2)  # seeing expressed in FWHM

    r_0 = 0.976 * wavelengths / epsilon_0.to(u.rad).value  # Fried parameter

    return epsilon_0.to("arcsec").value * np.sqrt(1 - 2.183 * (r_0 / wave_front.to(u.angstrom)) ** 0.356)


def save_plot(figure, outputfile_name):
    figure.savefig("%s.png" % outputfile_name, bbox_inches='tight')
    print("Saved plot %s in %s" % (outputfile_name, os.getcwd()))


def plot_regions(region_file, ax, fits_header=None):
    """Plot regions onto the image from a ds9 file. If the fits file is given the package pyregion is used which will draw text as well.

    Parameteres
    -----------
    region_file: str
        ds9 region file in Image coordinates
    ax: matplotlib.ax
        ax where the regions are to be drawn
    fits_header: astropy.io.fits.Header
        Fits header to convert image coordinates to fk5 (this will draw texts too and it's preferable than ds9)
    """
    if fits is None:
        regs = read_ds9(region_file)
        for reg in regs:
            print(reg)
            print("WARNING: I believe this method no longer exist. See https://astropy-regions.readthedocs.io/en/latest/api/regions.Region.html#regions.Region")
            reg.plot(ax=ax)
    else:
        # takes a string as input!
        reg_string = open(region_file, "r").read()
        regs = pyregion.parse(reg_string).as_imagecoord(fits_header)
        print(regs)
        patch_list, artist_list = regs.get_mpl_patches_texts()

        for p in patch_list:
            ax.add_patch(p)

        for t in artist_list:
            ax.add_artist(t)


def format_axis(ax, labelsize=17):
    ax.minorticks_on()
    # major ticks
    ax.tick_params(length=10, width=2, direction='in', labelsize=labelsize)
    # minor ticks
    ax.tick_params(which='minor', length=5)


def create_linestyle_array(data_length):
    """Create an array of colors given the length of a dataset. Useful for plots where a unique color is needed for each dataset.

    The returned colors come from the jet map.

    Parameters
    ----------
    data_length : The length of your data for the color array creation.

    """
    print("Creating linestyle array for %i datasets" % data_length)
    m = np.array(['solid', 'dotted', 'dashed', 'dashdot'])

    while data_length > len(m):
        m = np.concatenate((m, m))
    return m


def create_color_array(data_length, cmap='jet'):
    """Create an array of colors given the length of a dataset. Useful for plots where a unique color is needed for each dataset.

    The returned colors come from the jet map.

    Parameters
    ----------
    data_length : The length of your data for the color array creation.

    """
    print("Creating color array for %i datasets" % data_length)
    setmap = plt.get_cmap(name=cmap)

    colors = setmap(np.linspace(0, 1, data_length))
    return colors


def get_markers_array(data_length):
    """Get an array of markers given the length of a dataset. Useful for plots where a unique marker is needed for each dataset.

    There are 17 different markers and after that they are repeated.

    Parameters
    ----------
    data_length : The length of your data for the marker array creation.

    """
    m = ['.', 'v', '*', 'x', '8', 's', 'p', 'P', '^', 'h', '+', 'X', 'D', 'd', '>', '<', 'o']

    while data_length > len(m):
        m.append(m)

    return m


def get_linestyle_array(data_length):
    """Get an array of markers given the length of a dataset. Useful for plots where a unique marker is needed for each dataset.

    There are 17 different markers and after that they are repeated.

    Parameters
    ----------
    data_length : The length of your data for the marker array creation.

    """
    m = ['-', '--', '-.', ':']

    while data_length > len(m):
        m.append(m)

    return m


def remove_legend_repetitions(ax):
    """Removes any repetead entries in the legend.
    Parameters
    ----------
    ax : The axis were the repeated entries are to be removed. """
    handles, labels = ax.get_legend_handles_labels()
    by_label = OrderedDict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='best', prop={'size': 14}, frameon=False)
