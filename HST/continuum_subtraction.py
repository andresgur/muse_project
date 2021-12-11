# @Author: Andrés Gúrpide <agurpide>
# @Date:   26-04-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 11-12-2021
import argparse
from photutils.aperture import aperture_photometry
from regions import read_ds9
from astropy.io import fits
from astropy import wcs
import hst_utils as hst_ut
import numpy as np
def region_flux(hdu_image, apertures):
    image_wcs = wcs.WCS(hdu_image[1].header)

    hst_filter = hst_ut.get_image_filter(hdu_image[0].header)
    detector = hdu_image[0].header["DETECTOR"] if "DETECTOR" in hdu_image[0].header else ""
    uv_filters = ["F200LP", "F300X", "F218W", "F225W", "F275W", "FQ232N", "FQ243N", "F280N"]
    if detector == "UVIS" and hst_filter in uv_filters:
        photflam = float(hdu_image[0].header["PHTFLAM1"])
    elif "PHOTFLAM" in hdu_image[0].header:
        photflam = float(hdu_image[0].header["PHOTFLAM"])
    elif "PHOTFLAM" in hdu_image[1].header:
        photflam = float(hdu_image[1].header["PHOTFLAM"])
    counts = np.array([aperture_photometry(hdu_image[1].data, region_aperture, wcs=image_wcs)["aperture_sum"] for region_aperture in apertures])
    fluxes = counts
    return np.median(fluxes), hst_filter


# read arguments
ap = argparse.ArgumentParser(description='Produced a continuum subtracted narrow band image')
ap.add_argument("-n", "--narrow_band", nargs=1, help="Narrow band image", required=True)
ap.add_argument("-c", "--continuum_band", nargs=1, help="Continuum band image", required=True)
ap.add_argument("-r", "--region", nargs=1, help="Regions where to extract the count rates to rescale the images to take into account the different bandwidth. More than one will take the average.", required=True)
args = ap.parse_args()
regions = read_ds9(args.region[0])
apertures = [hst_ut.region_to_aperture(region) for region in regions]

narrow_band_image = args.narrow_band[0]
narrow_band_fits = fits.open(narrow_band_image)
flux_narrow, filter_narrow = region_flux(narrow_band_fits, apertures)
continuum_band_image = args.continuum_band[0]
continuum_band_fits = fits.open(continuum_band_image)
flux_continuum, continuum_filter = region_flux(continuum_band_fits, apertures)

rescale_factor = flux_narrow / flux_continuum

print("Applying a rescaling factor of %.2f" % rescale_factor)
clean_image = fits.PrimaryHDU(data=(narrow_band_fits[1].data - continuum_band_fits[1].data * rescale_factor), header=narrow_band_fits[1].header)
clean_image.header['COMMENT'] = "Image created after subtracting %s continuum image from %s narrow band image" % (continuum_band_image, narrow_band_image)
outfile = "%s_%s.fits" % (filter_narrow, continuum_filter)
clean_image.writeto("%s" % outfile, overwrite=True)
print("Wrote output image to %s" % outfile)
