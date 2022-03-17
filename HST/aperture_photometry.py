# @Author: Andrés Gúrpide <agurpide>
# @Date:   19-10-2020
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 17-03-2022

import os
from regions import read_ds9
from photutils.aperture import aperture_photometry
from astropy.io import fits
import argparse
from numpy import log10, sqrt
import astropy.units as u
from astropy import wcs
import logging
import hst_utils as hst_ut


parser = argparse.ArgumentParser(description='Extracts fluxes from the given apertures.')
parser.add_argument("images", help="Image files where to extract fluxes", nargs='+', type=str)
parser.add_argument("-r", "--regions", type=str, help='Source (first) and background (second line, optional) extraction region file to use for the aperture photometry', nargs=1)
parser.add_argument("-a", "--aperture_correction", type=float, help='Aperture correction (see https://stsci.edu/hst/instrumentation/wfc3/data-analysis/photometric-calibration/uvis-encircled-energy and https://stsci.edu/hst/instrumentation/acs/data-analysis/aperture-corrections)', nargs="?", default=1)
args = parser.parse_args()
regions = read_ds9(args.regions[0])
source_reg = regions[0]
source_aperture = hst_ut.region_to_aperture(source_reg)
# if a background region was given
if len(regions) > 1:
    bkg_reg = regions[1]
    bkg_aperture = hst_ut.region_to_aperture(bkg_reg)
else:
    logging.warning("No background was given, no background correction will be performed.")
for image_file in args.images:

    if os.path.isfile(image_file):
        hst_hdul = fits.open(image_file)
        date = hst_hdul[0].header["DATE-OBS"]
        hst_filter = hst_ut.get_image_filter(hst_hdul[0].header)
        instrument = hst_hdul[0].header["INSTRUME"]
        units = hst_hdul[1].header["BUNIT"]
        exp_time = float(hst_hdul[0].header["EXPTIME"])
        detector = hst_hdul[0].header["DETECTOR"] if "DETECTOR" in hst_hdul[0].header else ""
        pivot_wavelength = float(hst_hdul[1].header["PHOTPLAM"])
        filter_bandwidth = float(hst_hdul[1].header["PHOTBW"])
        # if UV filter then https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2017/WFC3-2017-14.pdf
        # use phftlam1 keyword for UV filters
        uv_filters = ["F200LP", "F300X", "F218W", "F225W", "F275W", "FQ232N", "FQ243N", "F280N"]
        if detector == "UVIS" and filter in uv_filters:
            photflam = float(hst_hdul[0].header["PHTFLAM1"])
        elif "PHOTFLAM" in hst_hdul[0].header:
            photflam = float(hst_hdul[0].header["PHOTFLAM"])
        elif "PHOTFLAM" in hst_hdul[1].header:
            photflam = float(hst_hdul[1].header["PHOTFLAM"])
        photflam = photflam * u.erg / u.AA / u.s / u.cm**2
        print("PHOTFLAM keyword value: %.2E %s" % (photflam.value, photflam.unit))
        zero_point = float(hst_hdul[1].header["PHOTZPT"])
        image_data = hst_hdul[1].data
        hst_wcs = wcs.WCS(hst_hdul[1].header)
        phot_source = aperture_photometry(image_data, source_aperture, wcs=hst_wcs)
        if len(regions) > 1:
            phot_bkg = aperture_photometry(image_data, bkg_aperture, wcs=hst_wcs)

        aperture_keyword = "corrected_aperture_sum(%s)" % units
        # background correction
        if len(regions) > 1:
            phot_source[aperture_keyword] = phot_source["aperture_sum"] - phot_bkg["aperture_sum"] / bkg_aperture.area * source_aperture.area
            phot_source["corrected_aperture_err"] = sqrt(phot_source["aperture_sum"] + (sqrt(phot_bkg["aperture_sum"]) / bkg_aperture.area * source_aperture.area) ** 2)
        else:
            phot_source[aperture_keyword] = phot_source["aperture_sum"]
            phot_source["corrected_aperture_err"] = sqrt(phot_source["aperture_sum"])

        phot_source_conf = phot_source[aperture_keyword] + phot_source["corrected_aperture_err"]
        phot_source_conf_neg = phot_source[aperture_keyword] - phot_source["corrected_aperture_err"]
        counts_at_infinity = phot_source[aperture_keyword] / args.aperture_correction
        counts_at_infinity_conf = phot_source_conf / args.aperture_correction
        counts_at_infinity_neg = phot_source_conf_neg / args.aperture_correction

        # divide by the exposure time if needed
        if "/S" not in units:
            print("Units: %s. Applying exposure time correction" % units)
            counts_at_infinity /= exp_time
            counts_at_infinity_conf /= exp_time
            counts_at_infinity_neg /= exp_time

        flux_header = "flux(%s)" % photflam.unit
        phot_source[flux_header] = counts_at_infinity * photflam
        phot_source["flux_err"] = counts_at_infinity_conf * photflam - phot_source[flux_header]
        #https://www.stsci.edu/hst/instrumentation/acs/data-analysis/zeropoints
        phot_source["mag"] = -2.5 * log10(counts_at_infinity) - zero_point
        phot_source["mag_err_neg"] = -2.5 * log10(counts_at_infinity_conf) - zero_point - phot_source["mag"]
        phot_source["mag_err_pos"] = -2.5 * log10(counts_at_infinity_neg) - zero_point - phot_source["mag"]
        # formatting
        phot_source["xcenter"].info.format = '%.2f'
        phot_source["ycenter"].info.format = '%.2f'
        phot_source["aperture_sum"].info.format = '%.3f'
        phot_source[aperture_keyword].info.format = '%.2f'
        phot_source["corrected_aperture_err"].info.format = '%.1f'
        phot_source[flux_header].info.format = '%.3E'
        phot_source['flux_err'].info.format = '%.2E'
        phot_source["mag"].info.format = "%.2f"
        phot_source["mag_err_neg"].info.format = "%.2f"
        phot_source["mag_err_pos"].info.format = "%.3f"
        reg_basename = os.path.basename(args.regions[0]).replace('.reg', '')
        out_data_file = "aperture_phot_%s_%s.csv" % (hst_filter, reg_basename)
        phot_source.write(out_data_file, overwrite=True)
        out_info_file = image_file.replace(".fits", "%s_apt_info.txt" % reg_basename)
        f = open(out_info_file, "w+")
        f.write("Date:%s\nInstrument:%s\nFilter:%s\nDetector:%s\nExposure(s):%.2f\nPivot wavelength (A):%.1f\nRMS:%.1f\nPHOTFLAM:%s\nAperture correction:%.4f" % (date, instrument, hst_filter, detector, exp_time, pivot_wavelength, filter_bandwidth, photflam.value, args.aperture_correction))
        aperture_reg = type(source_aperture).__name__
        attributes = vars(source_aperture)
        att_list = ''.join("%s: %s" % item for item in attributes.items())
        f.write("\n###Aperture details###\n%s\n%s" % (aperture_reg, att_list))
        f.close()
        xspec_outfile = image_file.replace(".fits", "_to_xspec.ascii")
        f = open("%s" % ("%s" % xspec_outfile), "w+")
        f.write("%.2f %.2f %.3e %.3e" % (pivot_wavelength - filter_bandwidth / 2, pivot_wavelength + filter_bandwidth / 2, phot_source[flux_header].value, phot_source['flux_err'].value))
        f.close()
        print("Use 'ftflx2xsp infile=%s xunit=angstrom yunit=ergs/cm^2/s/A nspec=1 phafile=hst_%s.fits rspfile=hst_%s.rsp' to convert to XSPEC format" % (xspec_outfile, hst_filter, hst_filter))
        print("Warning: you may need to change the '-' sign in the ascii file for this to work")
        print("Output stored to %s and %s" % (out_info_file, out_data_file))
        #
        # print some useful info
