#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Author: Andrés Gúrpide <agurpide>
# @Date:   06-04-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 12-12-2024

import muse_utils as mu
import argparse
from mpdaf.obj import Image
import numpy as np
import pyregion
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser(description='Determines the peak and FWHM of a source or an image and write a region file with the result.')
ap.add_argument("input_file", nargs=1, help="Image or cube whose peak is to be determined", type=str)
ap.add_argument('-r', '--region', nargs='?', help="Region to crop the image", type=str, default=None)
ap.add_argument("--psf", help="Profile of the PSF fit. [moffat, gaussian] Default moffat.", choices=['moffat', 'gaussian'],
                default="moffat")
ap.add_argument("--circular", help="Flag for circular profile. Default False", action='store_true')
ap.add_argument("-b", "--background", nargs="?", help="The background value, otherwise it will be fitted.",
                type=float)
ap.add_argument("--ext", help="Extension to read from the fits file. Default will try 0",
                default=0, type=int)
args = ap.parse_args()
img = Image(args.input_file[0], ext=args.ext)
python_argument = "%s %s -r %s --psf %s --ext %d" % (__file__, args.input_file[0], args.region,
                                            args.psf, args.ext)
if args.region is not None:
    mask = mu.region_to_mask(img, args.region)
    img.mask[:] = mask
    img.crop()
    img.plot(zscale=True, use_wcs=True, scale="log")
    plt.savefig("cropped_image.png", dpi=100)

center_region = pyregion.open(args.region)[0]
fit_center = (center_region.coord_list[1], center_region.coord_list[0])

print("Circular profile: %s" % args.circular)

fwhm = img.fwhm_gauss()
print("fwhm", fwhm)
try:
    if args.psf == "moffat":
        print("Fitting Moffat profile")
        profile_fit = img.moffat_fit(cont=0, fit_back=True, weight=True, full_output=False, flux=np.max(img.data), peak=True, circular=args.circular,
                                     factor=1, center=fit_center, fwhm=fwhm, fit_n=False, n=2)

        rotation = profile_fit.rot
    elif args.psf == "gaussian":
        print("Fitting Gaussian profile")
        if args.background is None:
            fit_background = True
            cont = 0
        else:
            fit_background = False
            cont = args.background
        profile_fit = img.gauss_fit(center=fit_center, flux=np.max(img.data), fwhm=fwhm, circular=args.circular,
                                    cont=cont, fit_back=fit_background, peak=True, factor=1, weight=True, verbose=True, full_output=True)
        rotation = profile_fit.rot
    peak = {"x": profile_fit.center[1], "y": profile_fit.center[0]}
    fwhm = profile_fit.fwhm

except ValueError:
    print("Couldn't fit the %s profile. FWHM and peak will be rough estimates" % args.psf)
    peak = img.peak()
    rotation = 0
print("Peak:")
print(peak)
print("FWHM:")
print(fwhm)
print("Rotation:")
print(rotation)
if img.primary_header["TELESCOP"] == "2.5m":
    rotation += img.primary_header["SPA"]
elif "ESO-VLT" in img.primary_header["TELESCOP"]:
    print("WARNING THE ROTATION VALUE in the final ds9 region MIGHT BE INVALID FOR ESO TELESCOPES!")
out_region = "%s_fit_%s.reg" % (args.psf, args.region.replace(".reg",""))
f = open(out_region, "w+")
f.write("# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
if args.circular:
    f.write("fk5 \ncircle(%.6f, %.6f, %.4f\") # color=green width=2 \n#python %s" % (peak["x"], peak["y"], fwhm[1], python_argument))
else:
    f.write("fk5 \nellipse(%.6f, %.6f, %.4f\", %.4f\", %.2f) # color=green width=2 \n#python %s" % (peak["x"], peak["y"], fwhm[1], fwhm[0], rotation, python_argument))
if args.psf == "moffat":
    f.write("\n#n = %.1f" %profile_fit.n)
f.write("\n#cont = %.5f" %profile_fit.cont)
f.write("\n#(RA_err, DEC_err) = (%.5f, %.5f) arcsec" %(profile_fit.err_center[1] *3600, profile_fit.err_center[0] * 3600))
f.write("\n#(FWHM_err, FWHM_err) = (%.3f, %.3f) arcsec" %(profile_fit.err_fwhm[1], profile_fit.err_fwhm[0]))

f.close()
print("Results written to %s" % out_region)
print("Use 'ds9 %s -region %s' to check the results" %(args.input_file[0], out_region))
