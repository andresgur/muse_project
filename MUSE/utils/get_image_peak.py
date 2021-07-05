# @Author: Andrés Gúrpide <agurpide>
# @Date:   06-04-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 05-05-2021
# Script to derive the peak and the fwhm of an image and write the result to a region file

import muse_utils as mu
import argparse
from mpdaf.obj import Image
import numpy as np
import pyregion

ap = argparse.ArgumentParser(description='Determines the peak and FWHM of an image and write a region file with the result.')
ap.add_argument("input_file", nargs=1, help="Image or cube whose peak is to be determined")
ap.add_argument('-r', '--region', nargs='?', help="Region to crop the image", type=str, default=None)
ap.add_argument("--psf", help="Profile of the PSF fit", choices=['moffat', 'gaussian'], default="moffat")
args = ap.parse_args()
try:
    img = Image(args.input_file[0])
# error with the data extension not being the first one
except ValueError:
    img = Image(args.input_file[0], ext=1)
python_argument = "%s %s -r %s --psf %s" % (__file__, args.input_file[0], args.region, args.psf)
if args.region is not None:
    mask = mu.region_to_mask(img, args.region)
    img.mask[np.where(mask)] = True
    img.crop()
fwhm = img.fwhm()
center_region = pyregion.open(args.region)[0]
fit_center = (center_region.coord_list[1], center_region.coord_list[0])
try:

    if args.psf == "moffat":
        print("Fitting Moffat profile")
        profile_fit = img.moffat_fit(cont=1, fit_back=True, weight=True, full_output=False, flux=np.max(img.data), peak=True, circular=True,
                                     factor=1, center=fit_center, fwhm=fwhm, fit_n=True, n=1)

        rotation = profile_fit.rot
    elif args.psf == "gaussian":
        print("Fitting Gaussian profile")
        profile_fit = img.gauss_fit(center=fit_center, flux=np.max(img.data), fwhm=fwhm, circular=True,
                                    cont=0, fit_back=True, peak=True, factor=1, weight=True, verbose=True, full_output=False)
        rotation = 0
    peak = {"x": profile_fit.center[1], "y": profile_fit.center[0]}
    fwhm = profile_fit.fwhm

except ValueError:
    print("Couldn't fit the Moffat profile. FWHM and peak will be rough estimates")
    peak = img.peak()
    rotation = 0
print(peak)
print(fwhm)
out_region = "peak_region.reg"
f = open(out_region, "w+")
f.write("# Region file format: DS9 version 4.1\nglobal color=green dashlist=8 3 width=1 font='helvetica 10 normal roman' select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n")
f.write("fk5 \nellipse(%.6f, %.6f, %.4f\", %.4f\", %.2f) # color=white width=2 \n#python ~%s" % (peak["x"], peak["y"], fwhm[0], fwhm[1], rotation, python_argument))
f.close()
print("Results written to %s" % out_region)
