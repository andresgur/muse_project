# @Author: Andrés Gúrpide <agurpide>
# @Date:   22-06-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 22-06-2021

import argparse
from regions import read_ds9
import numpy as np
from mpdaf.obj import Image
import muse_utils as mu
import pyregion

ap = argparse.ArgumentParser(description='Compute stats of input images over the whole FOV or certain regions.')
ap.add_argument("images", nargs='+', help='List of input images', type=str)
ap.add_argument("-r", "--regions", nargs='?', help='Region ds9 with one or more regions in it. Only annulus, circles and ellipses regions are valid', type=str)
args = ap.parse_args()
header = "#region\tmean\tstd\tmedian\tsum\tarea\n"
for image in args.images:
    img = Image(image)
    outstring = header
    if args.regions is not None:
        region_list = read_ds9(args.regions)
        pyregions = pyregion.open(args.regions)
        for region, (i, pyreg) in zip(region_list, enumerate(pyregions)):
            # mask everything outside the region
            sp = pyregion.core.ShapeList([pyreg])
            mask = sp.get_mask(img.get_data_hdu())
            img.mask = ~mask
            # mask nan values
            mask_sel = np.where(np.isnan(img.data))
            img.mask_selection(mask_sel)
            mean = np.ma.mean(img.data)
            median = np.ma.median(img.data)
            std = np.ma.std(img.data)
            sum = np.ma.sum(img.data)
            area = mu.region_to_aperture(region, img.wcs.wcs).area
            img.unmask()
            # use the name given to the region, otherwise number it
            if "label" in region.meta:
                region_name = region.meta["label"]
            else:
                region_name = "reg%d" % i
            outstring += "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" % (region_name, mean, std, median, sum, area)

    else:
        mean = np.ma.mean(img.data)
        median = np.ma.median(img.data)
        std = np.ma.std(img.data)
        sum = np.ma.sum(img.data)
        outstring += "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" % ("image", mean, std, median, sum, area)
    outfile = image.replace(".fits", ".dat")
    with open(outfile, "w+") as f:
        f.write(outstring)
    print("Output stored to %s" % outfile)
