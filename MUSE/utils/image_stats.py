# @Author: Andrés Gúrpide <agurpide>
# @Date:   22-06-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 10-02-2022

import argparse
from regions import Regions
import numpy as np
from mpdaf.obj import Image
import muse_utils as mu
import pyregion
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt


def get_stats(maskeddata):
    """Get the statistics of the masked data"""
    mean = np.ma.mean(maskeddata)
    median = np.ma.median(maskeddata)
    std = np.ma.std(maskeddata)
    validpixels = np.ma.count(maskeddata)

    std_err = std / (validpixels)**0.5
    sum = np.ma.sum(maskeddata)
    sum_2 = np.ma.sum(maskeddata ** 2)
    max = np.ma.max(maskeddata)
    min = np.ma.min(maskeddata)
    clipped_mean, clipped_median, clipped_std  = sigma_clipped_stats(maskeddata, sigma=10, maxiters=20)
    clipped_std_err = clipped_std / (validpixels)**0.5
    return mean, std, std_err, median, sum, sum_2, max, min, validpixels, clipped_mean, clipped_median, clipped_std, clipped_std_err

ap = argparse.ArgumentParser(description='Compute stats of input images over the whole FOV or certain regions.')
ap.add_argument("images", nargs='+', help='List of input images', type=str)
ap.add_argument("-r", "--region", nargs='?', help='ds9 region file with one or more regions in it. Only annulus, circles and ellipses regions are valid',
                type=str)
args = ap.parse_args()
header = "#region\tmean\tstd\tstd_err\tmedian\tsum\tsum_2\tmax\tmin\tclip_mean\tclip_median\tclip_std\tclip_std_err\tarea\tvalidpixels\n"
for image in args.images:
    try:
        img = Image(image)
    except ValueError:
        img = Image(image, ext=1)
    outstring = header
    if args.region is not None:
        region_list = Regions.read(args.region, format="ds9")
        pyregions = pyregion.open(args.region)
        for region, (i, pyreg) in zip(region_list, enumerate(pyregions)):
            # mask everything outside the region
            sp = pyregion.core.ShapeList([pyreg])
            mask = sp.get_mask(img.get_data_hdu())
            img.mask = ~mask
            # mask nan values
            mask_sel = np.where(np.isnan(img.data))
            img.mask_selection(mask_sel)
            mean, std, std_err, median, sum, sum_2, max, min, validpixels, clipped_mean, clipped_median, clipped_std, clipped_std_err = get_stats(img.data)
            try:
                area = mu.region_to_aperture(region, img.wcs.wcs).area
            except AttributeError:
                area = 0
            img.unmask()
            # use the name given to the region, otherwise number it
            if "label" in region.meta:
                region_name = region.meta["label"]
            else:
                region_name = "reg%d" % i
            outstring += "%s\t%.3f\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%d\n" % (region_name, mean, std, std_err, median, sum, sum_2, max, min, clipped_mean, clipped_median, clipped_std, clipped_std_err, area, validpixels)
            plt.hist(img.data[~img.mask], label=region_name, alpha=0.5, edgecolor='black')
            plt.axvline(np.percentile(img.data[~img.mask], 95), color='red', linestyle='--', label='95th percentile')
            plt.axvline(clipped_mean, color='green', linestyle='--', label='Clipped mean')
            
    else:
        mean, std, std_err, median, sum, sum_2, max, min, validpixels, clipped_mean, clipped_median, clipped_std, clipped_std_err = get_stats(img.data)

        outstring += "%s\t%.3f\t%.2f\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\t%.4f\t%.4f\t%.2f\t%d\t%d\n" % ("image", mean, std, std_err, median, sum, sum_2, max, min, clipped_mean, clipped_median, clipped_std, clipped_std_err, 0, validpixels)
        plt.hist(img.data[~img.mask], alpha=0.5, edgecolor='black')
        plt.axvline(np.percentile(img.data[~img.mask], 95), color='red', linestyle='--', label='95th percentile')
        plt.axvline(clipped_mean, color='green', linestyle='--', label='Clipped mean')
    outfile = image.replace(".fits", "_stats.dat")
    with open(outfile, "w+") as f:
        f.write(outstring)
    print("Output stored to %s" % outfile)
    plt.xlabel("Instances")
    plt.legend()
    plt.yscale("log")
    plt.savefig(image.replace(".fits", "_stats.png"))
