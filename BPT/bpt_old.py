# @Author: Andrés Gúrpide <agurpide>
# @Date:   01-03-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 28-04-2021
# Script to create BPT diagram from two given line ratio maps
# !/usr/bin/env python3

# imports
import argparse
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Image
import muse_utils as mu
from matplotlib.colors import ListedColormap


class BPT_region:
    """Wrapper for the regions of the BPT diagram.

    Parameters
    ----------
    index: int
        Index to keep track of each region (0, 1 or 2)
    name: str
        Name of the region (for plotting purposes mainly i.e. HII, LINER, etc)
    color: str
        Color to assign to this region be used in plots
    """
    def __init__(self, index, color, name):
        self.index = index
        self.name = name
        self.color = color


def plot_bpt_map(image, colormap=ListedColormap(["green", "blue", "orange"]), outfile="bpt_image.png", regions=None, contours=None):
    """ Create plot of the BPT diagram map

    Parameters
    ----------
    image: mpdaf.obj.Image
        Image file containing only 0s, 1s and 2s for each of the different BPT regions
    colormap:
        By default uses matplotlib.colors.ListedColormap(["green", "blue", "orange"])
    outfile: str
        Name of the output file
    region: str
        Region file to be drawn (can contain multiple regions)
    contours: str
        Path to fits file to overlay contours on the BPT diagrams
    """
    img_figure, ax = plt.subplots(1, subplot_kw={'projection': image.wcs.wcs})
    # do not draw the colorbar colorbar=None
    # extent = left, right, bottom, top
    image.plot(ax=ax, scale='linear', show_xlabel=False, show_ylabel=False, zscale=False, cmap=colormap, colorbar=None,
               extent=None)
    ax.set_xlabel('Ra', labelpad=0)
    ax.set_ylabel('Dec', labelpad=-2)
    if args.contours is not None:
        ctrs = fits.open(args.contours)
        min_data = np.nanpercentile(ctrs[0].data, 10)
        max_data = np.nanpercentile(ctrs[0].data, 90)
        levels = np.linspace(min_data, max_data, 4)
        cmp = ListedColormap(["black"])
        ax.contour(ctrs[0].data, levels=levels, alpha=0.6, origin="lower", cmap=cmp)

    if regions is not None:
        for region in regions:
            mu.plot_regions(region, ax, image.data_header)
    img_figure.savefig(outfile, format="png", bbox_inches="tight", pad_inches=0.4)


# read arguments
ap = argparse.ArgumentParser(description='Create BPT diagram from two given line ratio maps and the BPT diagram type (1, 2 or 3). The separation is based on Kewley et al. 2006')
ap.add_argument("-map_y", nargs=1, help="Map with y-axis line ratio", type=str, required=True)
ap.add_argument("-map_x", nargs=1, help="Map with x-axis line ratio", type=str, required=True)
ap.add_argument("-r", "--regions", nargs='*', help="Region files to be overlaid on the image", default=None)
ap.add_argument("-c", "--contours", nargs='?', help="Fits file to use to overlay contours", type=str)
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='bpt_diagrams', type=str)
ap.add_argument("--bpt", nargs="?", help="BPT diagram. Default 1. (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))", default=1, type=int)
args = ap.parse_args()
plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')
outdir = args.outdir
if not os.path.isdir(outdir):
    os.mkdir(outdir)
y_axis_map = fits.open(args.map_y[0])
x_axis_map = fits.open(args.map_x[0])

logy = np.log10(y_axis_map[0].data)
logx = np.log10(x_axis_map[0].data)
invalid_pixels = np.where((np.isnan(logx)) | (np.isnan(logy)))
# 0 for composite or LINER
bpt_data = np.zeros(y_axis_map[0].data.shape)

figure, ax = plt.subplots(1)
plt.ylabel("log([OIII]/H$_\\beta$)")
if args.bpt == 1:
    # the indexes are arbitrary
    agn_region = BPT_region(0, "orange", "AGN")
    composite_region = BPT_region(1, "blue", "Composite")
    star_forming_region = BPT_region(2, "green", "HII")
    regions = [star_forming_region, agn_region, composite_region]
    # same order as for the regions
    cmp = ListedColormap([agn_region.color, composite_region.color, star_forming_region.color])
    # star forming region Equation 1 and 4 in Kewley 2006
    pure_star_forming = 0.61 / (logx - 0.05) + 1.3
    star_forming_regions = np.where(logy < pure_star_forming)
    pure_agn = 0.61 / (logx - 0.47) + 1.19
    # Composite regions are above pure star forming and below AGN
    composite_regions = np.where((logy < pure_agn) & (logy > pure_star_forming))
    bpt_data[composite_regions] = composite_region.index
    ax.plot(logx, pure_agn, color="red")
    ax.plot(logx, pure_star_forming, color="black", ls="--")
    plt.xlabel("log([NII]/H$_\\alpha$)")
elif args.bpt == 2 or args.bpt == 3:
    liner_region = BPT_region(0, "yellow", "LINER")
    star_forming_region = BPT_region(1, "green", "HII")
    seyfert_region = BPT_region(2, "purple", "Seyfert")
    regions = [star_forming_region, seyfert_region, liner_region]
    cmp = ListedColormap([liner_region.color, star_forming_region.color, seyfert_region.color])
    if args.bpt == 2:
        # star forming region
        pure_star_forming = 0.72 / (logx - 0.32) + 1.3
        star_forming_regions = np.where(logy < pure_star_forming)
        # seyfert - Line line
        def seyfert_liner_line(logx): return 1.89 * logx + 0.76
        seyfert_regions = np.where((logy > pure_star_forming) & (logy > seyfert_liner_line(logx)))
        ax.plot(logx, pure_star_forming, color="red", zorder=1000)
        ax.plot(logx[logx > -0.3], seyfert_liner_line(logx[logx > -0.3]), color="black", ls="--", zorder=1000)
        plt.xlabel("log([SII]/H$_\\alpha$)")
    elif args.bpt == 3:
        # star forming region
        def pure_star_forming(logx): return 0.73 / (logx + 0.59) + 1.33
        star_forming_regions = np.where((logy < pure_star_forming(logx)) & (logx < -0.59))
        # Seyfert regions
        def seyfert_liner_line(logx): return 1.18 * logx + 1.3
        seyfert_regions = np.where((pure_star_forming(logx) < logy) & (seyfert_liner_line(logx) < logy))
        ax.plot(np.sort(logx[(logx < -0.85)]), pure_star_forming(np.sort(logx[(logx < -0.85)])), color="red", zorder=1000)
        ax.plot(np.sort(logx[logx > -1.12]), seyfert_liner_line(np.sort(logx[logx > -1.12])), color="black", ls="--", zorder=1000)
        plt.xlabel("log([OI]/H$_\\alpha$)")
    # Seyfert regions
    bpt_data[seyfert_regions] = seyfert_region.index

bpt_data[star_forming_regions] = star_forming_region.index
bpt_data[invalid_pixels] = np.nan
for region in regions:
    ax.plot(logx[bpt_data == region.index], logy[bpt_data == region.index], color=region.color, ls="None", marker=".", label=region.name)

y_axis_map[0].data = bpt_data
if y_axis_map[0].header["WCSAXES"] == 3:
    y_axis_map[0].header["WCSAXES"] = 2
y_axis_map[0].header['COMMENT'] = "BPT diagram %d (1:star forming, 2: AGN, 3 composite)" % args.bpt
outfile = "%s/BPT_%d.fits" % (outdir, args.bpt)
ax.legend()
y_axis_map.writeto(outfile, overwrite=True)
img = Image(outfile)
plot_bpt_map(img, cmp, "%s/bpt%d_image.png" % (outdir, args.bpt), args.regions, args.contours)

figure.savefig("%s/bpt_%d.png" % (outdir, args.bpt))
print("Results stored to %s" % outfile)
