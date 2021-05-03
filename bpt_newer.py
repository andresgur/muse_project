# @Author: Andrés Gúrpide <agurpide> and Maxime Parra
# @Date:   29-04-2021
# @Email:  maxime.parra@irap.omp.eu
# !/usr/bin/env python3

'''
This script automatically searches for the BPT line maps in the subdirectories of
the starting directory, then computes the newest version of the BPT diagram from
Law et al. 2021 (https://arxiv.org/abs/2011.06012)

Arguments can be used to overlay contours, ds9 regions, change the output directory
and manually specify the starting line ratio files instead of doing an automatic search
'''

# imports
import sys
import os
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/')
sys.path.append('/home/mparra/PARRA/Scripts/Python/')
import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Image
import muse_utils as mu
from matplotlib.colors import ListedColormap
from line_utils import ratio_maker
import bpt_config as bptcfg


class BPT_region:
    """Wrapper for the regions of the BPT diagram.
    Parameters
    ----------
    index: int
        Index to keep track of each region (0, 1, 2 or 3)
    name: str
        Name of the region (for plotting purposes mainly i.e. HII, LINER, etc)
    color: str
        Color to assign to this region be used in plots
    """
    def __init__(self, index, color, name):
        self.index = index
        self.name = name
        self.color = color

class BPT_diagram:
    """Wrapper for the BPT diagram.
    Parameters
    ----------
    index: int
        Index to keep track of each region (0, 1, 2 or 3)
    name: str
        Name of the region (for plotting purposes mainly i.e. HII, LINER, etc)
    color: str
        Color to assign to this region be used in plots
    """
    def __init__(self, index):
        self.index = index

        if self.index == 1:
            self.agnliner_inter = (-0.24, 0.5)
            self.int_inter = (-0.61, 1)
        elif self.index == 2:
            self.agnliner_inter = (-0.22, 0.3)
            self.int_inter = (-1.1, 1)
        elif self.index == 3:
            self.agnliner_inter = (-0.9, 0.3)
            self.int_inter = (-0.25, 0.65)


    def pure_starform_crv(self, logx):
        if self.index == 1:
            return 0.359 / (logx + 0.032) + 1.083
        elif self.index == 2:
            return 0.41 / (logx - 0.198) + 1.164
        elif self.index == 3:
            return 0.612 / (logx + 0.360) + 1.179

    def int_crv(self, logy):
        if self.index == 1:
            return -0.479 * logy ** 4 - 0.594 * logy ** 3 - 0.542 * logy**2-0.056*logy-0.143
        elif self.index == 2:
            return -0.943*logy**4-0.45*logy**3+0.408*logy**2-0.61*logy-0.025
        elif self.index == 3:
            return 18.664*logy**4-36.343*logy**3+22.238*logy**2-6.134*logy-0.283

    def agnliner_crv(self, logx):
        if self.index == 1:
            return 0.95 * logx + 0.56
        elif self.index == 2:
            return 1.89 * logx + 0.76
        elif self.index == 3:
            return 1.18 * logx + 1.3


def plot_bpt_map(image, colormap=ListedColormap(["blue", "green","pink","orange"]), outfile="bpt_image.png", regions=None, contours=None):
    """ Create plot of the BPT diagram map
    Parameters
    ----------
    image: mpdaf.obj.Image
        Image file containing only 0s, 1s, 2s and 3s for each of the different BPT regions
    colormap:
        By default uses matplotlib.colors.ListedColormap(["blue", "green","pink","orange"])
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


def bpt_single(map_1, map_2, regs, conts, out, bptype):

    if not os.path.isdir(out):
        os.mkdir(outdir)
    x_axis_map = fits.open(map_1)
    y_axis_map = fits.open(map_2)

    logx = np.log10(x_axis_map[0].data)
    logy = np.log10(y_axis_map[0].data)

    invalid_pixels = np.where((np.isnan(logx)) | (np.isnan(logy)))

    # 0 for composite or LINER
    bpt_data = np.zeros(y_axis_map[0].data.shape)

    figure, ax = plt.subplots(1)

    #since the curves lead to wrong auto plot adjustments, we do it ourselves
    out_x=abs(np.nanmax(logx)-np.nanmin(logx)) * 0.05
    out_y=abs(np.nanmax(logy)-np.nanmin(logy)) * 0.05
    ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
    plt.ylabel("log([OIII]/H$_\\beta$)")

    # definition of the regions. The indexes are arbitrary
    starform_region = BPT_region(0, "blue", "Star Formation")
    int_region= BPT_region(1, "green", "Int.")
    agn_region = BPT_region(2, "purple", "AGN")
    liner_region = BPT_region(3, "orange", "LI(N)ER")

    regions=[starform_region,int_region,agn_region,liner_region]

    bpt_diagram = BPT_diagram(bptype)
    #definition of the curve parameters depending on the diagram type
    if bptype == 1:

        plt.xlabel("log([NII]/H$_\\alpha$)")

    if bptype == 2:

        plt.xlabel("log([SII]/H$_\\alpha$)")

    if bptype == 3:

        plt.xlabel("log([OI]/H$_\\alpha$)")

    #defining each region
    pure_starform = bpt_diagram.pure_starform_crv(logx)
    int_curve = bpt_diagram.int_crv(logy)
    agnliner_curve = bpt_diagram.agnliner_crv(logx)
    starform_regions = np.where(logy < bpt_diagram.pure_starform_crv(logx))
    bpt_data[starform_regions] = starform_region.index

    int_regions = np.where((logy > pure_starform) & (logx < int_curve) \
                          & (logy>bpt_diagram.int_inter[0]) & (logy<bpt_diagram.int_inter[1]))
    bpt_data[int_regions] = int_region.index

    agn_regions = np.where((logx > int_curve) & (logy>agnliner_curve))
    bpt_data[agn_regions] = agn_region.index

    liner_regions=np.where((logx>int_curve) & (logy<agnliner_curve))
    bpt_data[liner_regions] = liner_region.index

    bpt_data[invalid_pixels] = np.nan

    #plotting the separations
    inter_starform = np.sort(np.reshape(logx, np.size(logx)))
    ax.plot(inter_starform, bpt_diagram.pure_starform_crv(inter_starform), color="black", zorder=100)

    inter_def = logy[logy > bpt_diagram.int_inter[0]]
    inter_def = np.sort(inter_def[inter_def< bpt_diagram.int_inter[1]])
    ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)

    agnliner_def = logx[logx > bpt_diagram.agnliner_inter[0]]
    agnliner_def = np.sort(agnliner_def[agnliner_def < bpt_diagram.agnliner_inter[1]])
    # zorder 100 so it stands above the points
    ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)


    #ploting the regions
    for region in regions:
        ax.plot(logx[bpt_data == region.index], logy[bpt_data == region.index],\
                color=region.color, ls="None", marker=".", label=region.name)
    ax.legend()

    y_axis_map[0].data = bpt_data
    if y_axis_map[0].header["WCSAXES"] == 3:
        y_axis_map[0].header["WCSAXES"] = 2
    y_axis_map[0].header['COMMENT'] = "BPT diagram %d (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))" % bptype
    outfile = "%s/BPT_%d.fits" % (out, bptype)

    y_axis_map.writeto(outfile, overwrite=True)
    img = Image(outfile)
    
    # colormap with identical order, but only for non empty region
    regionslist=[starform_regions,int_regions,agn_regions,liner_regions]
    colorange=[starform_region.color,int_region.color,agn_region.color,liner_region.color]
    colorlist=[]
    for i in range(len(regionslist)):
        if np.size(regionslist[i])!=0:
            colorlist.append(colorange[i])
    cmp=ListedColormap(colorlist)
    
    plot_bpt_map(img, cmp, "%s/bpt%d_image.png" % (out, bptype), regs, conts)

    figure.savefig("%s/bpt_%d.png" % (out, bptype))
    print("Results stored to %s" % outfile)


# read arguments
ap = argparse.ArgumentParser(description='Create BPT diagram from two given line ratio maps and the BPT diagram type. Separations based on Law et al. 2021')
ap.add_argument("-r", "--regions", nargs='*', help="Region files to be overlaid on the image", default=None)
ap.add_argument("-c", "--contours", nargs='?', help="Fits file to use to overlay contours", type=str)
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='bpt_diagrams_v2', type=str)
ap.add_argument("--config", action="store_true", help="Flag to run from config file or from Maxime noob mode")
args = ap.parse_args()

outdir = args.outdir

#line names used for the automatic computation
#nested arrays of uneven dimension lengths and python are a pain so we create it carefully

bpt_lines=[[['OIII5007'],['HBETA']],
           [['NII6583'],['HALPHA']],
           [['SII6716','SII6731'],['HALPHA']],
           [['OI6300'],['HALPHA']]]

if args.config:
    plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')
    bpt_type = bptcfg.type
    lineratiomaps = [bptcfg.lineratio_paths["n2_ha"], bptcfg.lineratio_paths["s2_ha"], bptcfg.lineratio_paths["oI_ha"]]
    for i, ratiomap in zip(range(1, 4), lineratiomaps):

        bpt_single(ratiomap, bptcfg.lineratio_paths["o3_hb"], args.regions, args.contours, outdir, i)

else:
    print("Running in maxime noob mode let's goooo")
    sdir='/home/mparra/PARRA/Observ/Andres/optical/'
    os.chdir(sdir)

    r_results,r_names=ratio_maker(bpt_lines,'flux',outdir)

    if r_results[0]=='Unavailable' or r_results[1]==r_results[2]==r_results[3]=='Unavailable':
        print("Can't compute the BPT, too many ratios missing.")
    else:
        for i in range(1, 4):
            if r_results[i]=='Done':
                print('BPT graph '+str(i)+' can be plotted.')
                bpt_single(r_names[i],r_names[0], args.regions, args.contours, outdir, i)
