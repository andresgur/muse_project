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
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/BPT')
import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpdaf.obj import Image
import muse_utils as mu
from matplotlib.colors import ListedColormap
from line_utils import ratio_maker
import bpt_config as bptcfg
from astropy.wcs import WCS


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

        if self.index in [1,2,3]:
            self.y_axis="log([OIII]/H$_\\beta$)"
            if self.index == 1:
                self.x_axis="log([NII]/H$_\\alpha$)"
                self.agnliner_inter = (-0.24, 0.5)
                self.int_inter = (-0.61, 1)
            elif self.index == 2:
                self.x_axis="log([SII]/H$_\\alpha$)"
                self.agnliner_inter = (-0.22, 0.3)
                self.int_inter = (-1.1, 1)
            elif self.index == 3:
                self.x_axis="log([OI]/H$_\\alpha$)"
                self.agnliner_inter = (-0.9, 0.3)
                self.int_inter = (-0.25, 0.65)

        if self.index=='proj':
            self.x_axis="P1 (0.77*N2+0.54*S2+0.33*R3)"
            self.y_axis="P2 (-0.57*N2+0.82*S2)"

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

    def proj_x(self,logN2,logS2,logR3):
        return 0.77 * logN2 + 0.54 * logS2 + 0.33* logR3

    def proj_y(self,logN2,logS2):
        return -0.57*logN2+0.82*logS2

    def cold_projcrv(self,logy):
        return 2.597 * logy**3 - 1.865 * logy**2 + 0.105 * logy - 0.435

    def warm_projcrv(self,logy):
        return 3.4 * logy**3 - 2.233 * logy**2 - 0.184* logy - 0.172


def plot_bpt_map(image, colormap=ListedColormap(["blue", "green", "pink", "orange"]), outfile="bpt_image.png", regions=None, contours=None):
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
        # sometimes the data is in the 1 or the 0 hdu list
        file_data = ctrs[1].data if ctrs[0].data is None else ctrs[0].data
        min_data = np.nanpercentile(file_data, 30)
        max_data = np.nanpercentile(file_data, 87)
        levels = np.linspace(min_data, max_data, 3)
        print(levels)
        cmp = ListedColormap(["black"])
        cs = ax.contour(file_data, levels=levels, alpha=0.95, origin="lower", cmap=cmp, zorder=2, linestyles=[":", "--", "solid"])
        #ax.clabel(cs, cs.levels, inline=True, fmt="%d", fontsize=20)

    if regions is not None:
        for region in regions:
            mu.plot_regions(region, ax, image.data_header)
    img_figure.savefig(outfile, format="png", bbox_inches="tight", pad_inches=0.4)


def bpt_single(map_1, logy, regs, conts, bptype, colormap, grid_ax=None):

    x_axis_map = fits.open(map_1)

    logx = np.log10(x_axis_map[0].data)

    invalid_pixels = np.where((np.isnan(logx)) | (np.isnan(logy)))

    # 0 for composite or LINER
    bpt_data = np.zeros(x_axis_map[0].data.shape)
    bpt_indexes = np.zeros(x_axis_map[0].data.shape)

    figure, ax = plt.subplots(1)

    #since the curves lead to wrong auto plot adjustments, we do it ourselves
    out_x=abs(np.nanmax(logx)-np.nanmin(logx)) * 0.05
    out_y=abs(np.nanmax(logy)-np.nanmin(logy)) * 0.05
    ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
    if grid_ax is not None:
        grid_ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
        grid_ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    bpt_diagram = BPT_diagram(bptype)

    ax.set_xlabel(bpt_diagram.x_axis)
    ax.set_ylabel(bpt_diagram.y_axis)

    if grid_ax is not None:
        grid_ax.set_xlabel(bpt_diagram.x_axis)
    # defining each region
    pure_starform = bpt_diagram.pure_starform_crv(logx)
    int_curve = bpt_diagram.int_crv(logy)
    agnliner_curve = bpt_diagram.agnliner_crv(logx)
    starform_regions = np.where(logy < bpt_diagram.pure_starform_crv(logx))
    bpt_data[starform_regions] = logx[starform_regions] + logy[starform_regions]
    bpt_indexes[starform_regions] = 0
    int_regions = np.where((logy > pure_starform) & (logx < int_curve))
    bpt_data[int_regions] = logx[int_regions] + logy[int_regions]
    bpt_indexes[int_regions] = 1
    agn_regions = np.where((logx > int_curve) & (logy > agnliner_curve))
    bpt_data[agn_regions] = logx[agn_regions] + logy[agn_regions]
    bpt_indexes[agn_regions] = 2
    liner_regions = np.where((logx > int_curve) & (logy < agnliner_curve))
    bpt_data[liner_regions] = logx[liner_regions] + logy[liner_regions]
    bpt_indexes[liner_regions] = 3
    bpt_data[invalid_pixels] = np.nan
    bpt_indexes[invalid_pixels] = np.nan

    #plotting the separations
    inter_starform = np.sort(np.reshape(logx, np.size(logx)))
    ax.plot(inter_starform, bpt_diagram.pure_starform_crv(inter_starform), color="black", zorder=100)

    inter_def = logy[logy > bpt_diagram.int_inter[0]]
    inter_def = np.sort(inter_def[inter_def < bpt_diagram.int_inter[1]])
    ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)

    agnliner_def = logx[logx > bpt_diagram.agnliner_inter[0]]
    agnliner_def = np.sort(agnliner_def[agnliner_def < bpt_diagram.agnliner_inter[1]])
    # zorder 100 so it stands above the points
    ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    regions = [starform_regions, int_regions, agn_regions, liner_regions]

    if grid_ax is not None:
        grid_ax.plot(inter_starform, bpt_diagram.pure_starform_crv(inter_starform), color="black", zorder=100)
        grid_ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)
        grid_ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    # ploting the regions
    for map, region in zip(colormap, regions):
        if bpt_data[region].size == 0:
            continue
        cmap = mpl.cm.get_cmap(map)
        # apply offset to the minimum avoid very light colors
        norm = mpl.colors.Normalize(vmin=np.nanmin(bpt_data[region]) - 0.2, vmax=np.nanmax(bpt_data[region]))
        ax.scatter(logx[region], logy[region], c=bpt_data[region], cmap=cmap, norm=norm, ls="None", marker=".")
        if grid_ax is not None:
            grid_ax.scatter(logx[region], logy[region], c=bpt_data[region], cmap=cmap, norm=norm, ls="None", marker=".")
    ax.legend()
    bpt_fits = fits.PrimaryHDU(data=bpt_data, header=x_axis_map[0].header)
    if "WCSAXES" in bpt_fits.header:
        if bpt_fits.header["WCSAXES"] == 3:
            bpt_fits.header["WCSAXES"] = 2
    bpt_fits.header['COMMENT'] = "BPT diagram %d (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))" % bptype
    return bpt_fits, bpt_indexes, figure


def bpt_proj(map_N2, map_S2,map_R3, regs, conts, colormap, grid_ax=None):

    N2_map= fits.open(map_N2)
    S2_map= fits.open(map_S2)
    R3_map= fits.open(map_R3)

    logN2=np.log10(N2_map[0].data)
    logS2=np.log10(S2_map[0].data)
    logR3=np.log10(R3_map[0].data)

    invalid_pixels = np.where((np.isnan(logN2)) | (np.isnan(logS2)) | (np.isnan(logR3)))

    # 0 for composite or LINER
    bpt_data = np.zeros(N2_map[0].data.shape)
    bpt_indexes = np.zeros(N2_map[0].data.shape)

    figure, ax = plt.subplots(1)

    bpt_diagram = BPT_diagram('proj')

    #defining x and y for the projection according to the curves written in the class
    logx=bpt_diagram.proj_x(logN2, logS2, logR3)
    logy=bpt_diagram.proj_y(logN2, logS2)

    #since the curves lead to wrong auto plot adjustments, we do it ourselves
    out_x=abs(np.nanmax(logx)-np.nanmin(logx)) * 0.05
    out_y=abs(np.nanmax(logy)-np.nanmin(logy)) * 0.05

    ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)

    if grid_ax is not None:
        grid_ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
        grid_ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)

    ax.set_xlabel(bpt_diagram.x_axis)
    ax.set_ylabel(bpt_diagram.y_axis)

    if grid_ax is not None:
        grid_ax.set_xlabel(bpt_diagram.x_axis)

    # defining each region

    cold_regions=np.where(logx<bpt_diagram.cold_projcrv(logy))
    int_regions=np.where((logx>bpt_diagram.cold_projcrv(logy)) & (logx<bpt_diagram.warm_projcrv(logy)))
    warm_regions=np.where(logx>bpt_diagram.warm_projcrv(logy))

    bpt_data[cold_regions] = logx [cold_regions] + logy [cold_regions]
    bpt_indexes[cold_regions] = 0

    bpt_data[int_regions] = logx [int_regions] + logy [int_regions]
    bpt_indexes[int_regions] = 1

    bpt_data[warm_regions] = logx [warm_regions] + logy [warm_regions]
    bpt_indexes[warm_regions] = 2

    bpt_data[invalid_pixels] = np.nan
    bpt_indexes[invalid_pixels] = np.nan

    #plotting the separations
    logy_arr=np.sort(np.reshape(logy,np.size(logy)))

    #cold-int
    ax.plot(bpt_diagram.cold_projcrv(logy_arr),logy_arr,color="black",zorder=1000)

    #int-warm
    ax.plot(bpt_diagram.warm_projcrv(logy_arr),logy_arr,color="black",zorder=1000)

    regions = [cold_regions, int_regions, warm_regions]

    if grid_ax is not None:
        grid_ax.plot(bpt_diagram.cold_projecrv(logy_arr),logy_arr,color="black",zorder=1000)
        grid_ax.plot(bpt_diagram.warm_projecrv(logy_arr),logy_arr,color="black",zorder=1000)

    # ploting the regions
    for map, region in zip(colormap, regions):
        if bpt_data[region].size == 0:
            continue
        cmap = mpl.cm.get_cmap(map)
        # apply offset to avoid very white regions
        norm = mpl.colors.Normalize(vmin=np.nanmin(bpt_data[region]) - 0.2, vmax=np.nanmax(bpt_data[region]))
        ax.scatter(logx[region], logy[region], c=bpt_data[region], cmap=cmap, norm=norm, ls="None", marker=".")
        if grid_ax is not None:
            grid_ax.scatter(logx[region], logy[region], c=bpt_data[region], cmap=cmap, norm=norm, ls="None", marker=".")
    ax.legend()

    bpt_fits = fits.PrimaryHDU(data=bpt_data, header=N2_map[0].header)
    if "WCSAXES" in bpt_fits.header:
        if bpt_fits.header["WCSAXES"] == 3:
            bpt_fits.header["WCSAXES"] = 2
    bpt_fits.header['COMMENT'] = "2D projection of BPT diagrams N2 (log(NII/Ha) vs log(OIII/Hb)) and S2 (log(OIII/Hb) vs log(SII/Ha))"
    return bpt_fits, bpt_indexes, figure

# read arguments
ap = argparse.ArgumentParser(description='Create BPT diagram from two given line ratio maps and the BPT diagram type. Separations based on Law et al. 2021')
ap.add_argument("-r", "--regions", nargs='*', help="Region files to be overlaid on the image", default=None)
ap.add_argument("-c", "--contours", nargs='?', help="Fits file to use to overlay contours", type=str)
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='bpt_diagrams_v2', type=str)
ap.add_argument("--config", action="store_true", help="Flag to run from config file or from Maxime noob mode")
args = ap.parse_args()

# create output dir
outdir = args.outdir
if not os.path.isdir(outdir):
    os.mkdir(outdir)
#line names used for the automatic computation
#nested arrays of uneven dimension lengths and python are a pain so we create it carefully

bpt_lines=[[['OIII5007'],['HBETA']],
           [['NII6583'],['HALPHA']],
           [['SII6716','SII6731'],['HALPHA']],
           [['OI6300'],['HALPHA']]]

if args.config:
    # plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')
    bpt_type = bptcfg.type
    lineratiomaps = [bptcfg.lineratio_paths["n2_ha"], bptcfg.lineratio_paths["s2_ha"], bptcfg.lineratio_paths["oI_ha"]]

    '''STANDARD DIAGRAMS'''

    ymap = Image(bptcfg.lineratio_paths["o3_hb"])
    logy = np.log10(ymap.data.data)
    colormap = ["Blues", "Greens", "Purples", "Oranges"]

    grid_figure, grid_axes = plt.subplots(1, 3, sharey=True, gridspec_kw={"wspace": 0.05}, figsize=(19.2, 8))
    grid_img_figure, grid_im_axes = plt.subplots(1, 3, sharey=True, gridspec_kw={"wspace": 0.05}, subplot_kw={'projection': ymap.wcs.wcs}, figsize=(19.2, 8))
    grid_axes[0].set_ylabel("log([OIII]/H$_\\beta$)")

    grid_im_axes[0].set_ylabel('Dec', labelpad=-2)

    for myax in grid_im_axes[1:]:
        myax.coords[1].set_auto_axislabel(False)
        myax.coords[1].set_ticklabel_visible(False)

    for i, ratiomap in enumerate(lineratiomaps):
        bpt_type = i + 1
        bpt_fits, bpt_indexes, figure = bpt_single(ratiomap, logy, args.regions, args.contours, bpt_type, colormap, grid_axes[i])

        outfile = "%s/coloredBPT_%d.fits" % (outdir, bpt_type)
        bpt_fits.writeto(outfile, overwrite=True)
        bpt_img = Image(outfile)

        img_figure, ax = plt.subplots(1, subplot_kw={'projection': bpt_img.wcs.wcs})
        for index, map in enumerate(colormap):
            mask = (bpt_indexes != index) | (np.isnan(bpt_indexes))
            region = np.ma.masked_array(bpt_img.data, mask)
            if region.size == 0:
                continue
            #bpt_img.mask_selection(region)
            cmap = mpl.cm.get_cmap(map)
            norm = mpl.colors.Normalize(vmin=np.nanmin(region) - 0.2, vmax=np.nanmax(region))
            ax.imshow(region, cmap=cmap, norm=norm, origin="lower")
            grid_im_axes[i].imshow(region, cmap=cmap, norm=norm, origin="lower")
            bpt_img.unmask()
        grid_im_axes[i].set_xlabel('Ra', labelpad=1)
        ax.set_xlabel('Ra', labelpad=0)
        ax.set_ylabel('Dec', labelpad=-2)
        if args.contours is not None:
            ctrs = fits.open(args.contours)
            # sometimes the data is in the 1 or the 0 hdu list
            file_data = ctrs[1].data if ctrs[0].data is None else ctrs[0].data
            min_data = np.nanpercentile(file_data, 30)
            max_data = np.nanpercentile(file_data, 87)
            levels = np.linspace(min_data, max_data, 3)
            print(levels)
            cmp = ListedColormap(["black"])
            cs = ax.contour(file_data, levels=levels, alpha=0.95, origin="lower", cmap=cmp, zorder=2, linestyles=[":", "--", "solid"])
            #ax.clabel(cs, cs.levels, inline=True, fmt="%d", fontsize=20)

        if args.regions is not None:
            for region in args.regions:
                mu.plot_regions(region, ax, bpt_img.data_header)
                mu.plot_regions(region, grid_im_axes[i], bpt_img.data_header)
        img_figure.savefig(outfile.replace(".fits", "image.png"), format="png", bbox_inches="tight", pad_inches=0.4)
    grid_figure.savefig("%s/grid_bpts.png" % outdir)
    grid_img_figure.savefig("%s/grid_img_bpt.png" % outdir, format="png", bbox_inches="tight", pad_inches=0.4)
    print("Results stored to %s" % outfile)

    '''PROJECTION DIAGRAM'''

    ymap_proj = bptcfg.lineratio_paths["o3_hb"]
    colormap_proj = ["Blues", "Greens", "Reds"]

    bpt_type = 'proj'
    bpt_fits, bpt_indexes, figure = bpt_proj(lineratiomaps[0],lineratiomaps[1],ymap_proj, args.regions, args.contours, colormap_proj)

    outfile = "%s/coloredBPT_%s.fits" % (outdir, bpt_type)
    bpt_fits.writeto(outfile, overwrite=True)
    bpt_img = Image(outfile)

    img_figure_proj, ax_proj = plt.subplots(1, subplot_kw={'projection': bpt_img.wcs.wcs})
    for index, map in enumerate(colormap_proj):
        region = np.ma.masked_array(bpt_img.data, bpt_indexes != index)
        if region.size == 0:
            continue
        #bpt_img.mask_selection(region)
        cmap = mpl.cm.get_cmap(map)
        #norm = mpl.colors.Normalize(vmin=np.nanmin(region), vmax=np.nanmax(region))
        ax_proj.imshow(region, cmap=cmap, origin="lower")
        bpt_img.unmask()

    ax_proj.set_xlabel('Ra', labelpad=0)
    ax_proj.set_ylabel('Dec', labelpad=-2)
    if args.contours is not None:
        ctrs = fits.open(args.contours)
        # sometimes the data is in the 1 or the 0 hdu list
        file_data = ctrs[1].data if ctrs[0].data is None else ctrs[0].data
        min_data = np.nanpercentile(file_data, 30)
        max_data = np.nanpercentile(file_data, 87)
        levels = np.linspace(min_data, max_data, 3)
        print(levels)
        cmp = ListedColormap(["black"])
        cs = ax.contour(file_data, levels=levels, alpha=0.95, origin="lower", cmap=cmp, zorder=2, linestyles=[":", "--", "solid"])
        #ax.clabel(cs, cs.levels, inline=True, fmt="%d", fontsize=20)

    if args.regions is not None:
        for region in args.regions:
            mu.plot_regions(region, ax, bpt_img.data_header)

    img_figure_proj.savefig(outfile.replace(".fits", "image.png"), format="png", bbox_inches="tight", pad_inches=0.4)

    print("Results stored to %s" % outfile)

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
                bpt_single(r_names[i],r_names[0], args.regions, args.contours, i, colormap)
