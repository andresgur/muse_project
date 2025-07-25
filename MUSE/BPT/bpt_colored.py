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
import argparse
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpdaf.obj import Image
import muse_utils as mu
from matplotlib.colors import ListedColormap
from matplotlib.patches import Circle
from astropy.wcs import WCS
from configparser import ConfigParser, ExtendedInterpolation
from bpt.bpt import BPT_1, BPT_2, BPT_3



class BPT_diagram:
    """Wrapper for the BPT diagram.
    Parameters
    ----------
    index: int
        Index to keep track of each region (1, 2 or 3)
    """
    def __init__(self, index):
        self.index = index

        self.y_axis="log([OIII]/H$_\\beta$)"

        if self.index=='proj':
            self.x_axis="P1 (0.77*N2+0.54*S2+0.33*R3)"
            self.y_axis="P2 (-0.57*N2+0.82*S2)"


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
    img_figure.savefig(outfile, bbox_inches="tight", pad_inches=0.4)


def bpt_single(map_1, logy, regs, conts, mplcolormaps, grid_ax=None, ext=0):

    """Create a BPT diagram from a single line ratio map and the y-axis values
    
    
    Returns
    --------
    array
        The BPT data color coded for the spatial map
    array
        The BPT indexes (0, 1, 2, 3)
    figure
        A matplotlib figure with the BPT diagram
    """

    x_axis_map = fits.open(map_1)
    logx = np.log10(x_axis_map[ext].data.data)
    invalid_pixels = np.where((np.isnan(logx)) | (np.isnan(logy)))

    #logx = logx[invalid_pixels]
    #logy = logy[invalid_pixels]

    # 0 for composite or LINER
    bpt_data = np.zeros(x_axis_map[ext].data.shape, dtype=float)
    bpt_indexes = np.zeros(x_axis_map[ext].data.shape) # cannot do int as we use np.nan later on

    figure, ax = plt.subplots(1)

    #since the curves lead to wrong auto plot adjustments, we do it ourselves
    offset = 0.04
    valid_indexes = ~((np.isnan(logx)) | (np.isnan(logy)) | (np.isinf(logx)) | (np.isinf(logy)))
    min_logx = np.nanmin(logx[valid_indexes])
    max_logx = np.nanmax(logx[valid_indexes])
    min_logy = np.nanmin(logy[valid_indexes])
    max_logy = np.nanmax(logy[valid_indexes])
    out_x=abs(max_logx-min_logx) * offset
    out_y=abs(max_logy-min_logy) * offset
    ax.set_xlim(min_logx-out_x,max_logx+out_x)
    ax.set_ylim(min_logy-out_y,max_logy+out_y)
    if grid_ax is not None:
        grid_ax.set_xlim(min_logx-out_x, max_logx+out_x)
        grid_ax.set_ylim(min_logy-out_y, max_logy+out_y)
    if bpt_type == 1:
        bpt_diagram = BPT_1()
    elif bpt_type==2:
        bpt_diagram = BPT_2()
    elif bpt_type==3:
        bpt_diagram = BPT_3()

    ax.set_xlabel(bpt_diagram.x_axis)
    ax.set_ylabel(bpt_diagram.y_axis)

    if grid_ax is not None:
        grid_ax.set_xlabel(bpt_diagram.x_axis)
    # defining each separation line
    pure_starform = bpt_diagram.pure_starform_crv(logx)
    int_curve = bpt_diagram.intermediate_crv(logy)
    agnliner_curve = bpt_diagram.agnliner_crv(logx)

    # classify

    bpt_data[invalid_pixels] = np.nan
    bpt_indexes[invalid_pixels] = np.nan
    # the star forming region diverges beyond limit, so we need to add this other condition
    int_regions = np.where(((logy > pure_starform) & (logx <= int_curve)) | (logx > bpt_diagram.limit))
    bpt_data[int_regions] = logx[int_regions] + logy[int_regions]
    bpt_indexes[int_regions] = 1
    # For logy < bpt_diagram.int_inter[1]), we check that we are above the intermediate-agn curve.
    # Elsewhere, this curve doesn't exist, so we check that we are above the starform curve.
    # In both cases, we check that we are above the LINER curve.
    agn_regions = np.where((((logx > int_curve) & (logy < bpt_diagram.int_inter[1])) | ((logy > pure_starform) & (logy > bpt_diagram.int_inter[1]))) & (logy > agnliner_curve))
    bpt_data[agn_regions] = logx[agn_regions] + logy[agn_regions]
    bpt_indexes[agn_regions] = 2
    # careful the int line does not work above a given limit
    liner_regions = np.where(((logy > pure_starform) | (logx > bpt_diagram.limit)) & (logy < agnliner_curve) & ((logx > int_curve) | (logy > bpt_diagram.int_inter[1])))
    bpt_data[liner_regions] = logx[liner_regions] + logy[liner_regions]
    bpt_indexes[liner_regions] = 3
    # the star forming region diverges beyond limit, so we need to add this other condition
    starform_regions = np.where((logy < pure_starform) & (logx < bpt_diagram.limit))
    bpt_data[starform_regions] = logx[starform_regions] + logy[starform_regions]
    bpt_indexes[starform_regions] = 0
    # plot the separations
    inter_starform = np.sort(np.reshape(logx, np.size(logx)))
    # avoid divergence and
    inter_good_values = inter_starform[np.where(inter_starform < bpt_diagram.limit)]
    ax.plot(inter_good_values, bpt_diagram.pure_starform_crv(inter_good_values), color="black", zorder=100)

    # set some limits over the lines so they look beautiful
    inter_def = logy[logy > bpt_diagram.int_inter[0]]
    inter_def = np.sort(inter_def[inter_def < bpt_diagram.int_inter[1]])
    ax.plot(bpt_diagram.intermediate_crv(inter_def), inter_def, color='black', zorder=100)

    agnliner_def = logx[logx > bpt_diagram.agnliner_inter[0]]
    agnliner_def = np.sort(agnliner_def[agnliner_def < bpt_diagram.agnliner_inter[1]])
    # zorder 100 so it stands above the points
    ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)

    if grid_ax is not None:
        grid_ax.plot(inter_good_values, bpt_diagram.pure_starform_crv(inter_good_values), color="black", zorder=100)
        grid_ax.plot(bpt_diagram.intermediate_crv(inter_def), inter_def, color='black', zorder=100)
        grid_ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    regions = [starform_regions, int_regions, agn_regions, liner_regions]
    # ploting the regions
    for cmap, region_index in zip(mplcolormaps, regions):
        # masking the nan and infinite values
        mask = ((np.isnan(bpt_data[region_index])) | (np.isinf(bpt_data[region_index])))
        region = np.ma.masked_array(bpt_data[region_index], mask)
        # getting a list of only unmasked values
        region_c = region.compressed()
        if region_c.size == 0:
            # no values for this region
            continue
        # apply offset to the minimum avoid very light colors
        norm = mpl.colors.Normalize(vmin=np.nanmin(region_c, 0) - 0.2,
                                     vmax=np.nanpercentile(region_c, 99))
        ax.scatter(logx[region_index], logy[region_index], c=region, cmap=cmap, norm=norm, ls="None", marker=".")
        if grid_ax is not None:
            grid_ax.scatter(logx[region_index], logy[region_index], c=region, cmap=cmap, norm=norm, ls="None", marker=".")
    ax.legend()

    # save fits files
    # colored maps
    return bpt_data, bpt_indexes, figure


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
    for map, region_index in zip(colormap, regions):
        # masking the nan and infinite values
        mask = ((np.isnan(bpt_data[region_index])) | (np.isinf(bpt_data[region_index])))
        region = np.ma.masked_array(bpt_data[region_index], mask)
        # getting a list of only unmasked values
        region_c = region.compressed()
        if region_c.size == 0:
            # no values for this region
            continue
        cmap = mpl.colormaps[map]
        # apply offset to the minimum avoid very light colors
        norm = mpl.colors.Normalize(vmin=np.nanmin(region_c, 0) - 0.2,
                                     vmax=np.nanpercentile(region_c, 99))
        ax.scatter(logx[region_index], logy[region_index], c=region, cmap=cmap, norm=norm, ls="None", marker=".")
        if grid_ax is not None:
            grid_ax.scatter(logx[region_index], logy[region_index], c=region, cmap=cmap, norm=norm, ls="None", marker=".")
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
ap.add_argument("--config", nargs=1, help="Input config file", required=True)
args = ap.parse_args()

# create output dir
outdir = args.outdir
if not os.path.isdir(outdir):
    os.mkdir(outdir)

# define colors for each BPT region
colormap = ["Blues", "Greens", "Purples", "Oranges"]
mplcolormaps = [mpl.colormaps[map] for map in colormap]


stylefile = "%s/.config/matplotlib/stylelib/paper.mplstyle" % os.environ["HOME"]
if os.path.isfile(stylefile):
    plt.style.use(stylefile)

# read config file
bptcfg = ConfigParser(interpolation=ExtendedInterpolation())
bptcfg.read("%s" %args.config[0])

bpt_type = bptcfg["bpt"]["diagram"]
lineratio_paths = bptcfg["Paths"]
lineratiomaps = [lineratio_paths["n2_ha"], lineratio_paths["s2_ha"], lineratio_paths["oI_ha"]]


'''STANDARD DIAGRAMS'''

ymap = Image(lineratio_paths["o3_hb"], ext=1)
print("WARNING: values below -1.0 for  [OIII]/Hb are being ignored")
ymap.data[np.where(np.log10(ymap.data) < -1.0)] = np.nan
logy = np.log10(ymap.data.data)
print(ymap.wcs.wcs)
grid_figure, grid_axes = plt.subplots(1, 3, sharey=True, gridspec_kw={"wspace": 0.05}, figsize=(25.6, 9.8))
grid_img_figure, grid_im_axes = plt.subplots(1, 3, sharey=True, gridspec_kw={"wspace": 0.05}, subplot_kw={'projection': ymap.wcs.wcs}, figsize=(25.6, 9.8))
grid_axes[0].set_ylabel("log([OIII]/H$_\\beta$)")

grid_im_axes[0].set_ylabel('Dec', labelpad=-2)

if bptcfg.has_section("ULX_position"):
    ra=float(bptcfg["ULX_position"]["ra"])
    dec=float(bptcfg["ULX_position"]["dec"])
    err=float(bptcfg["ULX_position"]["err"])/3600
    circle=(ra, dec, err)

    for i in range(len(grid_im_axes)):
        c = Circle(circle[:2],circle[2],
                    edgecolor='lime',
                    transform=grid_im_axes[i].get_transform('fk5'),
                    facecolor='none', lw=5)
        grid_im_axes[i].add_patch(c)

for myax in grid_im_axes[1:]:
    myax.coords[1].set_auto_axislabel(False)
    myax.coords[1].set_ticklabel_visible(False)

for i, ratiomap in enumerate(lineratiomaps):
    bpt_type = i + 1
    bpt_map, bpt_indexes, figure = bpt_single(ratiomap, logy, args.regions,
                                               args.contours, mplcolormaps, grid_axes[i], ext=1)
    

    outfile = "%s/coloredBPT_%d.fits" % (outdir, bpt_type)
    # map
    data = bpt_map
    bpt_fits = fits.PrimaryHDU(data=data, header=ymap.primary_header)
    if "WCSAXES" in bpt_fits.header:
        if bpt_fits.header["WCSAXES"] == 3:
            bpt_fits.header["WCSAXES"] = 2
    bpt_fits.header['COMMENT'] = "BPT diagram %d (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))" % bpt_type
    bpt_fits.writeto(outfile, overwrite=True)

    # single indexes maps
    outfile = "%s/BPT_%d.fits" % (outdir, bpt_type)
    data = bpt_indexes
    bpt_fits = fits.PrimaryHDU(data=data, header=ymap.primary_header)
    if "WCSAXES" in bpt_fits.header:
        if bpt_fits.header["WCSAXES"] == 3:
            bpt_fits.header["WCSAXES"] = 2
    bpt_fits.header['COMMENT'] = "BPT diagram %d (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))" % bpt_type
    bpt_fits.writeto(outfile, overwrite=True)
    # save figure
    figure.savefig("%s/coloredBPT_%d.png" % (outdir, bpt_type))
    plt.close(figure)

    img_figure, ax = plt.subplots(1, subplot_kw={'projection': ymap.wcs.wcs})

    if bptcfg.has_section("ULX_position"):
        c = Circle(circle[:2],circle[2],
                    edgecolor='lime',
                    transform=ax.get_transform('fk5'),
                    facecolor='none', lw=1.5)
        ax.add_patch(c)

    for index, cmap in enumerate(mplcolormaps):
        mask = ((bpt_indexes != index) | (np.isnan(bpt_indexes)))
        region = np.ma.masked_array(bpt_map.data, mask)
        mask_inf = ((np.isnan(region)) | (np.isinf(region)))
        region = np.ma.masked_array(region, mask_inf)
        # getting a list of only unmasked values
        region_c = region.compressed()
        if region_c.size == 0:
            # no values for this region
            continue
        # apply offset to the minimum avoid very light colors
        norm = mpl.colors.Normalize(vmin=np.nanmin(region_c, 0) - 0.2,
                                              vmax=np.nanpercentile(region_c, 99))
        ax.imshow(region, cmap=cmap, norm=norm, origin="lower", interpolation="nearest")
        grid_im_axes[i].imshow(region, cmap=cmap, norm=norm, origin="lower", interpolation="nearest")


    grid_im_axes[i].set_xlabel('Ra', labelpad=1)
    ax.set_xlabel('Ra', labelpad=1)
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
            mu.plot_regions(region, ax, ymap.primary_header)
            mu.plot_regions(region, grid_im_axes[i], ymap.primary_header)
    img_figure.savefig(outfile.replace(".fits", "image.png"), format="png", bbox_inches="tight", pad_inches=0.4, dpi=200)
    print("Results stored to %s" % outfile)
grid_figure.savefig("%s/grid_bpts.png" % outdir)
grid_img_figure.savefig("%s/grid_img_bpt.png" % outdir, bbox_inches="tight", pad_inches=0.4, dpi=200)

'''PROJECTION DIAGRAM'''
if False:
    ymap_proj = lineratio_paths["o3_hb"]
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
        cmap = mpl.colormarp[map]
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

    img_figure_proj.savefig(outfile.replace(".fits", "image.png"), bbox_inches="tight", pad_inches=0.4)

print("Results stored to %s" % outfile)
