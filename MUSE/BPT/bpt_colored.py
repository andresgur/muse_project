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
from bpt_utils import bpt_single,bpt_proj


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
        bpt_fits, bpt_indexes, figure = bpt_single(ratiomap, logy, args.regions, args.contours, bpt_type, colormap, grid_axes[i],figsize='big',bpt_labels=True)

        outfile = "%s/coloredBPT_%d.fits" % (outdir, bpt_type)
        
        figure.savefig(outfile.replace(".fits", ".png"), format="png", bbox_inches="tight", pad_inches=0.4)
        
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
        img_figure.savefig(outfile.replace(".fits", "image.png"), format="png", bbox_inches="tight", pad_inches=0.4, dpi=300)
        
    grid_axes[1].legend()
    grid_figure.savefig("%s/grid_bpts.png" % outdir)
    grid_img_figure.savefig("%s/grid_img_bpt.png" % outdir, format="png", bbox_inches="tight", pad_inches=0.4)
    print("Results stored to %s" % outfile)
    
    '''PROJECTION DIAGRAM'''
    
    ymap_proj = bptcfg.lineratio_paths["o3_hb"]
    colormap_proj = ["Blues", "Greens", "Reds"]

    bpt_type = 'proj'
    bpt_fits, bpt_indexes, figure = bpt_proj(lineratiomaps[0],lineratiomaps[1],ymap_proj, args.regions, args.contours, 
                                             colormap_proj,figsize='big',bpt_labels=True)

    outfile = "%s/coloredBPT_%s.fits" % (outdir, bpt_type)
    
    figure.savefig(outfile.replace(".fits", ".png"), format="png", bbox_inches="tight", pad_inches=0.4)
    
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
    
    img_figure_proj.savefig(outfile.replace(".fits", "image.png"), format="png", bbox_inches="tight")
    
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
