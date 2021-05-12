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

The bpt_geometry function creates a figure for each BPT region of each diagram, containing a close-up
of the region, colored according to the map to the right (masked to show only the corresponding region)
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

# sdir='/home/mparra/PARRA/Observ/Andres/optical/'
sdir='/home/mparra/PARRA/Observ/NGC5917/MUSE/'

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
            self.pure_starform_max=-0.032
            self.agnliner_inter = (-0.24, 0.5)
            self.int_inter = (-0.61, 1)
        elif self.index == 2:
            self.pure_starform_max=0.198
            self.agnliner_inter = (-0.22, 0.3)
            self.int_inter = (-1.1, 1.1)
        elif self.index == 3:
            self.pure_starform_max=-0.36
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
    img_figure, ax = plt.subplots(1, subplot_kw={'projection': image.wcs.wcs},figsize=(10,8))
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

    figure, ax = plt.subplots(1,figsize=(10,8))

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
    starform_regions = np.where(( logy < bpt_diagram.pure_starform_crv(logx)) &\
                                ( logx < bpt_diagram.pure_starform_max))
    bpt_data[starform_regions] = starform_region.index

    int_regions = np.where((logy > pure_starform) & (logx < int_curve) & (logx < bpt_diagram.pure_starform_max))
    int_regions=np.concatenate((int_regions,np.where((logx < int_curve) &\
                                                     (logx > bpt_diagram.pure_starform_max))),axis=1)
    int_regions=(np.array(int_regions[0]),np.array(int_regions[1]))
    bpt_data[int_regions] = int_region.index

    agn_regions = np.where((logx > int_curve) & (logy>agnliner_curve))
    
    bpt_data[agn_regions] = agn_region.index

    liner_regions=np.where((logx>int_curve) & (logy<agnliner_curve))
    bpt_data[liner_regions] = liner_region.index

    bpt_data[invalid_pixels] = np.nan

    #plotting the separations
    inter_starform = np.sort(np.reshape(logx, np.size(logx)))
    inter_starform=inter_starform[inter_starform<bpt_diagram.pure_starform_max]
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
                color=region.color, ls="None", marker="+", label=region.name)
    ax.legend()

    y_axis_map[0].data = bpt_data

    if 'WCSAXES' in y_axis_map[0].header:
        if y_axis_map[0].header["WCSAXES"] == 3:
            y_axis_map[0].header["WCSAXES"] = 2
            
    y_axis_map[0].header['COMMENT'] = "BPT diagram %d (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))" % bptype
    outfile = "%s/BPT_%d.fits" % (out, bptype)

    y_axis_map.writeto(outfile, overwrite=True)
    img = Image(outfile)
    # colormap with identical order
    regionslist=[starform_regions,int_regions,agn_regions,liner_regions]
    colorange=[starform_region.color,int_region.color,agn_region.color,liner_region.color]
    colorlist=[]
    for i in range(len(regionslist)):
        if np.size(regionslist[i])!=0:
            colorlist.append(colorange[i])
    #colormap with identical order, but only for non empty regions
    cmp=ListedColormap(colorlist)
    
    plot_bpt_map(img, cmp, "%s/bpt%d_image.png" % (out, bptype), regs, conts)

    figure.savefig("%s/bpt_%d.png" % (out, bptype))
    print("Results stored to %s" % outfile)

def bpt_geometry(map_1, map_2, out, bptype, regs=None,conts=None,color_geo='plasma'):

    '''computations - No change from bpt_single'''
    
    x_axis_map = fits.open(map_1)
    y_axis_map = fits.open(map_2)

    logx = np.log10(x_axis_map[0].data)
    logy = np.log10(y_axis_map[0].data)

    invalid_pixels = np.where((np.isnan(logx)) | (np.isnan(logy)))

    # 0 for composite or LINER
    bpt_data = np.zeros(y_axis_map[0].data.shape)

    # definition of the regions. The indexes are arbitrary
    starform_region = BPT_region(0, "blue", "Star Formation")
    int_region= BPT_region(1, "green", "Int.")
    agn_region = BPT_region(2, "purple", "AGN")
    liner_region = BPT_region(3, "orange", "LI(N)ER")

    regions=[starform_region,int_region,agn_region,liner_region]

    bpt_diagram = BPT_diagram(bptype)


    #defining each region
    pure_starform = bpt_diagram.pure_starform_crv(logx)
    int_curve = bpt_diagram.int_crv(logy)
    agnliner_curve = bpt_diagram.agnliner_crv(logx)
    starform_regions = np.where(( logy < bpt_diagram.pure_starform_crv(logx)) &\
                                ( logx < bpt_diagram.pure_starform_max))
    bpt_data[starform_regions] = starform_region.index

    int_regions = np.where((logy > pure_starform) & (logx < int_curve) & (logx < bpt_diagram.pure_starform_max))
    int_regions=np.concatenate((int_regions,np.where((logx < int_curve) &\
                                                     (logx > bpt_diagram.pure_starform_max))),axis=1)
    int_regions=(np.array(int_regions[0]),np.array(int_regions[1]))
    bpt_data[int_regions] = int_region.index

    agn_regions = np.where((logx > int_curve) & (logy>agnliner_curve))
    
    bpt_data[agn_regions] = agn_region.index

    liner_regions=np.where((logx>int_curve) & (logy<agnliner_curve))
    bpt_data[liner_regions] = liner_region.index

    bpt_data[invalid_pixels] = np.nan

    #computing the separations
    inter_starform = np.sort(np.reshape(logx, np.size(logx)))
    inter_starform=inter_starform[inter_starform<bpt_diagram.pure_starform_max]
    
    inter_def = logy[logy > bpt_diagram.int_inter[0]]
    inter_def = np.sort(inter_def[inter_def< bpt_diagram.int_inter[1]])

    agnliner_def = logx[logx > bpt_diagram.agnliner_inter[0]]
    agnliner_def = np.sort(agnliner_def[agnliner_def < bpt_diagram.agnliner_inter[1]])


    '''plots'''

    #we only plot a geometry graph if the corresponding region is not empty
    regslist=[starform_regions,int_regions,agn_regions,liner_regions]
    graph_list=[]
    for i in range(len(regslist)):
        if np.size(regslist[i])>=1:
            graph_list.append(i)
    
    figure_geo=np.array([None]*len(graph_list))
    ax_geo=np.array([[None]*2]*len(graph_list))
    
    for i in graph_list:
        
        #we first create an elevation map according to the distance to the central pixel, which will be the
        #basis for the colormap
        
        array_geo=np.copy(logx)
        cmap_geo=color_geo

        for j in range(np.size(logx,0)):
            for k in range(np.size(logx,1)):
                array_geo[j][k]=100-(abs(j-np.size(logx,0)/2)**2+abs(k-np.size(logx,0)/2)**2)**(1/2)
              
        #this 1d variable is for the scatter_plot
        array_geo_scat=array_geo[bpt_data == regions[i].index]
 
        #this 2d map is for the map plot
        array_geo_map=np.where(bpt_data==regions[i].index,array_geo,np.nan)
    
        #same process than plot_bpt_map here, but we don't use the function since we plot the figures
        #together
        y_axis_map[0].data = array_geo_map
        
        if 'WCSAXES' in y_axis_map[0].header:
            if y_axis_map[0].header["WCSAXES"] == 3:
                y_axis_map[0].header["WCSAXES"] = 2
                
        y_axis_map[0].header['COMMENT'] = "BPT diagram geometry %d (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))" % bptype
        outfile = "%s/BPT_%d_geo_%d.fits" % (out, bptype,i)
    
        y_axis_map.writeto(outfile, overwrite=True)

        img = Image(outfile)

        '''
        here we create the figure. Since it is not possible to create subplots with different subplot_kw
        keywords, we create both with the projection, delete the first and replace it with a non-projected
        subplot (note : it's not possible to do this the otherway because add_subplot doesn't take subplot_kw
        as an argument)
        '''
        
        figure_geo[i], ax_geo[i] = plt.subplots(1,2,figsize=(16,8),subplot_kw={'projection': img.wcs.wcs})
        
        ax_geo[i][0].remove()
        ax_geo[i][0]=figure_geo[i].add_subplot(1,2,1)
        
        figure_geo[i].suptitle('BPT geometry for the '+regions[i].name+' region in the BPT diagram '+str(bptype))
        
        logx_singreg=logx[bpt_data == regions[i].index]
        logy_singreg=logy[bpt_data == regions[i].index]
        
        #since the curves lead to wrong auto plot adjustments, we do it ourselves
        out_x=abs(np.nanmax(logx_singreg)-np.nanmin(logx_singreg)) * 0.05
        out_y=abs(np.nanmax(logy_singreg)-np.nanmin(logy_singreg)) * 0.05
        
        ax_geo[i][0].set_xlim(np.nanmin(logx_singreg)-out_x,
                              np.nanmax(logx_singreg)+out_x)
        ax_geo[i][0].set_ylim(np.nanmin(logy_singreg)-out_y,
                              np.nanmax(logy_singreg)+out_y)
        
        ax_geo[i][0].set_ylabel("log([OIII]/H$_\\beta$)")
        
        #definition of the curve parameters depending on the diagram type
        if bptype == 1:
    
            ax_geo[i][0].set_xlabel("log([NII]/H$_\\alpha$)")
    
        if bptype == 2:
    
            ax_geo[i][0].set_xlabel("log([SII]/H$_\\alpha$)")
    
        if bptype == 3:
    
            ax_geo[i][0].set_xlabel("log([OI]/H$_\\alpha$)")
        
        ax_geo[i][1].set_xlabel('Ra', labelpad=0)
        ax_geo[i][1].set_ylabel('Dec', labelpad=-2)
        
        #plots for the first subplot
        ax_geo[i][0].plot(inter_starform, bpt_diagram.pure_starform_crv(inter_starform), color="black", zorder=100)
        ax_geo[i][0].plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)
        ax_geo[i][0].plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    
        ax_geo[i][0].scatter(logx_singreg,logy_singreg,c=array_geo_scat, cmap=cmap_geo,ls="None", marker="+", label=regions[i].name)
        ax_geo[i][0].legend()  
        
        #plots for the second subplot
        img.plot(ax=ax_geo[i][1], scale='linear', show_xlabel=False, show_ylabel=False, cmap=cmap_geo, 
                   extent=None)
        
        if args.contours is not None:
            ctrs = fits.open(args.contours)
            min_data = np.nanpercentile(ctrs[0].data, 10)
            max_data = np.nanpercentile(ctrs[0].data, 90)
            levels = np.linspace(min_data, max_data, 4)
            cmp = ListedColormap(["black"])
            ax_geo[i][1].contour(ctrs[0].data, levels=levels, alpha=0.6, origin="lower", cmap=cmp)
    
        if regs is not None:
            for reg in regs:
                mu.plot_regions(reg, ax_geo[i][1], img.data_header)

        outfile_png = "%s/BPT_%d_geo_%d.png" % (out, bptype,i)
        
        figure_geo[i].savefig(outfile_png, format="png", bbox_inches="tight", pad_inches=0.4)
        
    
# read arguments
ap = argparse.ArgumentParser(description='Create BPT diagram from two given line ratio maps and the BPT diagram type (1, 2 or 3). The separation is based on Kewley et al. 2006')
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
        bpt_geometry(ratiomap,bptcfg.lineratio_paths["o3_hb"], outdir, i,args.regions, args.contours)
else:
    print("Running in auto mode.")
    
    os.chdir(sdir)

    r_results,r_names=ratio_maker(bpt_lines,'flux',outdir)

    if r_results[0]=='Unavailable' or r_results[1]==r_results[2]==r_results[3]=='Unavailable':
        print("Can't compute the BPT, too many ratios missing.")
    else:
        for i in range(1, 4):
            if r_results[i]=='Done':
                print('BPT graph '+str(i)+' can be plotted.')
                bpt_single(r_names[i],r_names[0], args.regions, args.contours, outdir, i)
                bpt_geometry(r_names[i],r_names[0], outdir, i,args.regions, args.contours)
