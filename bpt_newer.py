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
import sys, os
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

def bpt_single(map_1,map_2,regs,conts,out,bptype):
    
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
    out_x=abs(np.nanmax(logx)-np.nanmin(logx))*0.05
    out_y=abs(np.nanmax(logy)-np.nanmin(logy))*0.05
    ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
    plt.ylabel("log([OIII]/H$_\\beta$)")
    
    # definition of the regions. The indexes are arbitrary
    starform_region=BPT_region(0,"blue","Star Formation")
    int_region= BPT_region(1, "green", "Int.")
    agn_region = BPT_region(2, "purple", "AGN")
    liner_region = BPT_region(3, "orange", "LI(N)ER")
    
    regions=[starform_region,int_region,agn_region,liner_region]
    #colormap with identical order
    cmp=ListedColormap([starform_region.color,int_region.color,agn_region.color,liner_region.color])
    
    #definition of the curve parameters depending on the diagram type
    if bptype == 1:
        
        plt.xlabel("log([NII]/H$_\\alpha$)")
        #We define functions, variables and intervals. The function and intervals are used
        #for the separation plots, the variable for the regions
        def pure_starform_crv(x):
            return 0.359/(x+0.032)+1.083
        pure_starform=0.359/(logx+0.032)+1.083
        
        def int_crv(y):
            return -0.479*y**4-0.594*y**3-0.542*y**2-0.056*y-0.143
        int_curve=-0.479*logy**4-0.594*logy**3-0.542*logy**2-0.056*logy-0.143
        int_inter=[-0.61,1]
        
        def agnliner_crv(x):
            return 0.95*x+0.56
        agnliner_curve=0.95*logx+0.56
        agnliner_inter=[-0.24,0.5]

    if bptype == 2:
        
        plt.xlabel("log([SII]/H$_\\alpha$)")
        def pure_starform_crv(x):
            return 0.41/(x-0.198)+1.164
        pure_starform=0.41/(logx-0.198)+1.164
        
        def int_crv(y):
            return -0.943*y**4-0.45*y**3+0.408*y**2-0.61*y-0.025
        int_curve=-0.943*logy**4-0.45*logy**3+0.408*logy**2-0.61*logy-0.025
        int_inter=[-1.1,1]
        
        def agnliner_crv(x):
            return 1.89*x+0.76
        agnliner_curve=1.89*logx+0.76
        agnliner_inter=[-0.22,0.3]
        
    if bptype == 3:

        plt.xlabel("log([OI]/H$_\\alpha$)")
        def pure_starform_crv(x):
            return 0.612/(x+0.360)+1.179
        pure_starform=0.612/(logx+0.360)+1.179
        
        def int_crv(y):
            return 18.664*y**4-36.343*y**3+22.238*y**2-6.134*y-0.283
        int_curve=18.664*logy**4-36.343*logy**3+22.238*logy**2-6.134*logy-0.283
        int_inter=[-0.25,0.65]
        
        def agnliner_crv(x):
            return 1.18*x+1.3
        agnliner_curve=1.18*logx+1.3
        agnliner_inter=[-0.9,0.3]

        
    #defining each region
    starform_regions=np.where(logy<pure_starform)
    bpt_data[starform_regions] = starform_region.index
    
    int_regions=np.where((logy>pure_starform) & (logx<int_curve) \
                          & (logy>int_inter[0]) & (logy<int_inter[1]))
    bpt_data[int_regions]=int_region.index
    
    agn_regions=np.where((logx>int_curve) & (logy>agnliner_curve))
    bpt_data[agn_regions]=agn_region.index
    
    liner_regions=np.where((logx>int_curve) & (logy<agnliner_curve))
    bpt_data[liner_regions]=liner_region.index

    bpt_data[invalid_pixels] = np.nan
    
    #plotting the separations
    inter_starform=np.sort(np.reshape(logx,np.size(logx)))
    ax.plot(inter_starform, pure_starform_crv(inter_starform), color="black")
 
    inter_def=logy[logy>int_inter[0]]
    inter_def=np.sort(inter_def[inter_def<int_inter[1]])    
    ax.plot(int_crv(inter_def),inter_def,color='black')
    
    agnliner_def=logx[logx>agnliner_inter[0]]
    agnliner_def=np.sort(agnliner_def[agnliner_def<agnliner_inter[1]])
    ax.plot(agnliner_def,agnliner_crv(agnliner_def),color='red')
     

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
    plot_bpt_map(img, cmp, "%s/bpt%d_image.png" % (out, bptype), regs, conts)
    
    figure.savefig("%s/bpt_%d.png" % (out, bptype))
    print("Results stored to %s" % outfile)

# read arguments
ap = argparse.ArgumentParser(description='Create BPT diagram from two given line ratio maps and the BPT diagram type (1, 2 or 3). The separation is based on Kewley et al. 2006')
ap.add_argument("-map_y", nargs=1, help="Map with y-axis line ratio", type=str,default=None)
ap.add_argument("-map_x", nargs=1, help="Map with x-axis line ratio", type=str,default=None)
ap.add_argument("-r", "--regions", nargs='*', help="Region files to be overlaid on the image", default=None)
ap.add_argument("-c", "--contours", nargs='?', help="Fits file to use to overlay contours", type=str)
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='bpt_diagrams_v2', type=str)
ap.add_argument("--bpt", nargs="?", help="BPT diagram. Default 1. (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))", default=1, type=int)
args = ap.parse_args()

# plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

sdir='/home/mparra/PARRA/Observ/Andres/optical/'
outdir = args.outdir
    
#line names used for the automatic computation
#nested arrays of uneven dimension lengths and python are a pain so we create it carefully

bpt_lines=[[['OIII5007'],['HBETA']],
           [['NII6583'],['HALPHA']],
           [['SII6716','SII6731'],['HALPHA']],
           [['OI6300'],['HALPHA']]]

if args.map_y!=None and args.map_x!=None:
    
    bpt_single(args.map_x,args.max_y,args.regions,args.contours,outdir,args.bpt)
    
else :
    
    os.chdir(sdir)
    
    r_results,r_names=ratio_maker(bpt_lines,'flux',outdir)
    
    if r_results[0]=='Unavailable' or r_results[1]==r_results[2]==r_results[3]=='Unavailable':
        print("Can't compute the BPT, too many ratios missing.")
    else :
        for i in range(3):
            if r_results[i+1]=='Done':
                print('BPT graph '+str(i+1)+' can be plotted.')
                bpt_single(r_names[i+1],r_names[0],args.regions,args.contours,outdir,i+1)
