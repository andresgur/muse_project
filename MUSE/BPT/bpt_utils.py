# @Author: Andrés Gúrpide <agurpide> and Maxime Parra
# @Date:   07-07-2021
# @Email:  maxime.parra@irap.omp.eu
# !/usr/bin/env python3

import numpy as np
from astropy.io import fits
from mpdaf.obj import Image
import matplotlib.pyplot as plt
import matplotlib as mpl
import muse_utils as mu
from matplotlib.colors import ListedColormap
from mappings_utils import plot_shock_mod

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
            
            #region names
            self.region_names=["Star Formation","Int.","AGN","LI(N)ER"]
            
            #titles for graphs
            self.y_axis="log([OIII]/H$_\\beta$)"
            #column names for mappings tables
            self.y_column='OIII_Hb'
            
            if self.index == 1:
                self.x_axis="log([NII]/H$_\\alpha$)"
                self.x_column='NII_Ha'
                self.agnliner_inter = (-0.24, 0.5)
                self.int_inter = (-0.61, 1)
            elif self.index == 2:
                self.x_axis="log([SII]/H$_\\alpha$)"
                self.x_column='SII_Ha'
                self.agnliner_inter = (-0.22, 0.3)
                self.int_inter = (-1.1, 1)
            elif self.index == 3:
                self.x_axis="log([OI]/H$_\\alpha$)"
                self.x_column='OI_Ha'
                self.agnliner_inter = (-0.9, 0.3)
                self.int_inter = (-0.25, 0.65)
            
        if self.index=='proj':
            
            self.region_names=["dyn. cool","Int.","dyn. warm"]

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
    if contours is not None:
        ctrs = fits.open(contours)
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


def bpt_single(map_1, logy, regs, conts, bptype, colormap, grid_ax=None, title=None, shock_mods=None, mod_labels=None,figsize=None,bpt_labels=False):

    x_axis_map = fits.open(map_1)

    logx = np.log10(x_axis_map[0].data)

    invalid_pixels = np.where((np.isnan(logx)) | (np.isnan(logy)))

    # initialisating the maps
    bpt_data = np.zeros(x_axis_map[0].data.shape)
    bpt_indexes = np.zeros(x_axis_map[0].data.shape)

    #main bpt figure
    if figsize==None:
        figure, ax = plt.subplots(1)
    elif figsize=='big':
        figure, ax = plt.subplots(1,figsize=(10,8))
        
    #since the curves lead to wrong auto plot adjustments, we do it ourselves
    out_x=abs(np.nanmax(logx)-np.nanmin(logx)) * 0.05
    out_y=abs(np.nanmax(logy)-np.nanmin(logy)) * 0.05
    ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
    ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
    
    if title is not None:
        figure.suptitle(title)
    
    if grid_ax is not None:
        grid_ax.set_ylim(np.nanmin(logy)-out_y,np.nanmax(logy)+out_y)
        grid_ax.set_xlim(np.nanmin(logx)-out_x,np.nanmax(logx)+out_x)
        
    #introducing the BPT class for the axis names, limits, and curves
    bpt_diagram = BPT_diagram(bptype)
        
    # defining each region:
    #separations first
    pure_starform = bpt_diagram.pure_starform_crv(logx)
    int_curve = bpt_diagram.int_crv(logy)
    agnliner_curve = bpt_diagram.agnliner_crv(logx)
    
    #then the starforming region
    starform_regions = np.where(logy < bpt_diagram.pure_starform_crv(logx))
    bpt_data[starform_regions] = logx[starform_regions] + logy[starform_regions]
    bpt_indexes[starform_regions] = 0
    
    #the intermediaite region
    int_regions = np.where((logy > pure_starform) & (logx < int_curve))
    bpt_data[int_regions] = logx[int_regions] + logy[int_regions]
    bpt_indexes[int_regions] = 1
    
    #the agn region
    agn_regions = np.where((logx > int_curve) & (logy > agnliner_curve))
    bpt_data[agn_regions] = logx[agn_regions] + logy[agn_regions]
    bpt_indexes[agn_regions] = 2
    
    #and the liner region
    liner_regions = np.where((logx > int_curve) & (logy < agnliner_curve))
    bpt_data[liner_regions] = logx[liner_regions] + logy[liner_regions]
    bpt_indexes[liner_regions] = 3
    
    #the pixels of too low SNR stay at nan values to ensure white plotting in the maps
    bpt_data[invalid_pixels] = np.nan
    bpt_indexes[invalid_pixels] = np.nan

    #plotting the separations of the three curves
    #first one
    inter_starform = np.sort(np.reshape(logx, np.size(logx)))
    ax.plot(inter_starform, bpt_diagram.pure_starform_crv(inter_starform), color="black", zorder=100)

    #second one
    inter_def = logy[logy > bpt_diagram.int_inter[0]]
    inter_def = np.sort(inter_def[inter_def < bpt_diagram.int_inter[1]])
    ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)

    #third one
    agnliner_def = logx[logx > bpt_diagram.agnliner_inter[0]]
    agnliner_def = np.sort(agnliner_def[agnliner_def < bpt_diagram.agnliner_inter[1]])
    # zorder 100 so it stands above the points
    ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
    
    #array for the newly created regions
    regions = [starform_regions, int_regions, agn_regions, liner_regions]
    
    #mapping models plotting if inputted
    if shock_mods is not None:
        plot_shock_mod(shock_mods,ax,bpt_diagram,mod_labels)
            
    #setting the axis names for the main graphs (after the mapping models to superseed the axis names)
    ax.set_xlabel(bpt_diagram.x_axis)
    ax.set_ylabel(bpt_diagram.y_axis)
    
    #and other graphs if needed
    if grid_ax is not None:
        grid_ax.set_xlabel(bpt_diagram.x_axis)
        
    #plotting the separations in a supplementary graph if inputted
    if grid_ax is not None:
        grid_ax.plot(inter_starform, bpt_diagram.pure_starform_crv(inter_starform), color="black", zorder=100)
        grid_ax.plot(bpt_diagram.int_crv(inter_def), inter_def, color='black', zorder=100)
        grid_ax.plot(agnliner_def, bpt_diagram.agnliner_crv(agnliner_def), color='red', zorder=100)
        
    # ploting the regions in the maps
    for i, (map, region) in enumerate(zip(colormap, regions)):
        if bpt_data[region].size == 0:
            continue
        cmap = mpl.cm.get_cmap(map)
        norm = mpl.colors.Normalize(vmin=np.nanmin(bpt_data[region]) - 0.2, vmax=np.nanmax(bpt_data[region]))
        ax.scatter(logx[region], logy[region], c=bpt_data[region], cmap=cmap, norm=norm, ls="None", marker=".",
                   label=bpt_diagram.region_names[i] if bpt_labels==True else '')
        if grid_ax is not None:
            grid_ax.scatter(logx[region], logy[region], c=bpt_data[region], cmap=cmap, norm=norm, ls="None", marker=".",
                            label=bpt_diagram.region_names[i] if bpt_labels==True else '')
    ax.legend()
    bpt_fits = fits.PrimaryHDU(data=bpt_data, header=x_axis_map[0].header)
    
    if "WCSAXES" in bpt_fits.header:
        if bpt_fits.header["WCSAXES"] == 3:
            bpt_fits.header["WCSAXES"] = 2
    bpt_fits.header['COMMENT'] = "BPT diagram %d (1: log(OIII/Hb) vs log(NII/Ha), 2:log(OIII/Hb) vs log(SII/Ha), 3:log(OIII/Hb) vs log(OI/Ha))" % bptype
    return bpt_fits, bpt_indexes, figure

def bpt_proj(map_N2, map_S2,map_R3, regs, conts, colormap, grid_ax=None,figsize=None,title=None,bpt_labels=False):

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

    #main bpt figure
    if figsize==None:
        figure, ax = plt.subplots(1)
    elif figsize=='big':
        figure, ax = plt.subplots(1,figsize=(10,8))
    
    if title is not None:
        figure.suptitle(title)

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
    for i, (map, region) in enumerate(zip(colormap, regions)):
        if bpt_data[region].size == 0:
            continue
        cmap = mpl.cm.get_cmap(map)
        norm = mpl.colors.Normalize(vmin=np.nanmin(bpt_data[region]) - 0.2, vmax=np.nanmax(bpt_data[region]))
        ax.scatter(logx[region], logy[region], c=bpt_data[region], cmap=cmap, norm=norm, ls="None", marker=".",
                   label=bpt_diagram.region_names[i] if bpt_labels==True else '')
        if grid_ax is not None:
            grid_ax.scatter(logx[region], logy[region], c=bpt_data[region], cmap=cmap, norm=norm, ls="None", marker=".")
    ax.legend()
    
    bpt_fits = fits.PrimaryHDU(data=bpt_data, header=N2_map[0].header)
    if "WCSAXES" in bpt_fits.header:
        if bpt_fits.header["WCSAXES"] == 3:
            bpt_fits.header["WCSAXES"] = 2
    bpt_fits.header['COMMENT'] = "2D projection of BPT diagrams N2 (log(NII/Ha) vs log(OIII/Hb)) and S2 (log(OIII/Hb) vs log(SII/Ha))"
    return bpt_fits, bpt_indexes, figure

def bpt_geometry(map_1, map_2, out, bptype, regs=None,conts=None,color_geo='plasma'):

    '''computations - No change from bpt_single'''
    '''probably needs to be updated'''
    
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
        
        ax_geo[i][0].set_ylabel("log([OIII]/H$_\\beta$)",fontsize=12)
        
        #definition of the curve parameters depending on the diagram type
        if bptype == 1:
    
            ax_geo[i][0].set_xlabel("log([NII]/H$_\\alpha$)",fontsize=12)
    
        if bptype == 2:
    
            ax_geo[i][0].set_xlabel("log([SII]/H$_\\alpha$)",fontsize=12)
    
        if bptype == 3:
    
            ax_geo[i][0].set_xlabel("log([OI]/H$_\\alpha$)",fontsize=12)
        
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
        
        if conts is not None:
            ctrs = fits.open(conts)
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
