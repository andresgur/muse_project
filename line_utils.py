# @Author: Andrés Gúrpide <agurpide> and Maxime Parra
# @Date:   03-09-2019
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 26-04-2021
# Script to create line map ratios out of two line map images fits file

# imports
import os
import sys
import numpy as np
import logging
from astropy.io import fits
from mpdaf.obj import Image
import matplotlib.pyplot as plt
import copy as copy

def ratio_maker(line_names,datatype,outdir):
    
    '''
    automatically searches for linemaps in the subdirectories of the current directories,
    and if possible, computes the desired ratios in the "outdir" directory
    
    the shape of the "line_names" array to use is :
        line_names=[ratio_1,ratio_2,ratio_3,...]
    with : 
        ratio_x=[[num_1,num_2,...],[denum_1,denum_2,...]]
    all the elements must be standard camel line outputs (ex:'OII3727','HBETA',...)
    
    the datatypes can be any outputs of camel such as int, fux, wave,...
    if a single data type is given, it will be assumed for every ratio
    
    The function returns :
        -an array of length equal to the number of asked ratios, containing
        'Done' if the ratio could be computed, and 'Unavailable' if it couldn't
        -an array containing the newly created ratio file paths
    '''
    line_files=copy.deepcopy(line_names)
    ratio_results=np.array([None]*len(line_names))
    ratio_names=np.array([None]*len(line_names))
    if np.size(datatype)==1:
        datatype=np.array([datatype]*len(line_names))

    for i in range(len(line_names)):
        
        for j in [0,1]:
            for k in range(len(line_names[i][j])):
                    for root, dirs, files in os.walk(".", topdown=False):        
                        for name in files:
                            if name.endswith(line_names[i][j][k]+'.fits') and '_'+datatype[i]+'_' in name \
                                and 'cleaned_images' in root and 'smoothed2' not in root:
                                    line_files[i][j][k]=os.path.join(root,name)
        
        if any(j in line_files[i] for j in line_names[i]):
            ratio_results[i]='Unavailable'
        else:
            str_num='+'.join(line_names[i][0])                
            str_det='+'.join(line_names[i][1])
            ratio_names[i]=datatype[i]+'_'+str_num+'_'+str_det
            line_ratios(line_files[i][0],line_files[i][1],outdir,ratio_names[i])
            ratio_names[i]=outdir+'/'+ratio_names[i]+'.fits'
            ratio_results[i]='Done'
    return ratio_results,ratio_names
        
                                
  

def line_ratios(linemaps_numerator, linemaps_denominator,outdir,outname):
    
    #numerator : The line map(s) to be added in the numerator of the line ratio
    #denominator : The line map(s) to be added in the denominator of the line ratio
    #outdir : Output directory
    #outname : Output file name
    
    def add_maps(linemaps):
        """Add maps given a list of files
        Parameters
        ----------
        linemaps: list
            List of strings of maps to be added
        """
        added_maps_log = ' '.join(linemaps)
        added_data = []
        # add line maps
        logging.info('Adding %s maps together for numerator' % added_maps_log)
        
        
            
        for linemap in linemaps:
            if os.path.isfile(linemap):
                linemapfits = fits.open(linemap)
                if len(added_data) == 0:
                    added_data = linemapfits[0].data
                else:
                    added_data += linemapfits[0].data
            else:
                logging.warning("Line map %s not found." % linemap)
                continue
    
        return added_data, added_maps_log
    
    #avoid string parsing if there is only one line map in argument
    if type(linemaps_numerator)==str:
        linemaps_numerator=[linemaps_numerator]
    if type(linemaps_denominator)==str:
        linemaps_denominator=[linemaps_denominator]
    
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    numerator_data, added_maps_numerator = add_maps(linemaps_numerator)
    
    if len(numerator_data) == 0:
        logging.error("No line maps were found for the numerator")
        sys.exit()
    
    denominator_data, added_maps_denominator = add_maps(linemaps_denominator)
    
    if len(denominator_data) == 0:
        logging.error("No line maps were found for the denominator")
        sys.exit()
    # get the header from one of the files (or the next one in case it doesn't exist)
    if os.path.isfile(linemaps_denominator[0]):
        header = fits.open(linemaps_denominator[0])[0].header
    else:
        header = fits.open(linemaps_denominator[1])[0].header
    
    ratio_fits = fits.PrimaryHDU(data=numerator_data / denominator_data, header=header)
    ratio_fits.data[np.where(denominator_data is None)] = None
    
    ratio_fits.header['COMMENT'] = "Ratio of %s/%s line maps" % (added_maps_numerator, added_maps_denominator)
    if ratio_fits.header['WCSAXES'] == 3:
        ratio_fits.header['WCSAXES'] = 2  # set number of axes 2 for the image
    
    ratio_fits.writeto(outdir + "/" + outname + ".fits", overwrite=True)
    #print('Line ratio %s/%s written to %s/%s.fits' % (added_maps_numerator, added_maps_denominator, outdir, outname))
    
def map_plot(filename):    
    
    if os.path.exists(filename):
        
        img=Image(filename)
        img_figure, ax = plt.subplots(1, subplot_kw={'projection': img.wcs.wcs},figsize=(10,8))

        img.plot(ax=ax, scale='linear', show_xlabel=False, show_ylabel=False, zscale=True, colorbar='v',
                    extent=None)
        ax.set_xlabel('Ra', labelpad=0)
        ax.set_ylabel('Dec', labelpad=-2)
