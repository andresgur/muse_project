#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 11:29:25 2021

@author: mparra
"""
#general stuff
import os,sys
import numpy as np
import logging
import pandas as pd

#astropy tools
from astropy import units as u
from astropy.wcs import WCS
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord,Angle,Distance

#regions and photometry
from regions import read_ds9,CircleSkyRegion
from photutils.detection import DAOStarFinder
from photutils.aperture import CircularAperture,aperture_photometry

#from local file
from apt_phot import allreg_to_aperture
from apt_phot import calculate_values

#plots
from mpdaf.obj import Image
import matplotlib.pyplot as plt


sdir='/home/mparra/PARRA/Observ/NGC3583/HST/batch/stacked'

name_555='clean_F555W.fits'
name_814='clean_F814W.fits'


os.chdir(sdir)
hdul_555=fits.open(name_555)
hdul_814=fits.open(name_814)


img_555=hdul_555[1].data[500:1500,500:1500]
img_814=hdul_814[1].data[500:1500,500:1500]

mean_555, median_555, std_555 = sigma_clipped_stats(img_555, sigma=3.0)
print('For 555, mean=%f, median=%f, std=%f'%(mean_555, median_555, std_555)) 
daofind_555 = DAOStarFinder(fwhm=2.0, threshold=4.*std_555)  
sources_555 = daofind_555(img_555 - median_555)  
positions_555 = np.transpose((sources_555['xcentroid'], sources_555['ycentroid']))
apt_555 = CircularAperture(positions_555, r=4.)

ch_ctr=["11:14:16.3078 +48:18:34.038"]

ch_sky=SkyCoord(ch_ctr,frame='fk5',unit=(u.hourangle, u.deg))

#angle for the region in which to restrain the detection
rad=Angle(10,'arcsec')

ch_region_sky=CircleSkyRegion(ch_sky[0],radius=rad)
WCS_555=WCS(hdul_555[1].header)
ch_region_555=ch_region_sky.to_pixel(WCS_555)

ch_ctr=np.array([ch_sky.dec.value[0],ch_sky.ra.value[0]])

ch_region_555.center.x=ch_region_555.center.x-500
ch_region_555.center.y=ch_region_555.center.y-500


img_figure_555, ax_555 = plt.subplots(1,2,subplot_kw={'projection': WCS_555},figsize=(16,8))
img_figure_555.suptitle('HST WFC3 F555W filter zoom around the ULX and broader view of NGC3583')
image_555=Image(name_555)
image_555=image_555[500:1500,500:1500]
image_555.plot(vmin=0,vmax=0.1,cmap='magma',scale='sqrt')
ch_region_555.plot(ax=ax_555[1],color='lime',label='10" region around the chandra source')

apt_555[0].plot(color='green',lw=1.5,alpha=0.5,label=r'F555W >4 s$\sigma$ sources')
apt_555[1:].plot(color='green',lw=1.5,alpha=0.5)
ax_555[1].legend(loc=3)

image_555_crop=image_555.copy()
image_555_crop.mask_ellipse(center=ch_ctr,radius=10,inside=False,
                            posangle=image_555_crop.get_rot())
image_555_crop.crop()
image_555_crop.plot(ax=ax_555[0],vmin=0,vmax=0.1,cmap='magma',scale='sqrt')
plt.show()

'''F814W detection and plots'''

#detection
mean_814, median_814, std_814 = sigma_clipped_stats(img_814, sigma=3.0)
print('For 814, mean=%f, median=%f, std=%f'%(mean_814, median_814, std_814)) 
daofind_814 = DAOStarFinder(fwhm=2.0, threshold=4.*std_814)  
sources_814 = daofind_814(img_814 - median_814)  
positions_814 = np.transpose((sources_814['xcentroid'], sources_814['ycentroid']))
apt_814 = CircularAperture(positions_814, r=4.)

#limitting region definition ; the WCS for the 814 filter is slightly different
WCS_814=WCS(hdul_814[1].header)
ch_region_814=ch_region_sky.to_pixel(WCS_814)
ch_region_814.center.x=ch_region_814.center.x-500
ch_region_814.center.y=ch_region_814.center.y-500

img_figure_814, ax_814 = plt.subplots(1,2,subplot_kw={'projection': WCS_814},figsize=(16,8))
img_figure_814.suptitle('HST WFC3 F814W filter zoom around the ULX and broader view of NGC3583')
image_814=Image(name_814)
image_814=image_814[500:1500,500:1500]
image_814.plot(vmin=0,vmax=0.1,cmap='magma',scale='sqrt')
ch_region_814.plot(ax=ax_814[1],color='lime',label='10" region around the chandra source')

apt_814[0].plot(color='red',lw=1.5,alpha=0.8,label=r'F814W >4 s$\sigma$ sources')
apt_814[1:].plot(color='red',lw=1.5,alpha=0.8)
ax_814[1].legend(loc=3)

image_814_crop=image_814.copy()
image_814_crop.mask_ellipse(center=ch_ctr,radius=10,inside=False,
                            posangle=image_814_crop.get_rot())
image_814_crop.crop()
image_814_crop.plot(ax=ax_814[0],vmin=0,vmax=0.1,cmap='magma',scale='sqrt')
plt.show()

'''
should be improved by switching pixel matching to sky position matching
'''

#maximal centroid position difference between two sources in both filters to be considered as a single source
n_pix=3

sources_combined=[]
for i in range(np.size(positions_555,0)):
    for j in range(np.size(positions_814,0)):
        if (positions_555[i][0]-positions_814[j][0])**2+(positions_555[i][1]-positions_814[j][1])**2<n_pix**2:
                #10" region test, we only do it for the F555W
                if (positions_555[i][0]-ch_region_555.center.x)**2+(positions_555[i][1]-ch_region_555.center.y)**2<ch_region_555.radius**2:
                    sources_combined.append([i,j])

#now we can still have sources that are associated with several others.
#for those, we should only consiert the closest every time
#for this specific dataset the problem doesn't appear

sources_combined=np.array(sources_combined)
trsrc=np.transpose(sources_combined)

if np.any(np.unique(trsrc[0])!=np.sort(trsrc[0])) or np.any(np.unique(trsrc[1])!=np.sort(trsrc[1])):
    logging.error('Sourcess overlap. Improve the code you lazy fuck.')
    sys.exit()
    
'''
we want to replot in the left graph the "cleaned" sources, however circular aperture doesn't work well with
subplots so we recreate the circles ourselves.
For this, we create an array containing the sources, and for each sources the corresponding circle for each filter
'''

circles_combined=np.copy(sources_combined)
circles_combined=circles_combined.astype(object)

#we use the chandra region pixel computation to know the difference in pixel between the non cropped and the
#cropped image. Translating the aperture positions then allows us to plot everything at the right place.
#The reason we don't use standard aperture plotting is because the function doesn't support subplots

crop_555_sx=ch_region_555.center.x-ch_region_555.radius
crop_555_sy=ch_region_555.center.y-ch_region_555.radius
crop_814_sx=ch_region_814.center.x-ch_region_814.radius
crop_814_sy=ch_region_814.center.y-ch_region_814.radius

for i in range(np.size(sources_combined,0)):
    circles_combined[i][0]=plt.Circle((positions_555[sources_combined[i][0]][0]-crop_555_sx,
                                       positions_555[sources_combined[i][0]][1]-crop_555_sy),4,color='green',fill=False,
                                      label=r'F555W >4 $\sigma$ sources with a counterpart' if i==0 else '')
    circles_combined[i][1]=plt.Circle((positions_814[sources_combined[i][1]][0]-crop_555_sx,
                                       positions_814[sources_combined[i][1]][1]-crop_555_sy),4,color='red',fill=False,
                                      label=r'F814W >4 $\sigma$ counterparts' if i==0 else '')
    
    #one artist can only be added to a single figure so we copy the circles to be able to draw them on 
    #both figures
    copy_1=plt.Circle((positions_555[sources_combined[i][0]][0]-crop_814_sx,
                                       positions_555[sources_combined[i][0]][1]-crop_814_sy),4,color='green',fill=False,
                                      label=r'F555W >4 $\sigma$ counterparts' if i==0 else '')
    copy_2=plt.Circle((positions_814[sources_combined[i][1]][0]-crop_814_sx,
                                        positions_814[sources_combined[i][1]][1]-crop_814_sy),4,color='red',fill=False,
                                      label=r'F814W >4 $\sigma$ sources with a counterpart' if i==0 else '')
    
    ax_555[0].add_patch(circles_combined[i][0])
    ax_555[0].add_patch(circles_combined[i][1])
    ax_814[0].add_patch(copy_1)
    ax_814[0].add_patch(copy_2)
    
    ax_555[0].set_title('10 arcsec zoom around the ULX position')
    ax_555[1].set_title('Broad view of the lower-left part of the galaxy')
    ax_814[0].set_title('10 arcsec zoom around the ULX position')
    ax_814[1].set_title('Broad view of the lower-left part of the galaxy')


#arbitrary common position, we do not forget to invert the abritrary coordinates
bg_ctr_x=609
bg_ctr_y=353
bg_apt=CircularAperture([bg_ctr_x,bg_ctr_y],r=8)

#circles for display, we don't care much about the difference in sky coordinates between both images here
bg_circ_555=plt.Circle((bg_ctr_x-crop_555_sx,bg_ctr_y-crop_555_sy),8,color='orange',fill=False,label='common background')
bg_circ_814=plt.Circle((bg_ctr_x-crop_814_sx,bg_ctr_y-crop_814_sy),8,color='orange',fill=False,label='common background')

ax_555[0].add_patch(bg_circ_555)
ax_814[0].add_patch(bg_circ_814)

#%% photometry for the detections
sys.path.append('/home/mparra/PARRA/Scripts/Python/HST/')


mjd_555=hdul_555[0].header['EXPSTART']
mjd_814=hdul_814[0].header['EXPSTART']

mjd=[mjd_555,mjd_814]
filters=['F555W','F814W']
detector='uvis1'

#array containing the photometry for each source, in each filter
#we directly compute pos and neg aperture elements, corresponding to -1 and +1 std errors
phot_comb=np.copy(circles_combined)

#common background photometry array (one element per image)
phot_bg=np.array([None]*2)

phot_bg[0]=aperture_photometry(img_555,bg_apt)
phot_bg[1]=aperture_photometry(img_814,bg_apt)
    
#zeropoints for an aperture of 0.396, i.e. 10 pixels. The last value is the vega magnitude
zpt=np.array([None]*2)
zpt[0]=calculate_values(detector,filters[0],mjd[0],'0.396')
zpt[1]=calculate_values(detector,filters[1],mjd[1],'0.396')

#vega magnitudes array
mv_values=np.array([[[None]*3]*2]*np.size(sources_combined,0))

for i in range(np.size(sources_combined,0)):

    phot_comb[i][0]=aperture_photometry(img_555,apt_555[sources_combined[i][0]])
    phot_comb[i][1]=aperture_photometry(img_814,apt_814[sources_combined[i][1]])
    
    for j in [0,1]:
    
        #we don't distinguish the apt_555 and 814 areas since we chose equal values
        phot_comb[i][j]["corrected_aperture"] = phot_comb[i][j]["aperture_sum"] - phot_bg[j]["aperture_sum"]  \
                                                / bg_apt.area * apt_555[i].area
        phot_comb[i][j]["corrected_aperture_err"] = np.sqrt(phot_comb[i][j]["aperture_sum"] + (np.sqrt(phot_bg[j]["aperture_sum"]) \
                                                    / bg_apt.area * apt_555[i].area) ** 2)
        phot_comb[i][j]["corrected_aperture_pos"] = phot_comb[i][j]["corrected_aperture"]  \
                                                    + phot_comb[i][j]["corrected_aperture_err"]
        phot_comb[i][j]["corrected_aperture_neg"] = phot_comb[i][j]["corrected_aperture"] \
                                                    - phot_comb[i][j]["corrected_aperture_err"]
        
        if phot_comb[i][j]["corrected_aperture_neg"][0]>0:
            mv_values[i][j][0]=-2.5*np.log10(phot_comb[i][j]["corrected_aperture_neg"][0]/0.396)+zpt[j][-1]
        if phot_comb[i][j]["corrected_aperture"][0]>0:
            mv_values[i][j][1]=-2.5*np.log10(phot_comb[i][j]["corrected_aperture"][0]/0.396)+zpt[j][-1]
        if phot_comb[i][j]["corrected_aperture_pos"][0]>0:
            mv_values[i][j][2]=-2.5*np.log10(phot_comb[i][j]["corrected_aperture_pos"][0]/0.396)+zpt[j][-1]



mv_noerr=np.transpose(mv_values)[1]

'''
Photometry for the ULX + position on the initial graph
'''

#the star region in sky coordinates works for both filters
path_reg_ulx='/home/mparra/PARRA/Observ/NGC3583/HST/batch/stacked/star_F814W.reg'

reg_ulx=read_ds9(path_reg_ulx)[0]

apt_ulx=allreg_to_aperture(reg_ulx)

#here, we only convert to pixel for the area so we don't care about the filter used
apt_ulx_pix_555=apt_ulx.to_pixel(WCS_555)

apt_ulx_pix_814=apt_ulx.to_pixel(WCS_814)
#here we once again use the full image to avoid problems do to sky coordinates
phot_ulx=np.array([None]*2)

mv_ulx=np.array([None]*2)


for i in [0,1]:
    
    phot_ulx[i]=aperture_photometry(hdul_555[1].data,apt_ulx,wcs=WCS_555)
    
    phot_ulx[i]["corrected_aperture"]=phot_ulx[i]["aperture_sum"] - phot_bg[i]["aperture_sum"]  \
                                                    / bg_apt.area * apt_ulx_pix_555.area
    
    mv_ulx[i]=-2.5*np.log10(phot_ulx[i]["corrected_aperture"][0]/0.396)+zpt[i][-1]
    
pos_ulx =np.array([[apt_ulx_pix_555.positions[0],apt_ulx_pix_555.positions[1]],
                   [apt_ulx_pix_814.positions[0],apt_ulx_pix_814.positions[1]]])

circle_ulx=np.array([plt.Circle((pos_ulx[0][0]-500-crop_555_sx,pos_ulx[0][1]-500-crop_555_sy),radius=apt_ulx_pix_555.r,
                                color='cyan',fill=False,linewidth=1.5,label='ULX counterpart'),
                     plt.Circle((pos_ulx[1][0]-500-crop_814_sx,pos_ulx[1][1]-500-crop_814_sy),radius=apt_ulx_pix_555.r,
                                color='cyan',fill=False,linewidth=1.5,label='ULX counterpart')])
ax_555[0].add_patch(circle_ulx[0])
ax_814[0].add_patch(circle_ulx[1])

ax_555[0].legend(loc=3)
ax_814[0].legend(loc=3)

#the corresponding detection in the automatic algorithm has this centroid in the F555W filter: 
pos_ulx_auto=[641.9762277573907,431.0926792004869]

#so we search for the corresponding index in the sources_combined array 
index_ulx_auto=np.where(sources_combined==np.where(positions_555==641.9762277573907)[0][0])[0][0]


#%%
'''isochrones importation for bigger ages'''

iso=np.array([None]*4)

zrange=[0.01,0.03,0.1,0.3]
zrange=[0.01,0.01,0.01,0.01]

path_iso='/home/mparra/PARRA/Observ/Padova/isochrones/ULX_1e7_1e8_001_030_WFC3.dat'
path_iso_lowZ='/home/mparra/PARRA/Observ/Padova/isochrones/ULX_1e7_1e8_001_030_WFC3_lowZ.dat'

iso_file_highZ=pd.read_csv(path_iso,skiprows=12,header=0,delimiter=' ',comment='#')
iso_file_lowZ=pd.read_csv(path_iso_lowZ,skiprows=12,header=0,delimiter=' ',comment='#')

iso_file=iso_file_lowZ.append(iso_file_highZ)

iso_Z=np.array(iso_file['Zini'])
iso_age=np.array(iso_file['logAge'])

iso_555mag=np.array(iso_file['F555Wmag'])
iso_814mag=np.array(iso_file['F814Wmag'])

#logarithmic values for ages (here it's 1-10*1e7 years)
un_ages=np.unique(iso_age)

#metallicity values
un_z=np.unique(iso_Z)

list_555mag=np.array([[None]*len(un_ages)]*len(un_z))
list_814mag=np.array([[None]*len(un_ages)]*len(un_z))

for i in range(len(un_z)):
    for j in range(len(un_ages)):
        list_555mag[i][j]=iso_555mag[np.intersect1d(np.where(iso_Z==un_z[i]),np.where(iso_age==un_ages[j]))]
        list_814mag[i][j]=iso_814mag[np.intersect1d(np.where(iso_Z==un_z[i]),np.where(iso_age==un_ages[j]))]
        
d_gal=Distance(29.8*u.Mpc)


#%% graphs

import numpy as np
import matplotlib.pyplot as plt

fig_iso=np.array([None]*3)
ax_iso=np.array([None]*3)

#indexes corresponding to a decent range of metallicities
plot_z=[0,2,9,11,18,38]

color_range = plt.cm.hsv(np.linspace(0.1,1,10))

for i in [0,1,2]:
    
    fig_iso[i],ax_iso[i]=plt.subplots(1,2,figsize=(16,8))
    fig_iso[i].suptitle('Color-Magnitude diagram for NGC5383 in the (F555W-F814W) versus (F814W) band')
    for j in [0,1]:
        #metallicity array index for the corresponding subplot
        k=plot_z[2*i+j]
        ax_iso[i][j].set_title('Padova isochrones with Z= '
                              +str(un_z[plot_z[2*i+j]])+r' Z$_\odot$')
        ax_iso[i][j].set_xlabel(r'$m_{555W}-m_{814W}$ (mag)')
        ax_iso[i][j].set_ylabel(r'$m_{555W}$ (mag)')
        ax_iso[i][j].set_xlim(-5,5)
        ax_iso[i][j].set_ylim(20,30)
        
        #plotting the star points
        ax_iso[i][j].scatter(mv_noerr[0]-mv_noerr[1],mv_noerr[1],marker='.',color='black',label='detected objects')
        
        #ulx detected counterpart
        ax_iso[i][j].scatter(mv_noerr[0][index_ulx_auto]-mv_noerr[1][index_ulx_auto],mv_noerr[1][index_ulx_auto],marker='+',color='red',
                             linewidth=2,label='ULX')
        
        # ax_iso[i][j].scatter(mv_ulx[0]-mv_ulx[1],mv_ulx[1],marker='+',color='red',linewidth=2,label='ULX')
        #plotting the isochrones
        color_range = plt.cm.hsv(np.linspace(0.1,1,10))
        for l in range(len(un_ages)):
            ax_iso[i][j].plot(list_555mag[k][l]-list_814mag[k][l],
                              list_814mag[k][l]+d_gal.distmod.value,color=color_range[l],
                              label=str((l+1)*10)+' Myr' if l//2==l/2 else '')
        
        #inverting the y axis
        ax_iso[i][j].set_ylim(ax_iso[i][j].get_ylim()[::-1])
        
        ax_iso[i][j].legend()
    fig_iso[i].tight_layout()
