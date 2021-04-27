#%%H2 metallicity estimation : initialisation

import sys, os
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/')

import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from mpdaf import *
from mpdaf.obj import Image, WCS,dms2deg,deg2sexa

from line_ratios import line_ratios
import numpy as np
import os

#This part computes metallicity maps from the Pilyugin and Grebel 2016 (doi:10.1093/mnras/stw238)
#Valid only for HII regions
#There are three methods depending on the amount of lines available

sdir='/home/mparra/PARRA/Observ/Andres/optical/linemaps'
os.chdir(sdir)

#those computations are based on 4 line INTENSITY ratios.
#We will first try to find them in the subdirectories of the current directory
files_R2=np.array(['None']*3,dtype='object')
files_N2=np.array(['None']*3,dtype='object')
files_S2=np.array(['None']*3,dtype='object')
files_R3=np.array(['None']*3,dtype='object')


line_names=['HBETA','OII3727','OII3729','NII6548','NII6583','SII6716','SII6731','OIII4959','OIII5007']

for root, dirs, files in os.walk(".", topdown=False):        
        for name in files:
            for i in range(len(line_names)):
                if name.endswith(line_names[i]+'.fits') and '_int' in name and 'smoothed2' not in root \
                    and 'cleaned_images' in root:
                    if i==0:
                        files_R2[-1]=os.path.join(root,name)
                        files_N2[-1]=os.path.join(root,name)
                        files_S2[-1]=os.path.join(root,name)
                        files_R3[-1]=os.path.join(root,name)
                        print('\nHBeta line found \n')
                    else:
                        if i<3:
                            files_R2[i-1]=os.path.join(root,name)
                        elif i<5:
                            files_N2[i-3]=os.path.join(root,name)
                        elif i<7:
                            files_S2[i-5]=os.path.join(root,name)
                        elif i<9:
                            files_R3[i-7]=os.path.join(root,name)

#computing the ratios when they are available

if 'None' not in files_R2:
    print('OII/HBeta ratio usable. Computation...\n')
    line_ratios(files_R2[:-1],[files_R2[-1]],'H2_metal','R2_ratio')
    ratio_R2=fits.open('./H2_metal/R2_ratio.fits')
    ratio_R2=ratio_R2[0].data
    
if 'None' not in files_N2:
    print('NII/HBeta ratio usable \n')
    line_ratios(files_N2[:-1],[files_N2[-1]],'H2_metal','N2_ratio')
    ratio_N2=fits.open('./H2_metal/N2_ratio.fits')
    metal_map=ratio_N2.copy()
    ratio_N2=ratio_N2[0].data
    
if 'None' not in files_S2:
    print('SII/HBeta ratio usable \n')
    line_ratios(files_S2[:-1],[files_S2[-1]],'H2_metal','S2_ratio')
    ratio_S2=fits.open('./H2_metal/S2_ratio.fits')
    ratio_S2=ratio_S2[0].data
    
if 'None' not in files_R3:
    print('OIII/HBeta ratio usable \n')
    line_ratios(files_R3[:-1],[files_R3[-1]],'H2_metal','R3_ratio')
    ratio_R3=fits.open('./H2_metal/R3_ratio.fits')
    ratio_R3=ratio_R3[0].data
    
#Dichotomy of the possibilities
    
#too many lines are missing
if 'None' in files_N2:
    print('N2 ratio unavailable. Cannot compute')
    case='RIP'
else:
    
    #Everything is fine
    if 'None' not in files_R2 and 'None' not in files_R3:
            print('Computing using the standard formulae')
            case='R,3D'
    
    #R2 line is missing
    if 'None' in files_R2 and 'None' not in files_R3 and 'None' not in files_S2:
            print('R2 ratio unavailable.\n'
              'Computing using the S2 ratio instead.\n')
            case='S,3D'
            
    #two lines are missing (2 cases)
    if 'None' in files_S2 and 'None' in files_R3:
        print('Some of the ratios are unavailable.\n'
              'Attempting a two-dimensional computation for the higher values of the NII ratios.\n'
              'Using the R relation\n')
        case='R,2D'
        
    if 'None' in files_R2 and 'None' in files_R3:
        print('Some of the ratios are unavailable.\n'
              'Attempting a two-dimensional computation for the higher values of the NII ratios.\n'
              'Using the S relation\n')
        case='S,2D'
    

#defining the calibration functions. The arguments are ordered according to the lines wavelengths

def OH_R(R2,N2,R3):
    
    OH_RU=8.589+0.022*np.log(R3/R2)+0.399*np.log(N2)+(-0.137+0.164*np.log(R3/R2)+0.589*np.log(N2))*np.log(R2)-12
    OH_RL=7.932+0.944*np.log(R3/R2)+0.695*np.log(N2)+(+0.970-0.291*np.log(R3/R2)-0.019*np.log(N2))*np.log(R2)-12
    
    return np.where(np.log(N2)>=-0.6,OH_RU,OH_RL)

def OH_S(N2,S2,R3):
    
    OH_SU=8.424+0.030*np.log(R3/S2)+0.751*np.log(N2)+(-0.349+0.182*np.log(R3/S2)+0.508*np.log(N2))*np.log(S2)-12
    OH_SL=8.072+0.789*np.log(R3/S2)+0.726*np.log(N2)+(+1.069-0.170*np.log(R3/S2)+0.022*np.log(N2))*np.log(S2)-12
    
    return np.where(np.log(N2)>=-0.6,OH_SU,OH_SL)

def OH_R_2D(R2,N2):
    
    OH_RU_2D=8.589+0.329*np.log(N2)+(-0.205+0.549*np.log(N2))*np.log(R2)-12
    
    return np.where(np.log(N2)>=-0.6,OH_RU_2D,np.nan)

def OH_S_2D(N2,S2):
    
    OH_SU_2D=8.445+0.699*np.log(N2)+(-0.253+0.217*np.log(N2))*np.log(S2)-12
    
    return np.where(np.log(N2)>=-0.6,OH_SU_2D,np.nan)

#computing the results depending of the case

if case !='RIP' :
    os.chdir('H2_metal')

if case=='R,3D':
    metal_map[0].data=OH_R(ratio_R2,ratio_N2,ratio_R3)
    
if case=='S,3D':
    metal_map[0].data=OH_S(ratio_N2,ratio_S2,ratio_R3)
    
if case=='R,2D':
    metal_map[0].data=OH_S(ratio_N2,ratio_S2,ratio_R3)

if case=='S,2D':
    metal_map[0].data=OH_S(ratio_N2,ratio_S2,ratio_R3)
    
metal_map.writeto('metal_map.fits',overwrite=True)

metal_img=Image('metal_map.fits')

metal_fig=plt.figure(figsize=(10,8))
   
plt.suptitle('metal map computed with method '+case)
metal_img.plot(zscale=True,colorbar='v')
