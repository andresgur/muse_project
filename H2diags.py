# @Author: Maxime Parra
# @Date:   02-05-2021
# @Email:  maxime.parra@irap.omp.eu

'''
Strong line metallicity and ionization parameter estimations from Pilyugin et al. 2016 and Dors et al. 2017.
Since those methods are only usable on photo-ionized regions, we use the BPT diagnostic of OIII and NII as a limit : 
we only consider spaxels in the star formation and intermediate zones 
(following the logic of Dors et al. 2017, but with the newer BPT version of Law et al. 2021)
papers : 
Pilyugin et al. 2016 : https://academic.oup.com/mnras/article/457/4/3678/2589035
Dors et al. 2017 : doi.org/10.1093/mnras/stw3115
Law et al. 2021 : https://arxiv.org/abs/2011.06012

Note : this code will automatically search for automatic bpt outputs from bpt_newer.py
'''

#should be updated to use ratio_maker

import sys, os
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/')

import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import logging
import argparse
from mpdaf.obj import Image
from line_utils import line_ratios, ratio_maker

ap = argparse.ArgumentParser(description='Create BPT diagram from two given line ratio maps and the BPT diagram type. Separatios based on Law et al. 2021')
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='H2_diags', type=str)
ap.add_argument("-bpt", "--bptmap", nargs='?', help="BPT OIII/SII map file path", default=None, type=str)
args = ap.parse_args()

outdir=args.outdir

sdir='/home/mparra/PARRA/Observ/Andres/optical/'
os.chdir(sdir)

'''
First step : obtaining the BPT map of interest if it isn't provided
'''
if args.bptmap==None:
    for root, dirs, files in os.walk(".", topdown=False):        
                            for name in files:
                                if name.endswith('BPT_2.fits'):
                                        bpt_file=os.path.join(root,name)

else:
    bpt_file=args.bptmap
    
bpt_expmap=fits.open(bpt_file)
bpt_expmap=np.where(bpt_expmap[0].data<2,True,False)

'''
This part computes metallicity maps.
There are three methods depending on the amount of lines available.
Those computations are based on 4 line INTENSITY ratios.
We search for those lines using ratio_maker.
'''

# in order : R_2, N_2, S_2, R_3
pilyugin_lines=\
          [[['OII3727','OII3729'],['HBETA']],
           [['NII6548','NII6583'],['HBETA']],
           [['SII6716','SII6731'],['HBETA']],
           [['OIII4959','OIII5007'],['HBETA']]]

pilyugin_results,pilyugin_names=ratio_maker(pilyugin_lines,'int',outdir)
    
#We do not loop this to keep clear variable names
if pilyugin_results[0]=='Done':    
    ratio_R2=fits.open(pilyugin_names[0])
    ratio_R2=ratio_R2[0].data
    
if pilyugin_results[1]=='Done':    
    ratio_N2=fits.open(pilyugin_names[1])
    metal_file=ratio_N2.copy()
    ratio_N2=ratio_N2[0].data
    
if pilyugin_results[2]=='Done':    
    ratio_S2=fits.open(pilyugin_names[2])
    ratio_S2=ratio_S2[0].data
    
if pilyugin_results[3]=='Done':    
    ratio_R3=fits.open(pilyugin_names[3])
    ratio_R3=ratio_R3[0].data
    
#Dichotomy of the possibilities :
#too many lines are missing
if pilyugin_results[1]=='Unavailable':
    print('N2 ratio unavailable. Cannot compute')
    case='RIP'
    sys.exit()
    
else:

    #We got everything
    if pilyugin_results[0]==pilyugin_results[3]=='Done':
            print('Computing using the standard formulae')
            case='R,3D'
    
    #R2 line is missing
    if pilyugin_results[0]=='Unavailable' and pilyugin_results[2]==pilyugin_results[3]=='Done':
            print('R2 ratio unavailable.\n'
              'Computing using the S2 ratio instead.\n')
            case='S,3D'
            
    #two lines are missing (2 cases)
    if pilyugin_results[2]==pilyugin_results[3]=='Unavailable':
        print('Some of the ratios are unavailable.\n'
              'Attempting a two-dimensional computation for the higher values of the NII ratios.\n'
              'Using the R relation\n')
        case='R,2D'
        
    if pilyugin_results[0]==pilyugin_results[3]=='Unavailable':
        print('Some of the ratios are unavailable.\n'
              'Attempting a two-dimensional computation for the higher values of the NII ratios.\n'
              'Using the S relation\n')
        case='S,2D'
    

#defining the calibration functions. The arguments are ordered according to the lines wavelengths

def OH_R(R2,N2,R3):
    
    OH_RU=8.589+0.022*np.log(R3/R2)+0.399*np.log(N2)+(-0.137+0.164*np.log(R3/R2)+0.589*np.log(N2))*np.log(R2)
    OH_RL=7.932+0.944*np.log(R3/R2)+0.695*np.log(N2)+(+0.970-0.291*np.log(R3/R2)-0.019*np.log(N2))*np.log(R2)
    
    return np.where(np.log(N2)>=-0.6,OH_RU,OH_RL)

def OH_S(N2,S2,R3):
    
    OH_SU=8.424+0.030*np.log(R3/S2)+0.751*np.log(N2)+(-0.349+0.182*np.log(R3/S2)+0.508*np.log(N2))*np.log(S2)
    OH_SL=8.072+0.789*np.log(R3/S2)+0.726*np.log(N2)+(+1.069-0.170*np.log(R3/S2)+0.022*np.log(N2))*np.log(S2)
    
    return np.where(np.log(N2)>=-0.6,OH_SU,OH_SL)

def OH_R_2D(R2,N2):
    
    OH_RU_2D=8.589+0.329*np.log(N2)+(-0.205+0.549*np.log(N2))*np.log(R2)
    
    return np.where(np.log(N2)>=-0.6,OH_RU_2D,np.nan)

def OH_S_2D(N2,S2):
    
    OH_SU_2D=8.445+0.699*np.log(N2)+(-0.253+0.217*np.log(N2))*np.log(S2)
    
    return np.where(np.log(N2)>=-0.6,OH_SU_2D,np.nan)

#computing the results depending of the case

os.chdir(outdir)

metal_map=metal_file[0].data

if case=='R,3D':
    metal_map=OH_R(ratio_R2,ratio_N2,ratio_R3)
    
if case=='S,3D':
    metal_map=OH_S(ratio_N2,ratio_S2,ratio_R3)
    
if case=='R,2D':
    metal_map=OH_S(ratio_N2,ratio_S2,ratio_R3)

if case=='S,2D':
    metal_map=OH_S(ratio_N2,ratio_S2,ratio_R3)

if bpt_expmap.shape!=metal_map.shape:
    logging.warning("Can't compare the bpt map and the metal map : Their shapes are different.\n"
          "Using the entire metal map as a result. Be careful.")
    
for i in range(np.size(metal_map,0)):
    for j in range(np.size(metal_map[i])):
        if not bpt_expmap[i][j]:
            metal_map[i][j]=np.nan
            
metal_file[0].data=metal_map

metal_file.writeto('metal_map.fits',overwrite=True)

metal_img=Image('metal_map.fits')

metal_fig=plt.figure(figsize=(10,8))
   
plt.suptitle('metal map computed with method '+case)
metal_img.plot(zscale=True,colorbar='v')

'''
This part computes ionization parameter maps.
There are three methods depending on the amount of lines available.
Those computations are based on 4 line INTENSITY ratios.
We search for those lines using ratio_maker.
'''

#defining the calibration function. Here, Z is the metallicity normalized to solar metallicity

#there are significant uncertainties here, be careful

Z_star=8.69

def U_ion(Z,S2):
    a_ion=-0.26*Z-1.54
    b_ion=-3.69*Z**2+5.11*Z-5.26
    return a_ion*S2+b_ion

ion_par_file=metal_file.copy()

ion_par_map=ion_par_file[0].data

ion_par_map=U_ion(metal_map/Z_star,np.log(ratio_S2))

ion_par_file[0].data=ion_par_map

ion_par_file.writeto('ion_par.fits',overwrite=True)

ion_par_img=Image('ion_par.fits')

ion_par_fig=plt.figure(figsize=(10,8))

plt.suptitle('ionization map computed with method' +case)
ion_par_img.plot(zscale=True,colorbar='v')

