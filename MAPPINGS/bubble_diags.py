# @Author: Maxime Parra
# @Date:   06-05-2021
# @Email:  maxime.parra@irap.omp.eu

#NII/HBeta shock metallicity estimate
#comparison with the tables of Allen et al. 2008

import sys, os
sys.path.append('/home/mparra/PARRA/Scripts/Python/MUSE/')

from astropy.io import fits
import numpy as np
from line_utils import map_plot
import matplotlib.pyplot as plt
import logging

sdir='/home/mparra/PARRA/Observ/Andres/optical/bpt_diagrams_v2'

os.chdir(sdir)

'''bpt mask to keep only the bubble at the center'''

#no star forming regions in map 1
bpt_expmap1=fits.open('BPT_1.fits')
bpt_expmap1=np.where(bpt_expmap1[0].data>0,True,False)

#only AGN regions in map 2
bpt_expmap2=fits.open('BPT_2.fits')
bpt_expmap2=np.where(bpt_expmap2[0].data==2,True,False)

bpt_expmap=bpt_expmap1
for i in range(np.size(bpt_expmap,0)):
    for j in range(np.size(bpt_expmap[i])):
        bpt_expmap[i][j]=bpt_expmap1[i][j] and bpt_expmap2[i][j]
        
        #manual correction to take of the right pixels remaining and only keep
        #the bubble
        if j>62:
            bpt_expmap[i][j]=False
        
N2=fits.open('flux_NII6583_HALPHA.fits')
N2_out=N2.copy()
N2=N2[0].data

for i in range(np.size(N2,0)):
    for j in range(np.size(N2[i])):
        if not bpt_expmap[i][j]:
            N2[i][j]=np.nan
        
OIII=fits.open('flux_OIII5007_HBETA.fits')
OIII=OIII[0].data

os.chdir('..')
os.chdir('shock_metal/')

N2_out[0].data=np.log10(N2)

N2_out.writeto('NII_Halpha_flux_mask.fits',overwrite=True)

map_plot('NII_Halpha_flux_mask.fits',title='remaining bubble region and corresponding log(NII/Halpha) values')

#averaging for N2 (which isn't in log scale)
# the ~ is the complement function from C.
avg=np.average(N2[~np.isnan(N2)])

print('average N2 in the bubble : ',avg)



'''################################'''
'''comparison with the Allen tables'''
'''################################'''

os.chdir('/home/mparra/PARRA/prog/MAPPINGS3/tables/emission_line_ratios')

#namings for solar models in growing order of density
ordersolar=np.array(['T','U','M','V','L','S'])

names_LMC=[]
names_SMC=[]
names_solar_n1=[]
names_solar=[]
#obtaining all of the files for the SMC and LMC models (doc. in Allen et al. 2008)
for root, dirs, files in os.walk(".", topdown=False):        
    for name in files:
        if name.startswith('P'):
            names_SMC.append(name)
        if name.startswith('Q'):
            names_LMC.append(name)
        if name.startswith('M'):
            names_solar_n1.append(name)
    for letter in ordersolar:
        names_letter=[]
        for name in files:
            if name.startswith(letter):
                names_letter.append(name)
        names_solar.append(names_letter)

#order for elements. b0 is actually equal to b0.0001
#all those values are not used for every model but using a single array is easier
order=np.array([
    'b0_','b0_0001_','b0_001_','b0_01_','b0_05_','b0_1_','b0_2_','b0_32_','b0_4_','b0_5_','b0_632_','b1_','b1_0_',
    'b1_26_','b1_58_','b2_','b2_0_','be_','b3_16_','b4_','b4_0_','b5_','b5_0_','b6_32_','b10_','b10_2_','b12_65_',
    'b15_8_','b16_','b20_','b30_','b32_','b40_','b50_','b63_','b100_','b126_','b160_','b316_','b1000_'])

ordermod=np.array(['_s_','_p_','_sp_'])

#function to order the 1D-array in a more functional multi-dimensional one
def order_names(names):
    liste=np.array([[None]*int(len(names)/len(ordermod))]*len(ordermod))
    for elem in names:
        for mod in ordermod:
            if mod in elem:
                #position of the line where the elemen will be assigned
                l1=np.where(ordermod==mod)[0][0]
                
                #creation of a list of maximal size
                orderlist=np.array([None]*len(order))
                
                #filling the list with the existing values for the model
                for value in order:
                    for otherelem in names:
                        #here we add the strings to avoid doublons for values such as b_0 and b_0_5
                        if value+mod[1:] in otherelem:
                            orderlist[np.where(order==value)[0][0]]=otherelem
                
                #cleaning the remaining empty elements
                orderlist=orderlist[np.where(orderlist!=None)]
                
                #test to check if the remaining elements correspond to the initial number of elements 
                #(i.e. if all the elements have been found)
                if len(names)/len(ordermod)==len(orderlist):
                    for l2 in range(len(orderlist)):
                        liste[l1][l2]=orderlist[l2]
                else:
                    logging.error('Problem encountered while ordering.\n')
                    print('array:',names,'\n')
                    print('ordered and cut subarray:',orderlist,'\n')
                    sys.exit()
    return liste

   
list_LMC=order_names(names_LMC)
list_SMC=order_names(names_SMC)

list_solar_n1=order_names(names_solar_n1)

list_solar=np.array([None]*len(ordersolar))
for i in range(len(ordersolar)):
    list_solar[i]=order_names(names_solar[i])

'''
#obtaining the model parameters. For each array element (model), we search for the lines
in wavelengths and take the ratio values for the n lowest FWHM 
(the FWHM in our dataset is lower than 150)
'''

#wavelengths for OIII, Halpha and NII lines.
wavelengths=np.array(['5006.7700','6562.8000','6583.3400'])

#number of FWHM rows taken
n=10

def mod_loader(liste):
    mod_list=np.array([[[[None]*n]*len(wavelengths)]*np.size(liste,1)]*np.size(liste,0))
    for i in range(np.size(liste,0)):
        for j in range(np.size(liste[i])):
            currtxt=np.loadtxt(liste[i][j],dtype=str)
            for k in range(len(currtxt)):
                for line in wavelengths:
                    if line in currtxt[k]:
                        mod_list[i][j][np.where(wavelengths==line)[0][0]]=currtxt[k][4:4+10]
    return mod_list

'''
computing the ratios that will be used in the BPT
here the division is manually chosen from the wavelengths array sicne we want NII/Halpha as first
row (abscissa) and the OIII/Hbeta ratio as ordinate 
(which is already computed since the ratios given in the tables are in regard to Hbeta)
'''

def bpt_adapter(mod_list):
    tr_list=np.transpose(mod_list.astype(float),[2,0,1,3])
    ratio_list=np.array([tr_list[2]/tr_list[1],tr_list[0]])
    ratio_list=np.transpose(ratio_list,[1,0,2,3])
    return ratio_list
    
mod_LMC=mod_loader(list_LMC)
mod_SMC=mod_loader(list_SMC)

ratio_LMC=bpt_adapter(mod_LMC)
ratio_SMC=bpt_adapter(mod_SMC)

lratio_LMC=np.log10(ratio_LMC)
lratio_SMC=np.log10(ratio_SMC)

mod_solar_n1=mod_loader(list_solar_n1)
ratio_solar_n1=bpt_adapter(mod_solar_n1)
lratio_solar_n1=np.log10(ratio_solar_n1)

mod_solar=np.array([None]*len(ordersolar))
ratio_solar=np.array([None]*len(ordersolar))
lratio_solar=np.array([None]*len(ordersolar))

for i in range(len(ordersolar)):
    mod_solar[i]=mod_loader(list_solar[i])
    ratio_solar[i]=bpt_adapter(mod_solar[i])
    lratio_solar[i]=np.log10(ratio_solar[i])

    
#%%graphs


title=['Shocks-only Model', 'Precursor-only Model','Shocks+Precursor Model']

fwhm_list=[100,125,150,175,200,225,250,275,300,325,350,375,400,425,450,475,500,525,550,600,625,650,675,725,750,775,
       850,925,950]

titlend='FWHMs up to '+str(fwhm_list[n-1])+' km/s'

'''n= 1 /cm³ model comparison'''

for i in range(3):
    fig_allen,ax_allen=plt.subplots(1,figsize=(10,8))
    fig_allen.suptitle(title[i]+' for different metallicites and a fixed density, '+titlend)
    ax_allen.set_xlim(-1.5,0.8)
    ax_allen.set_ylim(-1,1.3)
    
    ax_allen.set_xlabel("log([NII]/H$_\\alpha$)")
    ax_allen.set_ylabel("log([OIII]/H$_\\beta$)")
    
    ax_allen.scatter(np.log10(N2),np.log10(OIII),color='lightgrey',marker='+',label='ULX bubble')
    
    ax_allen.plot(lratio_SMC[i][0],lratio_SMC[i][1],color='lightgreen')
    ax_allen.plot(np.transpose(lratio_SMC[i][0]),np.transpose(lratio_SMC[i][1]),color='lightgreen',linestyle=':')
    ax_allen.plot(10,10,color='lightgreen',label='SMC metallicity, n=1 cm$^{-3}$')
    
    ax_allen.plot(lratio_LMC[i][0],lratio_LMC[i][1],color='darkgreen')
    ax_allen.plot(np.transpose(lratio_LMC[i][0]),np.transpose(lratio_LMC[i][1]),color='darkgreen', linestyle=':')
    ax_allen.plot(10,10,color='darkgreen',label='LMC metallicity, n=1 cm$^{-3}$')

    ax_allen.plot(lratio_solar_n1[i][0],lratio_solar_n1[i][1],color='red')
    ax_allen.plot(np.transpose(lratio_solar_n1[i][0]),np.transpose(lratio_solar_n1[i][1]),color='red',linestyle=':')
    ax_allen.plot(10,10,color='red',label=r'solar metallicity, n=1 cm$^{-3}$')
    
    
    ax_allen.legend()
    plt.tight_layout()

'''solar models comparison'''

solar_colors=['darkblue','lightblue','red','violet','orange','gold']

densities=['0.01','0.1','1','10','100','1000']

for i in range(3):
    fig_allen,ax_allen=plt.subplots(1,2,figsize=(16,8))
    fig_allen.suptitle(title[i]+' for different densities and solar metallicity, '+titlend)
    
    for j in range(2):
        ax_allen[j].set_xlim(-1.5,0.8)
        ax_allen[j].set_ylim(-1,1.3)
        
        ax_allen[j].set_xlabel("log([NII]/H$_\\alpha$)")
        ax_allen[j].set_ylabel("log([OIII]/H$_\\beta$)")
        
        ax_allen[j].scatter(np.log10(N2),np.log10(OIII),color='lightgrey',marker='+',label='ULX bubble')
        
        for k in range(3):
            l=j*3+k
            ax_allen[j].plot(lratio_solar[l][i][0],lratio_solar[l][i][1],color=solar_colors[l])
            ax_allen[j].plot(np.transpose(lratio_solar[l][i][0]),np.transpose(lratio_solar[l][i][1]),
                          color=solar_colors[l],linestyle=':')
            ax_allen[j].plot(10,10,color=solar_colors[l],
                          label=r'solar metallicity, n='+densities[l]+' cm$^{-3}$')
            
        ax_allen[j].legend()
        
    plt.tight_layout()

