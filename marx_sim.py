import os
from xspec import *
import numpy as np
from regions import *
from astropy.io import fits
from astropy.time import Time
#Don't forget to initialize heasoft and marx in your console before launching the script/spyder
#Don't forget to launch this on your ciao python interpreter
#marx needs to be installed

'''
Currently the detector is automatically set to ACIS-S, since it's the most
likely to be used and there is no information in the header to differentiate 
between ACIS-S and ACIS-I
'''

AllModels.clear()
AllData.clear()
#starting dir
sdir='/home/mparra/PARRA/Observ/Andres/2950/repro/'
os.chdir(sdir)
#starting file
# file='/home/mparra/PARRA/Observ/Andres/2950/repro/ACIS_min20.fits'

file='ngc1313x1_extraction.pi'
# file='/home/mparra/PARRA/Observ/Andres/2950/repro/ngc1313x1_fk5.pi'

#reading the region file, which is assumed to have the same name as the file
regfile=file[:file.find('.')]+'.reg'

source_reg=read_ds9(regfile)[0]

#coordinates of the region file centroid, will be used as position for the marx simulation
source_ra=str(source_reg.center.ra.deg)
source_dec=str(source_reg.center.dec.deg)


os.chdir(sdir)

fits_file=fits.open(file)

#coordinates of nominal aimpoint
ra_nom=str(fits_file[0].header["RA_NOM"])
dec_nom=str(fits_file[0].header["DEC_NOM"])
roll_nom=str(fits_file[0].header["ROLL_NOM"])

#date of observation
date=str(round(Time(fits_file[0].header["DATE-OBS"]).decimalyear,3))

#grating of observation
grating=fits_file[0].header["GRATING"]

#reading spectral file
sp=Spectrum(file)

expos=str(round(sp.exposure,1))

Plot.commands=()
Plot.xAxis='kev'
Plot.device='/xs'
Plot('data')
e_list=Plot.x()
c_list=Plot.y()

Plot.device='none'

#creating the spectrum file which will serve as an argument for marx
f_as=open('marx_arg.txt','w+')
for i in range(len(e_list)):
    f_as.write('%f     %f \n' %(e_list[i],c_list[i]))
f_as.close()

#SourceFlux=-1 tells marx to take the integrated flux of the Spectrum

#Verbose=yes to have diagnostic messages

os.system('punlearn marx')

params='pset marx SpectrumType="FILE" SourceFlux=-1 SpectrumFile=marx_arg.txt'\
          +' ExposureTime='+expos+' Tstart='+date+' Outputdir=marxsim GratingType='+grating\
          +' DetectorType="ACIS-S" RA_Nom='+ra_nom+' Dec_Nom='+dec_nom\
          +' Roll_Nom='+roll_nom+' SourceRA='+source_ra+' SourceDEC='+source_dec+' Verbose=yes mode=h'
os.system(params)

nsim=100
os.system('mkdir -p marxouts')

for i in range(100):
    os.system('marx')
        
    os.chdir(sdir+'marxsim')
    fitsname='marxsim_'+str(i)+'.fits'
    os.system('marx2fits ./ '+fitsname)
    os.system('mv '+fitsname+' ../marxouts')
    os.system('rm *')
    os.chdir('..')
