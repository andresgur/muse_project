#%% initialisation
import os
from xspec import *
import numpy as np
from regions import *
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt

#marx needs to be installed (https://space.mit.edu/cxc/marx/)
#saotrace needs to be installed (https://cxc.harvard.edu/cal/Hrma/Raytrace/SAOTrace.html)

#Don't forget to initialize heasoft, marx and saotrace in your console before launching the script/spyder
#Don't forget to launch this on your ciao python interpreter

#saotrace might not work if you have conflicts in your perl libraries : initialising saotrace might
#replace one (but not all) of your perl libraries environment variables. 
#You can manually setup the other variables in your saotrace starting script

'''
This version of the script can use the advanced modelisation from SAOTRACE and thus reprocesses
part of the data to create its parameter files. This means the starting directory should be
set to the reprocessed directory

Some aspects of the scripts are specifically tailored towards POINT SOURCES
(arf correction and marx paremeters especially). Modifications are needed to use this on
extended sources.

Please input real region files and not links

No background included as of now in the simulation parameters (it is included in the comparison after)
(why are you even checking your psf if your source isn't bright enough to stand out)

More  details about the comparison at the end of the script
'''

''' PARAMETERS '''

#Saotrace ('yes'), Marx ('no'), or choice ('maybe')
sao_arg='no'

#number of simulations
nsim=10000

#starting dir
# sdir='/home/mparra/PARRA/Observ/Andres/2950/repro/'
sdir='/home/mparra/PARRA/Observ/Andres/4748/repro/'

#starting region
# regfile='ngc1313x1_extraction.reg'
regfile='ngc1313x1_fk5.reg'

#background region (will be used later)
bgregfile='background.reg'

#number of annuli
nbins=10

#maximum annuli size ratio compared to the maximal dimension of the source region
max_factor=2

#chandra pixel to arsec ratio
pta=0.492

#radius affected by pile-up in the initial data in arsecs
# rlim_pileup=1.2
rlim_pileup=1.0
 
'''############'''

if sao_arg=='maybe':
    lock=0
    while lock==0:
        decis=input('Do you want to use Saotrace ("yes"), or Marx ("no") ?')
        if decis in ['yes','no']:
            lock=1
else:
    decis=sao_arg
    
if decis=='yes':
    decis_int=1
elif decis=='no':
    decis_int=0
    

os.chdir(sdir)

#ditherfile, there should only be one ditherfile for an observation
for file in os.listdir('./'):
        if file.endswith('_asol1.fits'):
            ditherfile=file
        if file.endswith('_evt2.fits'):
            evt2=file
            
#reading the data file, which is assumed to have the same name as the region file
name=regfile[:regfile.find('.')]
data_file=name+'.pi'

#reading the region file
source_reg=read_ds9(regfile)[0]

#sky coordinates of the region file centroid, will be used as position for the marx simulation
source_ra=str(source_reg.center.ra.deg)
source_dec=str(source_reg.center.dec.deg)

'''
the region doesn't actually enclose the full counts of the source, so we can manually
correct by comparison to a full region, which here has been manually evaluated
to a 7" circle region centered on the other region's centroid for the 2950 observation
'''

if regfile=='ngc1313x1_extraction.reg':
    cor_factor=0.94
else:
    cor_factor=1

#reading the event file to get some more parameters
fits_evt2=fits.open(evt2)

#instrument used
instrum=fits_evt2[0].header['INSTRUME']

#coordinates of nominal aimpoint
ra_nom=str(fits_evt2[1].header["RA_NOM"])
dec_nom=str(fits_evt2[1].header["DEC_NOM"])
roll_nom=str(fits_evt2[1].header["ROLL_NOM"])

#date of observation
date=str(round(Time(fits_evt2[1].header["DATE-OBS"]).decimalyear,3))

#grating of observation
grating=fits_evt2[1].header["GRATING"]

#file reprocessing to recreate evenly binned event files and arf to get the best results possible
os.system('dmstat "'+evt2+'[sky=region('+regfile+')][bin sky=1]" centroid=yes')

#computing the source center
stdout=Popen('pget dmstat out_cntrd_phys', shell=True, stdout=PIPE).stdout
ctr_output=str(stdout.read())

#those are the estimated coordinates of the source center, based on the brightest pixel in the source region
ctr_x=ctr_output[2:ctr_output.find(',')]
ctr_y=ctr_output[ctr_output.find(',')+1:-3]

#computing some more informations about the source center
os.system('punlearn dmcoords')
os.system('dmcoords '+evt2+' op=sky x='+ctr_x+' y='+ctr_y)
stdout=Popen('pget dmcoords chip_id', shell=True, stdout=PIPE).stdout
chip_output=str(stdout.read())
chip_id=int(chip_output[2:-3])

#such as the center position in sky coordinates
stdout=Popen('pget dmcoords ra dec', shell=True, stdout=PIPE).stdout
coord_output=str(stdout.read())[2:-3]

ctr_ra=coord_output[:coord_output.find('\\n')]
ctr_dec=coord_output[coord_output.find('\\n')+2:]
ctr_coord=SkyCoord(ctr_ra,ctr_dec,unit=(u.hourangle,u.deg),frame='fk5')
ctr_ra=str(ctr_coord.ra.deg)
ctr_dec=str(ctr_coord.dec.deg)

#We now also know which part of the instrument is used for the detection
if instrum=='ACIS':
    if chip_id in [0,1,2,3]:
        detector='ACIS-I'
        subsys='ACIS-I'+str(chip_id)
    elif chip_id in [4,5,6,7,8,9]:
        detector='ACIS-S'
        subsys='ACIS-S'+str(chip_id-4)
    
if instrum=='HRC':
    detector=input('Since the CCD number overlap in this instrument, please confirm'\
                   +' if the detecor used is HRC-I or HRC-S')
        
os.system('asphist infile="'+ditherfile+'" outfile="marx.asp" evtfile="'+evt2+'" clobber=yes')

#here we only use mkarf to create the response file, a more accurate version can be computed with mkwarf
#but is much more complex to create

binning='0.1'

os.system('mkarf detsubsys="'+subsys+'" grating="'+grating+'" outfile="marx.arf" obsfile="'+evt2\
          +'" asphistfile="marx.asp" sourcepixelx='+ctr_x+' sourcepixely='+ctr_y+' engrid="0.3:10.0:'+binning
          +'" maskfile=NONE pbkfile=NONE dafile=NONE verbose=1 mode=h clobber=yes')

#not done in the marx tutorial but should be done for a point source...right ? 
os.system('arfcorr infile="'+evt2+'[sky=region('+regfile+')][bin sky]" region="region('+regfile+\
          ')" x='+ctr_x+' y='+ctr_y+' arf=marx.arf outfile=marx_corr.arf clobber=yes')
    
#if you want to create an image file to check if the detection is correctly centered    
# os.system('dmcopy "'+evt2+'"[EVENTS][bin x='+str(int(float(ctr_x)-100))+':'+str(int(float(ctr_x)+100))\
#           +':1,y='+str(int(float(ctr_y)-100))+':'+str(int(float(ctr_y)+100))\
#           +':1" image_src.fits option=image clobber=yes')

#source energy computation
os.system('dmtcalc "'+evt2+'[EVENTS][sky=region('+regfile+')]" marx_evt2_energy.fits'\
          +' expr="energy=(float)pi*0.0149" clobber=yes')
    
#bg energy computation. Should be negligible, but included nonetheless
os.system('dmtcalc "'+evt2+'[EVENTS][sky=region('+bgregfile+')]" marx_bg_energy.fits'\
          +' expr="energy=(float)pi*0.0149" clobber=yes')
    
#net histogram computation with real energy binning
os.system('dmextract "marx_evt2_energy.fits[bin energy=0.3:9.999:'+binning+']"'\
          +' bkg="marx_bg_energy.fits" outfile=marx_extract.spec clobber=yes opt=generic')
    
#reading the newly created files
arf_corr=fits.open('marx_corr.arf')
arf_data=arf_corr[1].data
extract=fits.open('marx_extract.spec')
extract_data=extract[1].data

#obtaining the exposure
expos=str(round(extract_data['BG_EXPOSURE'][0],1))

'''
Now we create the spectrum file which will serve as an argument for marx/sao
For Marx the (received) flux column is supposed to be a normalized flux DENSITY so the net rate column needs 
to be divided by the binning, on top of the effective area using the arf and 
the correction factor to account for the loss percentage.
'''

if decis=='yes':
    
    f_sao=open('marx_sao.rdb','w+')
    
    flux=extract_data['NET_RATE']/(arf_data['SPECRESP']*cor_factor)
    
    f_sao.write('#argument file for SAOTrace, manually created from marx_extract.spec and marx_corr.arf\n')
    f_sao.write('ENERG_LO\tENERG_HI\tFLUX\n')
    for i in range(len(flux)):
        
        f_sao.write('%f\t%f\t%f\n' %(arf_data['ENERG_LO'][i],arf_data['ENERG_HI'][i],flux[i]))
    
    f_sao.close()
    
    #sao parameters file
    
    fpar_sao=open('marx_sao_par.lua','w+')
    fpar_sao.write('ra_pnt='+ra_nom+'\ndec_pnt='+dec_nom+'\nroll_pnt='+roll_nom+'\n\n'\
                  +'dither_asol_chandra{file="'+ditherfile+'", ra=ra_pnt, dec=dec_pnt, roll=roll_pnt}\n\n'
                  +'point{position={ra='+ctr_ra+', dec='+ctr_dec+',ra_aimpt=ra_pnt, dec_aimpt=dec_pnt,},\n'
                  +'spectrum={{file="marx_sao.rdb", units="photons/s/cm2", scale=1, format="rdb", emin="ENERG_LO",'
                  +'emax="ENERG_HI",flux="FLUX"}}}')
    fpar_sao.close()
    
    #sao not working as of now
    
    os.system('trace-nest tag=saotrace srcpars=marx_sao_par.lua tstart='+date+' limit='+expos+' limit_type=sec')
    
    os.system('punlearn marx')
    
    #marx parameters : 
    #SourceFlux=-1 tells marx to take the integrated flux of the Spectrum
    #Verbose=yes to have diagnostic messages
    
    par_marx='pset marx SourceType="SAOSAC" SourceFlux=-1 SAOSACFile=saotrace.fits'\
              +' ExposureTime='+expos+' TStart='+date+' OutputDir=sao_sim GratingType='+grating\
              +' DetectorType="ACIS-S" RA_Nom='+ra_nom+' Dec_Nom='+dec_nom\
              +' DitherModel="FILE" DitherFile='+ditherfile\
              +' Roll_Nom='+roll_nom+' SourceRA='+source_ra+' SourceDEC='+source_dec+' Verbose=yes mode=h'    

elif decis=='no':
    
    f_marx=open('marx_arg.tbl','w+')
    
    fluxdens=extract_data['NET_RATE']/(arf_data['SPECRESP']*cor_factor*float(binning))
    
    f_marx.write('#argument file for larx, manually created from marx_extract.spec and marx_corr.arf\n')
    for i in range(len(fluxdens)):
        
        f_marx.write('%f %f\n' %(arf_data['ENERG_HI'][i],fluxdens[i]))
    
    f_marx.close()
    
    os.system('punlearn marx')
    
    #marx parameters : 
    #SourceFlux=-1 tells marx to take the integrated flux of the Spectrum
    #Verbose=yes to have diagnostic messages

    par_marx='pset marx SpectrumType="FILE" SourceFlux=-1 SpectrumFile=marx_arg.tbl'\
              +' ExposureTime='+expos+' TStart='+date+' OutputDir=marx_sim GratingType='+grating\
              +' DetectorType="ACIS-S" RA_Nom='+ra_nom+' Dec_Nom='+dec_nom\
              +' DitherModel="FILE" DitherFile='+ditherfile\
              +' Roll_Nom='+roll_nom+' SourceRA='+source_ra+' SourceDEC='+source_dec+' Verbose=yes mode=h'

os.system(par_marx)

#%% simulations

os.chdir(sdir)

#output directory
outdirs=['marx_outs','sao_outs']

#temporary outputs directory, will be deleted
simdirs=['marx_sim','sao_sim']

fitsname=['marxsim_','saosim_']

os.system('mkdir -p '+outdirs[decis_int])
os.system('mkdir -p '+simdirs[decis_int])

for i in range(1,nsim+1):
    os.system('marx')
    
    os.chdir(simdirs[decis_int])
    curr_fitsname=fitsname[decis_int]+str(i)+'.fits'
    
    #running a pileup command to emulate the loss of counts due to pileup
    os.system('marxpileup MarxOutputDir="./"')
    os.chdir('pileup')
    
    #creating a fits file from the outputs. This is the only output we're interested in
    #so we move it and start again (marx and marxpileup both overwrite files by default)
    os.system('marx2fits ./ '+curr_fitsname)
    os.system('mv '+curr_fitsname+' ../../'+outdirs[decis_int])
    os.chdir('../..')

#deleting all the temporary files    
os.chdir(simdirs[decis_int]+'/pileup')
os.system('rm *')
os.chdir('..')
os.system('rmdir pileup')
os.system('rm *')
os.chdir('..')
os.system('rmdir '+simdirs[decis_int])

#moving the initialisation files 
os.system('mv *marx* '+outdirs[decis_int])

#%%comparison with the data

'''
We  divide our regions in nbins (default 10) annuli ranging from the source center to 
twice the maximal dimension of the source region (which should be an ellipse or a circle in sky region)

At the end we renormalize the outer parts of the PSF to reduce errors. 
The inner limit results from manual computations of the proportion of the initial data affected by pile-up.

We also use our own computation of the source center instead of the source centroid
 
-> You'll need to change a few things to use ellipsoid annuli

Moreover, the code expects a background region containing a single region (in sky coordinates). 
If you use composite background regions you'll need to adapt the very last part
'''

if type(source_reg).__name__=='EllipseSkyRegion':
    largestsize=np.maximum(source_reg.width,source_reg.height)
elif type(source_reg).__name__=='CircleSkyRegion':
    largestsize=source_reg.radius
    
maxsize=largestsize*max_factor

#we will write the command line in pixels to avoid complications
maxsize_pixel=round(maxsize.value/pta,2)
step=maxsize_pixel/nbins

#First we extract the annuli for the initial file, this time with a background file
#however, chandra can't read directly the background file in our case so we rewrite it in an
#acceptable format

f_bg=open(bgregfile,'r')
bg_lines=f_bg.readlines()

#sometimes ds9 doesn't bother including the format (fk5 for example) when saving regions, so we test
#the length of the third line. If it's long, it's probably the region line instead of a format
if len(bg_lines[2])<20:
    bgreg_line=bg_lines[3]
else :
    bgreg_line=bg_lines[2]

f_bg_marx=open('marx_bg.reg','w+')
f_bg_marx.write(bgreg_line[:bgreg_line.find('#')])
f_bg_marx.close()
    
#string used in the command line
str_annulus='[bin sky=annulus('+str(ctr_x)+','+str(ctr_y)+',0:'\
              +str(maxsize_pixel)+':'+str(step)+')]'
              
#With this we can extract the entire concentric annuli in one go
os.system('dmextract "'+evt2+str_annulus+'" bkg="'+evt2+'[bin sky=@marx_bg.reg]"'\
          +' outfile=curr_profile.fits opt=generic clobber=yes')

#after which we store the surface brillance values and associated 1 sigma errors in arrays
with fits.open('curr_profile.fits') as curr_profile :
    sbdata=curr_profile[1].data['SUR_BRI']
    sbdata_err=curr_profile[1].data['SUR_BRI_ERR']

#3 sigma intervals on the initial data (not necessary for the plot but can be called after the script 
#for the values)
sig1_data=np.transpose(np.array([sbdata-1*sbdata_err,sbdata+1*sbdata_err]))
sig3_data=np.transpose(np.array([sbdata-3*sbdata_err,sbdata+3*sbdata_err]))

'''computations for the renormalization'''

#conversion of the renormalisation annuli's inner and outer radius
rlim_pixel=round(rlim_pileup/pta,2)
largestsize_pixel=round(largestsize.value,2)/pta

#second version of the command line string
str_annulus_renorm='[bin sky=annulus('+str(ctr_x)+','+str(ctr_y)+','+str(rlim_pixel)+':'\
              +str(largestsize_pixel)+':'+str(largestsize_pixel-rlim_pixel)+')]'

#extracting the annulus for the data
os.system('dmextract "'+evt2+str_annulus_renorm+'" bkg="'+evt2+'[bin sky=@marx_bg.reg]"'\
          +' outfile=curr_profile.fits opt=generic clobber=yes')

#now we can get the renormalization data value
with fits.open('curr_profile.fits') as curr_profile :
    sbdata_norm=curr_profile[1].data['SUR_BRI'][0]

#we finish by moving the products to the out directory
os.system('mv curr_profile.fits marx_bg.reg '+outdirs[decis_int])

'''evaluations for the simulations'''

#we start by creating a single region with every single annuli region inside for visualisation
#since write_ds9's format isn't usable by Chandra, we write the regions manually
os.chdir(sdir+outdirs[decis_int])

radii=np.linspace(0,maxsize,num=nbins+1)

ra_tuple=(int(ctr_coord.ra.hms[0]),int(ctr_coord.ra.hms[1]),str(round(ctr_coord.ra.hms[2],4)))
dec_tuple=(int(ctr_coord.dec.dms[0]),abs(int(ctr_coord.dec.dms[1])),str(abs(round(ctr_coord.dec.dms[2],4))))

dec_init_tuple=(int(ctr_coord.dec.dms[0]),int(abs(ctr_coord.dec.dms[1])),\
                str(abs(round(ctr_coord.dec.dms[2],4))),str(round(radii[1].value,4)))

f_reg=open('annuli.reg','w+')

f_reg.write(
'# Region file format: DS9 version 4.1\n'\
+'#global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1'\
+' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')

for  i in range(nbins):

    if i==0:
        curr_tuple=ra_tuple+dec_init_tuple
        f_reg.write('circle(%i:%i:%s,%i:%i:%s,%s")\n' %curr_tuple)
    else:
        curr_tuple=ra_tuple+dec_tuple+(str(round(radii[i].value,4)),str(round(radii[i+1].value,4)))
        f_reg.write('annulus(%i:%i:%s,%i:%i:%s,%s",%s")\n' %curr_tuple)

f_reg.close()

os.system('punlearn dmextract')


#Now, we extract the renormalization annuli and the concentric annuli in a similar manner than with the data, 
#adding the renormalization

sblist=np.array([[None]*nbins]*nsim)
sblist_err=np.array([[None]*nbins]*nsim)
sblist_norm=np.array([None]*nsim)
             

for i in range(1,nsim+1):
    currfile=fitsname[decis_int]+str(i)+'.fits'
    
    os.system('dmextract "'+currfile+str_annulus_renorm+'" outfile=curr_profile.fits opt=generic clobber=yes')
    with fits.open('curr_profile.fits') as curr_profile:
        sblist_norm[i-1]=curr_profile[1].data['SUR_BRI'][0]
        
    os.system('dmextract "'+currfile+str_annulus+'" outfile=curr_profile.fits opt=generic clobber=yes')
    with fits.open('curr_profile.fits') as curr_profile:
        sblist[i-1]=curr_profile[1].data['SUR_BRI']*sbdata_norm/sblist_norm[i-1]
        sblist_err[i-1]=curr_profile[1].data['SUR_BRI_ERR']*sbdata_norm/sblist_norm[i-1]
    

#ordering and transposing the surface brillance list for easier manipulation    
ord_sblist=np.sort(np.transpose(sblist))

argord_sblist=np.argsort(np.transpose(sblist))

#we will order it according to the other list's order in the next loop
ord_sblist_err=np.transpose(sblist_err)

#1 and 3 sigma values
sig1_val=0.6823
sig3_val=0.9973

#we could use ceil here but since the indexation begins at 0 we need to go back 1 index
#anyway so we might aswell directly use floor
ind1_up=int(np.floor(nsim*(1+sig1_val)/2))
ind1_lo=int(np.floor(nsim*(1-sig1_val)/2))
ind3_up=int(np.floor(nsim*(1+sig3_val)/2))
ind3_lo=int(np.floor(nsim*(1-sig3_val)/2))

#statistical analysis on the arrays
sig1_list=np.array([[None,None]]*nbins)
sig3_list=np.array([[None,None]]*nbins)

'''Wierdly enough, it seems that the images from Marx are often identical'''

for i in range(nbins):
    
    ord_sblist_err[i]=ord_sblist_err[i][argord_sblist[i]]
    
    sig1_list[i][0]=ord_sblist[i][ind1_lo]-1*ord_sblist_err[i][ind1_lo]
    sig1_list[i][1]=ord_sblist[i][ind1_up]+1*ord_sblist_err[i][ind1_up]
    
    sig3_list[i][0]=ord_sblist[i][ind3_lo]-3*ord_sblist_err[i][ind3_lo]
    sig3_list[i][1]=ord_sblist[i][ind3_up]+3*ord_sblist_err[i][ind3_up]

os.chdir('..')   

#%% graphs

fig_marx,ax_marx=plt.subplots(constrained_layout=True,figsize=(10,8))

fig_marx.suptitle(name+' radial extension simulation with '+str(nsim)+' point source Marx simulations')

ax_marx.set_xlabel('source size fraction')
ax_marx.set_ylabel(r'net surface brightness ([cts/s/cm$^2$])')

ax_marx.set_yscale('log')


#transposition and taking off negative values to have an easier time plotting
tr_sig1_list=np.transpose(sig1_list).astype(float)
tr_sig3_list=np.transpose(sig3_list).astype(float)

tr_sig3_data=np.transpose(sig3_data).astype(float)

limitor=np.where(tr_sig3_data[0]<0,tr_sig3_data[0][0],tr_sig3_data[0])

ax_marx.set_ylim(min(limitor)*0.9,max(limitor)*1.1)
                      
# for i in range(nbins):
#     if tr_sig1_list[0][i]<=0:
#         tr_sig1_list[0][i]=1e-10
#     if tr_sig1_list[1][i]<0:
#          tr_sig1_list[1][i]=1e-10
         
#     if tr_sig3_list[0][i]<=0:
#         tr_sig3_list[0][i]=1e-10
#     if tr_sig3_list[1][i]<=0:
#          tr_sig3_list[1][i]=1e-10
        
x_axis=np.linspace(max_factor/nbins,max_factor,num=nbins)

ax_marx.fill_between(x_axis, tr_sig3_list[0],tr_sig1_list[0],
                 color='blue',alpha=0.3,label=r'3 $\sigma$ marx simulations interval')

ax_marx.fill_between(x_axis, tr_sig1_list[1],tr_sig3_list[1],
                 color='blue',alpha=0.3,label='')

ax_marx.fill_between(x_axis, tr_sig1_list[0],tr_sig1_list[1],
                 color='blue',alpha=0.6,label=r'1 $\sigma$ marx simulations interval')

ax_marx.errorbar(x_axis,sbdata,yerr=3*sbdata_err,xerr=None,
             color='black',label=r'3 $\sigma$ data measurements')

ax_marx.errorbar(x_axis,sbdata,yerr=sbdata_err,xerr=None,
             color='red',label=r'1 $\sigma$ data measurements')

ax_marx.legend()

def ratiotorad(x):
    return x * largestsize.value

def radtoratio(x):
    return x / largestsize.value

sec_ax=ax_marx.secondary_xaxis('top',functions=(ratiotorad,radtoratio))
sec_ax.set_xlabel('angle (")')

'''ks test'''
from scipy.stats import ks_2samp

med_list=np.array([ord_sblist[i][int(np.floor(nsim/2))] for i in range(nbins)])

ks, p = ks_2samp(med_list,sbdata)

print('Result of the ks test between the data and the median marx simulation :')
print('ks =',ks)
print('p value for both samples coming from the same distribution =',p)
