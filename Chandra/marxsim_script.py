# @Author: Maxime Parra and Andrés Gúrpide <agurpide>
# @Date:   21-04-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   mparra
# @Last modified time: 29-06-2021



#%% initialisation
import os, shutil
import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from subprocess import Popen, PIPE
import argparse
import glob
from regions import read_ds9
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

Format for the energy range argument : min-max

Please input real region files and not links
No background included as of now in the simulation parameters (it is included in the comparison after)
(why are you even checking your psf if your source isn't bright enough to stand out)
More  details about the comparison at the end of the script


IF marxpileup.par is not found
try to use "locate marxpileup.par" and delete it if it is empty
'''

''' PARAMETERS '''
parser = argparse.ArgumentParser(description='Parameters for the MARX simulattions. To be launched from the directory where all Chandra data is stored (e.g. repro dir)')
parser.add_argument("--sao", help='Flag to use sao trace in the simulations (only when high degree of accuracy is needed, see https://cxc.harvard.edu/ciao/PSFs/raytracers.html)', action="store_true")
parser.add_argument("-r", "--region", help='Source region file', type=str, nargs=1, required=True)
parser.add_argument("-b", "--background", help='Background region file', type=str, nargs=1, required=False)
parser.add_argument("-n", "--nsimulations", help='Number of MARX simulations to perform', type=int, nargs='?', default=500)
parser.add_argument("-e", "--energy_range", help='energy range considered for the simulation', type=str, nargs=1,
                    default='0.3-10')
args = parser.parse_args()
sao_arg = args.sao

# number of simulations
nsim = args.nsimulations

e_range=args.energy_range[0]
e_min=e_range[:e_range.find('-')]
e_max=e_range[e_range.find('-')+1:]

#clean the ardlib parameter files (notably containing paths to bad pixel files)
os.system('punlearn ardlib')

# input region
regfile = args.region[0]


# background region
if args.background is not None:
    bgregfile = args.background[0]
else:
    bgregfile = "None"
    
python_argument = "%s -n %d --background %s -r %s --sao %s" % (__file__, nsim, bgregfile, regfile, args.sao)

# chandra ACIS pixel to arsec ratio
pta = 0.492

'''############'''

if sao_arg:
    simdir = 'sao_sim'
    fitsname = 'saosim_'
    outdir = "sao_outs"
else:
    simdir = 'marx_sim'
    fitsname = "marxsim_"
    outdir = "marx_outs"

#ditherfile, there should only be one ditherfile for an observation
#os.chdir("%s" % sdir)

ditherfiles = glob.glob('./*_asol1.fits')
if len(ditherfiles) > 1:
    ditherfile = "full_asol.fits"
    print("%d aspect solution files found. These files will be merged into a single aspect solution file named %s" % (len(ditherfiles), ditherfile))
    # join them after sorting chronologically
    ditherfile_string = ",".join(np.sort(ditherfiles))
    # merge them into a single file https://cxc.cfa.harvard.edu/ciao/why/asol.html
    os.system('dmmerge "%s" %s clobber=yes' % (ditherfile_string, ditherfile))

else:
    ditherfile = ditherfiles[0]
# in case there are multiple dither files, merge them onto a single string (https://cxc.cfa.harvard.edu/ciao/why/asol.html)
#ditherfile_string = "pcadf181476117N003_asol1.fits, pcadf181476245N003_asol1.fits"
evt2 = glob.glob('./*_evt2.fits')[0]

#reading the data file, which is assumed to have the same name as the region file
name = regfile[:regfile.find('.')]
data_file = name + '.pi'
# IF THE USER IS NOT DUMB WE COULD DO data_file = regfile.replace(".reg", ".pi")
#reading the region file
source_reg = read_ds9(regfile)[0]
# sky coordinates of the region file centroid, will be used as position for the marx simulation
source_ra = str(source_reg.center.ra.deg)
source_dec = str(source_reg.center.dec.deg)

#reading the event file to get some more parameters
fits_evt2 = fits.open(evt2)

#instrument used
instrum = fits_evt2[0].header['INSTRUME']

#coordinates of nominal aimpoint
ra_nom=str(fits_evt2[1].header["RA_NOM"])
dec_nom=str(fits_evt2[1].header["DEC_NOM"])
roll_nom=str(fits_evt2[1].header["ROLL_NOM"])

#date of observation
date = str(round(Time(fits_evt2[1].header["DATE-OBS"]).decimalyear, 3))

#grating of observation
grating = fits_evt2[1].header["GRATING"]

#file reprocessing to recreate evenly binned event files and arf to get the best results possible
os.system('punlearn dmcoords')
os.system('dmstat "'+evt2+'[sky=region('+regfile+')][bin sky=1]" centroid=yes')

#computing the source center

stdout = Popen('pget dmstat out_cntrd_phys', shell=True, stdout=PIPE).stdout
ctr_output=str(stdout.read())

#those are the estimated coordinates of the source center, based on the brightest pixel in the source region
ctr_x = ctr_output[2:ctr_output.find(',')]
ctr_y = ctr_output[ctr_output.find(',')+1:-3]

# computing some more informations about the source center
os.system('punlearn dmcoords')
os.system('dmcoords '+evt2+' op=sky x='+ctr_x+' y='+ctr_y)
stdout = Popen('pget dmcoords chip_id', shell=True, stdout=PIPE).stdout
chip_output = str(stdout.read())
chip_id = int(chip_output[2:-3])

# such as the center position in sky coordinates
stdout = Popen('pget dmcoords ra dec', shell=True, stdout=PIPE).stdout
coord_output = str(stdout.read())[2:-3]

ctr_ra = coord_output[:coord_output.find('\\n')]
ctr_dec = coord_output[coord_output.find('\\n') + 2:]
ctr_coord = SkyCoord(ctr_ra, ctr_dec, unit=(u.hourangle, u.deg), frame='fk5')
ctr_ra = str(ctr_coord.ra.deg)
ctr_dec = str(ctr_coord.dec.deg)

# We now also know which part of the instrument is used for the detection
if instrum == 'ACIS':
    if chip_id < 4:
        detector = 'ACIS-I'
        subsys = 'ACIS-I' + str(chip_id)
    elif chip_id > 3:
        detector = 'ACIS-S'
        subsys = 'ACIS-S' + str(chip_id - 4)
print("Detector used %s" % detector)

if instrum == 'HRC':
    detector=input('Since the CCD number overlap in this instrument, please confirm'\
                   +' if the detecor used is HRC-I or HRC-S')

os.system('asphist infile="%s" outfile="marx.asp" evtfile=%s clobber=yes' % (ditherfile, evt2))

#here we only use mkarf to create the response file, a more accurate version can be computed with mkwarf
#but is much more complex to create

binning = 0.1

os.system('mkarf detsubsys="'+subsys+'" grating="'+grating+'" outfile="marx.arf" obsfile="'+evt2
          +'" asphistfile="marx.asp" sourcepixelx='+ctr_x+' sourcepixely='+ctr_y+' engrid="'+e_min+':'+e_max
          +':%.1f" maskfile=NONE pbkfile=NONE dafile=NONE verbose=1 mode=h clobber=yes' % binning)

#not done in the marx tutorial but should be done for a point source...right ?
print("Creating ARF corrected file for point like source")
os.system('arfcorr infile="'+evt2+'[sky=region('+regfile+')][bin sky]" region="region('+regfile+\
          ')" x='+ctr_x+' y='+ctr_y+' arf=marx.arf outfile=marx_corr.arf clobber=yes')

#if you want to create an image file to check if the detection is correctly centered
# os.system('dmcopy "'+evt2+'"[EVENTS][bin x='+str(int(float(ctr_x)-100))+':'+str(int(float(ctr_x)+100))\
#           +':1,y='+str(int(float(ctr_y)-100))+':'+str(int(float(ctr_y)+100))\
#           +':1" image_src.fits option=image clobber=yes')


#source energy computation
os.system('dmtcalc "' + evt2 + '[EVENTS][sky=region('+regfile+')]" marx_evt2_energy.fits'\
          +' expr="energy=(float)pi*0.0149" clobber=yes')
    
if args.background is not None:
    #bg energy computation. Should be negligible, but included nonetheless
    print("Creating energy histogram WITH background subtraction")
    os.system('dmtcalc "'+evt2+'[EVENTS][sky=region('+bgregfile+')]" marx_bg_energy.fits'\
              +' expr="energy=(float)pi*0.0149" clobber=yes')

    # net histogram computation with real energy binning, only in the specified energy range
    os.system('dmextract "marx_evt2_energy.fits[bin energy='+e_min+':'+str(float(e_max)-0.001)+':%.1f]"'
              +' bkg=marx_bg_energy.fits outfile=marx_extract.spec clobber=yes opt=generic' % binning)
else:

    print("Creating energy histogram WITHOUT background subtraction")
    # source histogram computation with real energy binning, only in the specified energy range
    os.system('dmextract "marx_evt2_energy.fits[bin energy='+e_min+':'+str(float(e_max)-0.001)+':%.1f]"'
              +' outfile=marx_extract.spec clobber=yes opt=generic' % binning)

# reading the newly created files
arf_corr = fits.open('marx_corr.arf')
arf_data = arf_corr[1].data
extract = fits.open('marx_extract.spec')
extract_data = extract[1].data

# obtaining the exposure time
expos = "%.1f" % extract[1].header['EXPOSURE']
if args.background is not None:
    rate_keyword = "NET_RATE"
else:
    rate_keyword = "COUNT_RATE"
    
'''
Now we create the spectrum file which will serve as an argument for marx/sao
For Marx the (received) flux column is supposed to be a normalized flux DENSITY so the net rate column needs
to be divided by the binning, on top of the effective area using the arf and
the correction factor to account for the loss percentage.
'''

if sao_arg:
    print("Running SAO TRACE simulation")

    flux = extract_data[rate_keyword] / (arf_data['SPECRESP'])
    f_sao = open('marx_sao.rdb', 'w+')
    f_sao.write('#argument file for SAOTrace, manually created from marx_extract.spec and marx_corr.arf\n')
    f_sao.write('ENERG_LO\tENERG_HI\tFLUX\n')
    for i in range(len(flux)):

        f_sao.write('%f\t%f\t%f\n' % (arf_data['ENERG_LO'][i], arf_data['ENERG_HI'][i], flux[i]))

    f_sao.close()

    #sao parameters file

    fpar_sao=open('marx_sao_par.lua','w+')
    fpar_sao.write('ra_pnt='+ra_nom+'\ndec_pnt='+dec_nom+'\nroll_pnt='+roll_nom+'\n\n'\
                  +'dither_asol_chandra{file="%s", ra=ra_pnt, dec=dec_pnt, roll=roll_pnt}\n\n'
                  +'point{position={ra='+ctr_ra+', dec='+ctr_dec+',ra_aimpt=ra_pnt, dec_aimpt=dec_pnt,},\n'
                  +'spectrum={{file="marx_sao.rdb", units="photons/s/cm2", scale=1, format="rdb", emin="ENERG_LO",'
                  +'emax="ENERG_HI",flux="FLUX"}}}' % (ditherfile))
    fpar_sao.close()

    #sao not working as of now

    os.system('trace-nest tag=saotrace srcpars=marx_sao_par.lua tstart='+date+' limit='+expos+' limit_type=sec')

    os.system('punlearn marx')

    #marx parameters :
    #SourceFlux=-1 tells marx to take the integrated flux of the Spectrum
    #Verbose=yes to have diagnostic messages

    par_marx='pset marx SourceType="SAOSAC" SourceFlux=-1 SAOSACFile=saotrace.fits'\
              +' ExposureTime='+expos+' TStart='+date+' OutputDir=sao_sim GratingType='+grating\
              +' DetectorType=' + detector + ' RA_Nom='+ra_nom+' Dec_Nom='+dec_nom\
              +' DitherModel="FILE" DitherFile='+ditherfile\
              +' Roll_Nom='+roll_nom+' SourceRA='+source_ra+' SourceDEC='+source_dec+' Verbose=yes mode=h'

else:
    print("Running marx without sao trace")
    fluxdens = extract_data[rate_keyword] / (arf_data['SPECRESP'] * binning)

    f_marx = open('marx_arg.tbl', 'w+')
    f_marx.write('#argument file for Marx, manually created from marx_extract.spec and marx_corr.arf\n')
    for i in range(len(fluxdens)):

        f_marx.write('%f %f\n' % (arf_data['ENERG_HI'][i], fluxdens[i]))

    f_marx.close()

    os.system('punlearn marx')

    #marx parameters :
    #SourceFlux=-1 tells marx to take the integrated flux of the Spectrum
    #Verbose=yes to have diagnostic messages
    # MARX only ACCEPTS 1 single dither file so we need to merge them
    par_marx='pset marx SpectrumType="FILE" SourceFlux=-1 SpectrumFile=marx_arg.tbl'\
              +' ExposureTime='+expos+' TStart='+date+' OutputDir=marx_sim GratingType='+grating\
              +' DetectorType=' + detector +'  RA_Nom='+ra_nom+' Dec_Nom='+dec_nom\
              +' DitherModel="FILE" DitherFile='+ditherfile\
              +' Roll_Nom='+roll_nom+' SourceRA='+source_ra+' SourceDEC='+source_dec+' Verbose=yes mode=h'

os.system(par_marx)

# temporary outputs directory, will be deleted

totoutdir=outdir+'_'+e_range

os.system('mkdir -p ' +totoutdir)
with open("%s/python_command.txt" % totoutdir, "w+") as f:
    f.write(python_argument)
os.system('mkdir -p ' + simdir)

for i in range(1, nsim + 1):
    print("**************************************\n"+
          "***Performing MARX simulation %d/%d***\n"% (i, nsim + 1)+
          "**************************************\n")
    os.system('marx')

    os.chdir(simdir)
    curr_fitsname = '%s%d.fits' % (fitsname, i)

    # running a pileup command to emulate the loss of counts due to pileup
    os.system('marxpileup MarxOutputDir="./"')
    os.chdir('pileup')

    #creating a fits file from the outputs. This is the only output we're interested in
    #so we move it and start again (marx and marxpileup both overwrite files by default)
    os.system('marx2fits ./ '+ curr_fitsname)
    os.system('mv '+ curr_fitsname+' ../../' + totoutdir)
    os.chdir('../..')

#deleting all the temporary files
os.chdir(simdir + '/pileup')
os.system('rm *')
os.chdir('..')
os.system('rmdir pileup')
os.system('rm *')
os.chdir('..')
os.system('rmdir ' + simdir)

#moving the initialisation files
for files in os.listdir('./'):
    if 'marx' in files and not os.path.isdir(files):
        shutil.move(files,totoutdir)

print("Written results to %s" % totoutdir)

print("Run 'python check_extension.py -r %s --background %s -e %s -c %s' to analyse the profile of the PSF" % (regfile, bgregfile, e_range, 'your_pileup_value'))
