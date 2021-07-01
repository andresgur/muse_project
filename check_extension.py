# @Author: Andrés Gúrpide <agurpide>
# @Date:   21-04-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   mparra
# @Last modified time: 01-07-2021

'''

Marx PSF error : 3 sigma
data error : 1 sigma

We  divide our regions in nbins (default 10) annuli ranging from the source center to
twice the maximal dimension of the source region (which should be an ellipse or a circle in sky region)

The "core" argument represents the arsec radius of the central part of the data supposedly affected by pileup.
At the end of the script, we renormalize the outer parts of the Marx PSF (outside of this limit) to the data, to reduce errors.
Moreover, the PSF comparison and the graphs only considers radiuses beyond the core limit.

The core limit has to result from manual computations of the proportion of the initial data affected by pile-up.

We also use our own computation of the source center instead of the source centroid
-> You'll need to change a few things to use ellipsoid annuli

Moreover, the code expects a background region containing a single region (in sky coordinates).
If you use composite background regions you'll need to adapt the very last part

If used in read only, the script gives the arguments of the simulation folder and simply plots the comparison.
'''

#%% initialisation
import os,sys
import numpy as np
import logging
from regions import read_ds9
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import glob
import argparse
from scipy.stats import ks_2samp


def write_annuli(ctr_coord, radii):

    out_reg = '# Region file format: DS9 version 4.1\n'\
        + '#global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1'\
        + ' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n'

    out_reg += 'annulus(%d:%d:%.4f,%d:%d:%.4f,' % (ctr_coord.ra.hms[0], ctr_coord.ra.hms[1],
                      ctr_coord.ra.hms[2], ctr_coord.dec.dms[0], np.abs(ctr_coord.dec.dms[1]), np.abs(ctr_coord.dec.dms[2]))

    for radius in radii[1:-1]:
        out_reg += '%.3f",' % radius
    out_reg += '%.3f")' % radii[-1]
    f_reg = open('annuli.reg', 'w+')
    f_reg.write(out_reg)
    f_reg.close()

#marx needs to be installed (https://space.mit.edu/cxc/marx/)
#saotrace needs to be installed (https://cxc.harvard.edu/cal/Hrma/Raytrace/SAOTrace.html)

#Don't forget to initialize heasoft, marx and saotrace in your console before launching the script/spyder
#Don't forget to launch this on your ciao python interpreter

#saotrace might not work if you have conflicts in your perl libraries : initialising saotrace might
#replace one (but not all) of your perl libraries environment variables.
#You can manually setup the other variables in your saotrace starting script

parser = argparse.ArgumentParser(description='Check if a point like source is extended by comparing the data with PSF simulations.')
parser.add_argument("-b", "--bins", type=int, help='Number of radial bins for the PSF/data comparison', default=10, nargs='?')
parser.add_argument("-r", "--region", help='Source region file', type=str, nargs=1, required=True)
parser.add_argument("-bg", "--background", help='Background region file', type=str, nargs=1, required=False)
parser.add_argument("-e", "--energy_range", help='energy range considered for the simulation', type=str, nargs=1,
                    default='0.3-10')
parser.add_argument("-o", "--outdir", help='Output director', type=str, nargs="?", default="psf_comparison")
parser.add_argument("-f", "--factor", help='Factor to which the PSF computation will be performed. f * the extent of the input region. (2 by default)', type=float, default=2)
parser.add_argument("-c", "--core", help='Radius to which exclude the PSF core (in arcsec) for the flux normalization (useful in case of pile up see Lehmann et al. 2005)',
                    type=float, nargs="?", default=0)
parser.add_argument("-p", "--plot_only", help='plot only mode to read previously done simulations',type=bool, nargs=1,required=False)

args = parser.parse_args()

#starting dir
# plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')

#e_min and e_max come from the simdir, but we add the old 0.3-10 interval for compatibility with the older version 
e_range=args.energy_range[0]

if e_range=='old': 
    e_min=300
    e_max=10000
    simdir='marx_outs'
else:
   e_min=str(int(float(e_range[:e_range.find('-')])*1000))
   e_max=str(int(float(e_range[e_range.find('-')+1:])*1000))
   simdir='marx_outs'+'_'+e_range
   
outdir = args.outdir+'_'+e_range


# source region
regfile = args.region[0]

#reading the region file
source_reg = read_ds9(regfile)[0]

#defining the biggest dimension of the source region depending on its shape
if type(source_reg).__name__ == 'EllipseSkyRegion':
    largestsize = np.maximum(source_reg.width, source_reg.height)
elif type(source_reg).__name__ == 'CircleSkyRegion':
    largestsize = source_reg.radius
    
#now we do different things if in read-only mode or not
if args.plot_only is True:
    
    os.chdir(outdir)
    
    #getting the arguments back from the command line file
    with open("%s/python_command.txt" % outdir, "r") as file_arg:
        arg_str=file_arg.read()
        nbins=arg_str[arg_str.find('-b')+3:arg_str.find('-r')-2]
        regfile=arg_str[arg_str.find('-r')+3:arg_str.find('-bg')-2]
        bgregfile=arg_str[arg_str.find('-bg')+4:arg_str.find('-e')-2]
        e_range=arg_str[arg_str.find('-e')+3:arg_str.find('-o')-2]
        max_factor=arg_str[arg_str.find('-f')+3:arg_str.find('-c')-2]
        rlim_pileup=arg_str[arg_str.find('-c')+3:]
    
    #as well as the number of simulations
    psf_files = sorted(glob.glob("psf_profile*.fits"))
    nsim=len(psf_files)/2
    
    #printing everything
    print('Simulation arguments :\nnumber of bins: %s \nregion file: %s \nbackground region file: %s \nenergy range : %s keV\n\
          maximum psf factor: %.2f \nexcluded core radius: %.2f\nnumber of simulations: %s"'% (nbins, regfile, bgregfile, e_range, max_factor, rlim_pileup,nsim))
    
    #switching back datatypes to avoid problems later
    nbins=int(nbins)
    max_factor=float(max_factor)
    
    #reading the data and psf files for the surface brillance
    sblist = np.empty((nsim, nbins))
    sblist_err = np.empty((nsim, nbins))
    sblist_norm = np.empty((nsim, nbins))
    
    with fits.open('data_profile.fits') as data_profile:
        sbdata = data_profile[1].data['SUR_BRI']
        sbdata_err = data_profile[1].data['SUR_BRI_ERR']
        if bgregfile != "None":
            bgbri=data_profile[1].data['BG_SUR_BRI']

    with fits.open('datanorm_profile.fits') as datanorm_profile:
        sbdata_norm = datanorm_profile[1].data['SUR_BRI'][0]
    
    for i in range(1, nsim + 1):
   
        with fits.open("psf_profile_%d_normed.fits" % (i)) as marxnorm_profile:
            sblist_norm[i - 1] = marxnorm_profile[1].data['SUR_BRI'][0]
    
        with fits.open('psf_profile_%d.fits' % (i)) as marx_profile:
            sblist[i - 1] = marx_profile[1].data['SUR_BRI'] * sbdata_norm / sblist_norm[i - 1]
            sblist_err[i - 1] = marx_profile[1].data['SUR_BRI_ERR'] * sbdata_norm / sblist_norm[i - 1]
            
    #creating the rest of the parameters used for the graphs
    maxsize = largestsize * max_factor
    radii = np.linspace(0, maxsize, num=nbins + 1)

else:
    
    if args.region is None :
        logging.error('Error : missing region argument')
        sys.exit()
        
    #number of annuli
    nbins = args.bins
    
    #maximum annuli size ratio compared to the maximal dimension of the source region
    max_factor = args.factor
    
    #chandra pixel to arsec ratio
    pta = 0.492
    
    #radius affected by pile-up in the initial data in arsecs
    rlim_pileup = args.core
    
    #reading the event file to get some more parameters
    evt2 = glob.glob('./*_evt2.fits')[0]
    
    if args.background is not None:
        #background region (will be used later)
        bgregfile = args.background[0]
    else:
        bgregfile = "None"
    
    
    #sky coordinates of the region file centroid, will be used as position for the marx simulation
    source_ra=str(source_reg.center.ra.deg)
    source_dec=str(source_reg.center.dec.deg)
    
    
    #file reprocessing to recreate evenly binned event files and arf to get the best results possible
    os.system('dmstat "'+evt2+'[sky=region('+regfile+')][bin sky=1]" centroid=yes')
    
    #computing the source center
    stdout=Popen('pget dmstat out_cntrd_phys', shell=True, stdout=PIPE).stdout
    ctr_output=str(stdout.read())
    
    #those are the estimated coordinates of the source center, based on the brightest pixel in the source region
    ctr_x=ctr_output[2:ctr_output.find(',')]
    ctr_y=ctr_output[ctr_output.find(',') + 1:-3]
    
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
    ctr_coord=SkyCoord(ctr_ra,ctr_dec, unit=(u.hourangle, u.deg), frame='fk5')
    ctr_ra=str(ctr_coord.ra.deg)
    ctr_dec=str(ctr_coord.dec.deg)
    
    #creating the output directory if it doesn't exist already
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    #we start by putting a copy of the command to launch this iteration of the script inside
    python_argument = "%s -b %d -r %s -bg %s -e %s -o %s -f %.2f -c %.2f" % (__file__, nbins, regfile, bgregfile, e_range, outdir, max_factor, rlim_pileup)
    
    with open("%s/python_command.txt" % outdir, "w+") as f:
        f.write(python_argument)
    
    #see largestsize before the loop
    maxsize = largestsize * max_factor
    
    #we will write the command line in pixels to avoid complications
    maxsize_pixel = round(maxsize.value / pta, 2)
    step = maxsize_pixel / nbins
    
    #First we extract the annuli for the initial file, this time with a background file
    #however, chandra can't read directly the background file in our case so we rewrite it in an
    #acceptable format
    
    print("Extracting radial profile of the source up to %.2f %s" % (maxsize.value, maxsize.unit))
    
    #string used in the command line
    str_annulus='[bin sky=annulus(' + str(ctr_x)+',' + str(ctr_y)+',0:'\
                  + str(maxsize_pixel)+':'+str(step)+')]'
    
    #string used for the energy 
    str_ener='[energy='+e_min+':'+str(float(e_max)-0.001)+']'
    
    #With this we can extract the entire concentric annuli in one go
    # outfile="' + outdir + '/data_profile.fits" did not work so we create it here and move it
    
    if args.background is not None:
        f_bg = open(bgregfile, 'r')
        bg_lines = f_bg.readlines()
    
        #sometimes ds9 doesn't bother including the format (fk5 for example) when saving regions, so we test
        #the length of the third line. If it's long, it's probably the region line instead of a format
        if len(bg_lines[2]) < 20:
            bgreg_line = bg_lines[3]
        else:
            bgreg_line = bg_lines[2]
    
        f_bg_marx = open('%s/marx_bg.reg' % outdir, 'w+')
        f_bg_marx.write(bgreg_line[:bgreg_line.find('#')])
        f_bg_marx.close()
        print("Extracting source brightness profile WITH background subtraction.")
        os.system('dmextract "' + evt2 + str_ener + str_annulus + '" bkg="' + evt2 + '[bin sky=@' + outdir + '/marx_bg.reg]"' \
             + ' outfile="' + outdir + '/data_profile.fits" opt=generic clobber=yes')
    else:
        print("Extracting source brightness profile WITHOUT background subtraction.")
        os.system('dmextract "' + evt2 + str_ener + str_annulus + '" outfile="' + outdir + '/data_profile.fits" opt=generic clobber=yes')
    
    #after which we store the surface brillance values and associated 1 sigma errors in arrays
    with fits.open('%s/data_profile.fits' % outdir) as data_profile:
        sbdata = data_profile[1].data['SUR_BRI']
        sbdata_err = data_profile[1].data['SUR_BRI_ERR']
        if args.background is not None:
            bgbri=data_profile[1].data['BG_SUR_BRI']
    
    '''computations for the renormalization'''
    
    #conversion of the renormalisation annuli's inner and outer radius
    rlim_pixel = round(rlim_pileup / pta, 2)
    largestsize_pixel = round(largestsize.value, 2) / pta
    
    print("Extracting radial profile of the source excising the inner %.2f arcsec" % rlim_pileup)
    
    #second version of the command line string (for marx, take off the central part of the PSF)
    str_annulus_renorm = '[bin sky=annulus(' + str(ctr_x) + ',' + str(ctr_y) + ',' + str(rlim_pixel)+':'\
                  + str(largestsize_pixel) + ':' + str(largestsize_pixel - rlim_pixel) + ')]'
                  
    if args.background is not None:
    
        #extracting the annulus for the data
        os.system('dmextract "' + evt2 + str_ener + str_annulus_renorm + '" bkg="' + evt2 + '[bin sky=@' + outdir + '/marx_bg.reg]"'\
                  + ' outfile=' + outdir + '/datanorm_profile.fits opt=generic clobber=yes bkgerr=gehrels')
    else:
        os.system('dmextract "' + evt2 + str_ener + str_annulus_renorm + ' outfile=' + outdir + '/datanorm_profile.fits opt=generic clobber=yes bkgerr=gehrels')
    
    # now we can get the renormalization data value
    with fits.open('%s/datanorm_profile.fits' % outdir) as datanorm_profile:
        sbdata_norm = datanorm_profile[1].data['SUR_BRI'][0]
    
    
    radii = np.linspace(0, maxsize, num=nbins + 1)
    
    write_annuli(ctr_coord, radii.value)
    
    os.system('punlearn dmextract')
    
    #Now, we extract the renormalization annuli and the concentric annuli in a similar manner than with the data,
    #adding the renormalization
    simfiles = sorted(glob.glob("%s/marxsim_*.fits" % simdir))
    # simfiles won't be sorted so better not go iterate through them
    nsim = len(simfiles)
    sblist = np.empty((nsim, nbins))
    sblist_err = np.empty((nsim, nbins))
    sblist_norm = np.empty((nsim, nbins))
    
    for i in range(1, nsim + 1):
        os.system('dmextract "%s/marxsim_%d.fits%s" outfile=%s/psf_profile_%d_normed.fits opt=generic clobber=yes' % (simdir, i, str_annulus_renorm, outdir, i))
    
        with fits.open("%s/psf_profile_%d_normed.fits" % (outdir, i)) as marxnorm_profile:
            sblist_norm[i - 1] = marxnorm_profile[1].data['SUR_BRI'][0]
    
        os.system('dmextract "%s/marxsim_%d.fits%s" outfile=%s/psf_profile_%d.fits opt=generic clobber=yes' % (simdir, i, str_annulus, outdir, i))
        with fits.open('%s/psf_profile_%d.fits' % (outdir, i)) as marx_profile:
            sblist[i - 1] = marx_profile[1].data['SUR_BRI'] * sbdata_norm / sblist_norm[i - 1]
            sblist_err[i - 1] = marx_profile[1].data['SUR_BRI_ERR'] * sbdata_norm / sblist_norm[i - 1]

#%% graphs
fig_marx, ax_marx = plt.subplots()

fig_marx.suptitle('PSF extension comparison in the '+str(e_range)+' keV range, with '+str(nsim)+' simulations')
ax_marx.set_xlabel('Radius (arcsec)')
if args.background is not None:
    ax_marx.set_ylabel(r'Net surface brightness (cts/s/cm$^2$)')
else:
    ax_marx.set_ylabel(r'Surface brightness (cts/s/cm$^2$)')
ax_marx.set_yscale('log')

x_axis = np.linspace(max_factor / nbins, max_factor, num=nbins)
p = np.percentile(sblist, [50, 99.73, 100 - 99.73], axis=0)

ax_marx.errorbar(radii[1:].value, p[0], yerr=[p[1] - p[0], p[0] - p[2]], color='black', label=r'MARX PSF')
#ax_marx.plot(x_axis, p[0] - p[1], color='black', label=r'MARX PSF')
ax_marx.errorbar(radii[1:].value, sbdata, yerr=sbdata_err,
                 color='red', label=r'Observed')

if bgregfile != "None":
    ax_marx.axhline(y=bgbri[0],xmin=0,xmax=1,color='black',alpha=0.5,linestyle='--',linewidth=2,label='Background level')

# peform ks test
ks, p = ks_2samp(p[0], sbdata)

print('Result of the ks test between the data and the median marx simulation :')
print('ks = %.2f' % ks)
print('p value for both samples coming from the same distribution = %.2f' % p)
with open("%s/ks_text.txt" % outdir, "w+") as f:
    f.write("#ks\tp-value\n%.3f\t%.3f" % (ks, p))
print("Written results to %s" % outdir)

ax_marx.legend()

fig_marx.savefig("%s/psf_comparison_%d.png" % (outdir, nsim))


def ratiotorad(x):
    return x * largestsize.value


def radtoratio(x):
    return x / largestsize.value

#sec_ax=ax_marx.secondary_xaxis('top',functions=(ratiotorad,radtoratio))
#sec_ax.set_xlabel('Radius (arcsec)')
