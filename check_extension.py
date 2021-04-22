# @Author: Andrés Gúrpide <agurpide>
# @Date:   21-04-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 22-04-2021



#%% initialisation
import os
import numpy as np
from regions import read_ds9
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
import glob
import argparse
from scipy.stats import ks_2samp
import re
from pathlib import Path

#marx needs to be installed (https://space.mit.edu/cxc/marx/)
#saotrace needs to be installed (https://cxc.harvard.edu/cal/Hrma/Raytrace/SAOTrace.html)

#Don't forget to initialize heasoft, marx and saotrace in your console before launching the script/spyder
#Don't forget to launch this on your ciao python interpreter

#saotrace might not work if you have conflicts in your perl libraries : initialising saotrace might
#replace one (but not all) of your perl libraries environment variables.
#You can manually setup the other variables in your saotrace starting script

parser = argparse.ArgumentParser(description='Parameters to check extension of a give point like source.')
parser.add_argument("-b", "--bins", type=int, help='Number of radial bins for the PSF/data comparison', default=10, nargs='?')
parser.add_argument("--background", help='Background region file', type=str, nargs=1, required=True)
parser.add_argument("-r", "--region", help='Source region file', type=str, nargs=1, required=True)
parser.add_argument("-s", "--simdir", help='MARX simulation directory', type=str, nargs=1, required=True)
parser.add_argument("-e", "--event_file", help='Chandra event file of the observation', type=str, nargs=1, required=True)
parser.add_argument("-o", "--outdir", help='Output director', type=str, nargs="?", default="psf_comparison")
parser.add_argument("-c", "--core", help='Radius to which exclude the PSF core (in arcsec) for the flux normalization (useful in case of pile up see Lehmann et al. 2005)',
                    type=float, nargs="?", default=0)
args = parser.parse_args()

#starting dir
plt.style.use('/home/agurpide/.config/matplotlib/stylelib/paper.mplstyle')
# source region
regfile = args.region[0]

# PSF simulations dir
simdir = args.simdir[0]

outdir = args.outdir

#background region (will be used later)
bgregfile = args.background[0]

#number of annuli
nbins = args.bins

#maximum annuli size ratio compared to the maximal dimension of the source region
max_factor=2

#chandra pixel to arsec ratio
pta=0.492

#radius affected by pile-up in the initial data in arsecs
# rlim_pileup=1.2
rlim_pileup = args.core

#reading the event file to get some more parameters
evt2 = args.event_file[0]

python_argument = "%s -b %d --background %s -r %s -e %s -s %s -c %.2f" % (__file__, nbins, bgregfile, regfile, evt2, simdir, rlim_pileup)

#reading the region file
source_reg = read_ds9(regfile)[0]

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
if not os.path.isdir(outdir):
    os.mkdir(outdir)

with open("%s/python_command.txt" % outdir, "w+") as f:
    f.write(python_argument)

if type(source_reg).__name__ == 'EllipseSkyRegion':
    largestsize = np.maximum(source_reg.width, source_reg.height)
elif type(source_reg).__name__ == 'CircleSkyRegion':
    largestsize = source_reg.radius

maxsize = largestsize * max_factor

#we will write the command line in pixels to avoid complications
maxsize_pixel = round(maxsize.value / pta, 2)
step = maxsize_pixel / nbins

#First we extract the annuli for the initial file, this time with a background file
#however, chandra can't read directly the background file in our case so we rewrite it in an
#acceptable format

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

print("Extracting radial profile of the source")

#string used in the command line
str_annulus='[bin sky=annulus(' + str(ctr_x)+',' + str(ctr_y)+',0:'\
              + str(maxsize_pixel)+':'+str(step)+')]'

#With this we can extract the entire concentric annuli in one go
# outfile="' + outdir + '/data_profile.fits" did not work so we create it here and move it
os.system('dmextract "' + evt2+str_annulus + '" bkg="' + evt2 + '[bin sky=@' + outdir + '/marx_bg.reg]"' \
         + ' outfile="' + outdir + '/data_profile.fits" opt=generic clobber=yes')

#after which we store the surface brillance values and associated 1 sigma errors in arrays
with fits.open('%s/data_profile.fits' % outdir) as data_profile:
    sbdata = data_profile[1].data['SUR_BRI']
    sbdata_err = data_profile[1].data['SUR_BRI_ERR']

'''computations for the renormalization'''

#conversion of the renormalisation annuli's inner and outer radius
rlim_pixel = round(rlim_pileup / pta, 2)
largestsize_pixel = round(largestsize.value, 2) / pta

print("Extracting radial profile of the source excising the inner %.2f arcsec" % rlim_pileup)

#second version of the command line string
str_annulus_renorm = '[bin sky=annulus(' + str(ctr_x) + ',' + str(ctr_y) + ',' + str(rlim_pixel)+':'\
              + str(largestsize_pixel) + ':' + str(largestsize_pixel - rlim_pixel) + ')]'

#extracting the annulus for the data
os.system('dmextract "' + evt2 + str_annulus_renorm + '" bkg="' + evt2 + '[bin sky=@' + outdir + '/marx_bg.reg]"'\
          + ' outfile=' + outdir + '/curr_profile.fits opt=generic clobber=yes bkgerr=gehrels')

# now we can get the renormalization data value
with fits.open('%s/curr_profile.fits' % outdir) as curr_profile:
    sbdata_norm = curr_profile[1].data['SUR_BRI'][0]


radii = np.linspace(0, maxsize, num=nbins + 1)

ra_tuple = (int(ctr_coord.ra.hms[0]), int(ctr_coord.ra.hms[1]), str(round(ctr_coord.ra.hms[2], 4)))
dec_tuple = (int(ctr_coord.dec.dms[0]), abs(int(ctr_coord.dec.dms[1])), str(abs(round(ctr_coord.dec.dms[2], 4))))

dec_init_tuple = (int(ctr_coord.dec.dms[0]), int(abs(ctr_coord.dec.dms[1])), \
                str(abs(round(ctr_coord.dec.dms[2], 4))), str(round(radii[1].value, 4)))

f_reg = open('annuli.reg', 'w+')

f_reg.write(
'# Region file format: DS9 version 4.1\n'\
+'#global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1'\
+' highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\nfk5\n')

for i in range(nbins):

    if i == 0:
        curr_tuple = ra_tuple + dec_init_tuple
        f_reg.write('circle(%i:%i:%s,%i:%i:%s,%s")\n' % curr_tuple)
    else:
        curr_tuple = ra_tuple + dec_tuple + (str(round(radii[i].value, 4)), str(round(radii[i + 1].value, 4)))
        f_reg.write('annulus(%i:%i:%s,%i:%i:%s,%s",%s")\n' % curr_tuple)

f_reg.close()

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

    with fits.open("%s/psf_profile_%d_normed.fits" % (outdir, i)) as curr_profile:
        sblist_norm[i - 1] = curr_profile[1].data['SUR_BRI'][0]

    os.system('dmextract "%s/marxsim_%d.fits%s" outfile=%s/psf_profile_%d.fits opt=generic clobber=yes' % (simdir, i, str_annulus, outdir, i))
    with fits.open('%s/psf_profile_%d.fits' % (outdir, i)) as curr_profile:
        sblist[i - 1] = curr_profile[1].data['SUR_BRI'] * sbdata_norm / sblist_norm[i - 1]
        sblist_err[i - 1] = curr_profile[1].data['SUR_BRI_ERR'] * sbdata_norm / sblist_norm[i - 1]

#%% graphs
fig_marx, ax_marx = plt.subplots()

ax_marx.set_xlabel('Radius (arcsec)')
ax_marx.set_ylabel(r'Net surface brightness (cts/s/cm$^2$)')
ax_marx.set_yscale('log')

x_axis = np.linspace(max_factor / nbins, max_factor, num=nbins)
p = np.percentile(sblist, [50, 99.73, 100 - 99.73], axis=0)

ax_marx.errorbar(radii[1:].value, p[0], yerr=[p[1] - p[0], p[0] - p[2]], color='black', label=r'MARX PSF')
#ax_marx.plot(x_axis, p[0] - p[1], color='black', label=r'MARX PSF')
ax_marx.errorbar(radii[1:].value, sbdata, yerr=sbdata_err,
                 color='red', label=r'Observed')

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
