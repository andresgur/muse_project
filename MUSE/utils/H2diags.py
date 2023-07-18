# @Author: Maxime Parra and Andres Gurpide
# @Date:   02-05-2023
# @Email:  maxime.parra@univ-grenoble-alpes.fr

'''
Strong line metallicity and ionization parameter estimations from Pilyugin et al. 2016 and Dors et al. 2017.
Since those methods are only usable on photo-ionized regions, we use the BPT diagnostic of OIII and NII as a limit :
we only consider spaxels in the star formation and intermediate zones
(following the logic of Dors et al. 2017, but with the newer BPT version of Law et al. 2021)
papers :
Pilyugin et al. 2016 : https://academic.oup.com/mnras/article/457/4/3678/2589035
Dors et al. 2017 : doi.org/10.1093/mnras/stw3115
Law et al. 2021 : https://arxiv.org/abs/2011.06012
'''
import os
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Image
import numpy as np
import logging
import argparse
from configparser import ConfigParser, ExtendedInterpolation


def plot_map(filename, outfile):
    image = Image(filename)

    img_figure, ax = plt.subplots(1, subplot_kw={'projection': image.wcs.wcs},figsize=(10,8))

    img_figure.suptitle(filename)

    image.plot(ax=ax, scale='linear', show_xlabel=False, show_ylabel=False, extent=None,
                colorbar='v', cmap="cool", vmin=np.nanpercentile(image.data.data, 0.5),
                vmax=np.nanpercentile(image.data.data, 99.5))
    ax.set_xlabel('Ra', labelpad=0)
    ax.set_ylabel('Dec', labelpad=-2)

    img_figure.savefig("%s.png" % outfile)

#defining the calibration functions. The arguments are ordered according to the lines wavelengths
def OH_R(R2,N2,R3):

    OH_RU=8.589+0.022*np.log10(R3/R2)+0.399*np.log10(N2)+(-0.137+0.164*np.log10(R3/R2)+0.589*np.log10(N2))*np.log10(R2)
    OH_RL=7.932+0.944*np.log10(R3/R2)+0.695*np.log10(N2)+(+0.970-0.291*np.log10(R3/R2)-0.019*np.log10(N2))*np.log10(R2)

    return np.where(np.log10(N2)>=-0.6,OH_RU,OH_RL)

def OH_S(N2,S2,R3):

    OH_SU=8.424+0.030*np.log10(R3/S2)+0.751*np.log10(N2)+(-0.349+0.182*np.log10(R3/S2)+0.508*np.log10(N2))*np.log10(S2)
    OH_SL=8.072+0.789*np.log10(R3/S2)+0.726*np.log10(N2)+(+1.069-0.170*np.log10(R3/S2)+0.022*np.log10(N2))*np.log10(S2)

    return np.where(np.log10(N2)>=-0.6,OH_SU,OH_SL)

def OH_R_2D(R2,N2):

    OH_RU_2D=8.589+0.329*np.log10(N2)+(-0.205+0.549*np.log10(N2))*np.log10(R2)

    return np.where(np.log10(N2)>=-0.6,OH_RU_2D,np.nan)

def OH_S_2D(N2,S2):

    OH_SU_2D=8.445+0.699*np.log10(N2)+(-0.253+0.217*np.log10(N2))*np.log10(S2)

    return np.where(np.log10(N2)>=-0.6,OH_SU_2D,np.nan)

ap = argparse.ArgumentParser(description='Create metallicity map filtering out non-H2 regions (if the BPT diagram is provied)')
ap.add_argument("-o", "--outdir", nargs='?', help="Output dir", default='H2_diags', type=str)
ap.add_argument("-bpt", "--bptmap", nargs='?', help="BPT OIII/SII map file path. If not provided the data won't be filtered", default=None, type=str)
ap.add_argument("--config", nargs='?', help="Config file with line ratio maps", type=str, required=True)
args = ap.parse_args()

outdir=args.outdir

if not os.path.isdir(outdir):
    os.mkdir(outdir)

# read config file
bptcfg = ConfigParser(interpolation=ExtendedInterpolation())
bptcfg.read("%s" %args.config)

lineratio_paths = bptcfg["lineratio_paths"]
# we have N2
if bptcfg.has_option("lineratio_paths", "N_2"):
    ratio_N2 = fits.open(lineratio_paths["N_2"])
    metal_file = ratio_N2.copy()
    ratio_N2 = ratio_N2[0].data
    # we have R3
    if bptcfg.has_option("lineratio_paths", "R_3"):
        ratio_R3 = fits.open(lineratio_paths["R_3"])
        ratio_R3 = ratio_R3[0].data
        #We got everything
        if bptcfg.has_option("lineratio_paths", "R_2"):
            print('Computing using the standard formulae')
            ratio_R2 = fits.open(lineratio_paths["R_2"])
            ratio_R2 = ratio_R2[0].data
            case='R_3D'
            metal_map=OH_R(ratio_R2,ratio_N2,ratio_R3)
        # we have N2 R2 and S2
        elif bptcfg.has_option("lineratio_paths", "S_2"):
             ratio_S2 = fits.open(lineratio_paths["S_2"])
             ratio_S2 = ratio_S2[0].data
             print('R2 ratio unavailable.\n'
                   'Computing using the S2 ratio instead.\n')
             case='S_3D'
             metal_map = OH_S(ratio_N2,ratio_S2,ratio_R3)
    # there is no R3 but there is S2
    elif bptcfg.has_option("lineratio_paths", "S_2"):
             print('Some of the ratios are unavailable.\n'
                   'Attempting a two-dimensional computation for the higher values of the NII ratios.\n'
                   'Using the S relation\n')
             case='S_2D'
             ratio_S2 = fits.open(lineratio_paths["S_2"])
             ratio_S2 = ratio_S2[0].data
             metal_map = OH_S_2D(ratio_N2,ratio_S2)
# there is no N2 but there is R2 and S2
elif bptcfg.has_option("lineratio_paths", "R_2") and bptcfg.has_option("lineratio_paths", "S_2"):
    print('Some of the ratios are unavailable.\n'
       'Attempting a two-dimensional computation for the higher values of the NII ratios.\n'
       'Using the S relation\n')
    ratio_R2 = fits.open(lineratio_paths["R_2"])
    ratio_R2 = ratio_R2[0].data
    ratio_S2 = fits.open(lineratio_paths["S_2"])
    ratio_S2 = ratio_S2[0].data
    case='R_2D'
    metal_map = OH_R_2D(ratio_R2,ratio_S2)
    metal_file = ratio_S2.copy()
else:
    raise ValueError("Too few line ratios, impossible to proceed!")

print("Metallicity computed using the %s method" % case)

#computing the results depending of the case


if args.bptmap is not None:
    bpt_file = args.bptmap
    print("Excluding regions =>1 based on BPT diagram 2 (%s)" % bpt_file)
    bpt_map = fits.open(bpt_file)
    if bpt_map[0].data.shape!=metal_map.shape:
        logging.warning("Can't compare the bpt map and the metal map : Their shapes are different.\n"
              "Using the entire metal map as a result. Be careful.")
    metal_map[bpt_map[0].data>1] = np.nan
    #metal_file[0].data[bpt_map[0].data<1] = np.nan

metal_file[0].data = metal_map

metal_file[0].header['COMMENT'] = "Metal map computed with method %s" % case

outfile = "%s/metal_map_%s" % (outdir, case)

metal_file.writeto("%s.fits" % outfile, overwrite=True)

plot_map('%s.fits' % outfile, outfile=outfile)

'''
This part computes ionization parameter map.
'''
#defining the calibration function. Here, Z is the metallicity normalized to solar metallicity

#there are significant uncertainties here, be careful
Z_star=8.69

def U_ion(Z,S2):
    a_ion=-0.26*Z-1.54
    b_ion=-3.69*Z**2+5.11*Z-5.26
    return a_ion*S2+b_ion

ion_par_file = metal_file.copy()

ion_par_file[0].data =U_ion(metal_map/Z_star,np.log10(ratio_S2))

ion_par_file[0].header['COMMENT'] = "Metal map computed with method %s" % case

outfile = "%s/ion_map_%s" % (outdir, case)

ion_par_file.writeto('%s.fits' % outfile, overwrite=True)

plot_map('%s.fits' % outfile, outfile="%s" % outfile)

print("Results stored to %s" % outdir)
