# @Author: Andrés Gúrpide <agurpide>
# @Date:   20-05-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 15-06-2021
# Script to compute extinction map

from dust_extinction.parameter_averages import CCM89, O94, F04
from astropy.io import fits
import argparse
import configparser
import os
import numpy as np
import astropy.units as u
import glob

ap = argparse.ArgumentParser(description='Calculate an extinction map using the extinction curves for a given value of R')
ap.add_argument("-r", "--rootname", nargs='?', help="Root name of the input line files", type=str, default="cleancamel")
ap.add_argument("-c", "--curve", choices=["cardelli, odonnell, fitzpatrick"], help="Extinction curve to use", default="cardelli")
ap.add_argument("-R", "--ratio", help="Ratio of total to selective extinction Av/E(B-V)", type=float, default=3.1, nargs="?")
ap.add_argument("-i", "--intrinsic", help="Intrinsic Balmer decrement ratio", type=float, default=2.85, nargs="?")
ap.add_argument("-d", "--decrement", help="Ratio between halpha and hbeta flux maps", type=str, required=True, nargs=1)
args = ap.parse_args()
halpha = 6563.0 * u.angstrom
hbeta = 4861.0 * u.angstrom

outdir = "deredden"
lines = ['OII3727', 'OII3729', 'HBETA', 'NII6548', 'NII6583', 'SII6716', 'SII6731', 'OIII4959', 'OIII5007', 'HALPHA']


if not os.path.isdir(outdir):
    os.mkdir(outdir)
# read config file
Rv = args.ratio
# Rv = E(B - V)/Av
balmer_dec_map = fits.open(args.decrement[0], extname="IMAGE")
observed_ratios = balmer_dec_map[0].data
# get E(beta-alpha)
color_excess = 2.5 * np.log10(observed_ratios / args.intrinsic)
color_excess_map = fits.PrimaryHDU(data=color_excess, header=balmer_dec_map[0].header)
color_excess_map.writeto("%s/color_excess.fits" % outdir, overwrite=True)
if args.curve == "cardelli":
    extinction_function = CCM89(Rv)
elif args.curve == "odonnell":
    extinction_function = O94(Rv)
elif args.curve == "fitzpatrick":
    extinction_function = F04(Rv)

# E(alpha - beta) / A(V)
curve_color_excess = extinction_function.evaluate(hbeta, Rv) - extinction_function.evaluate(halpha, Rv)
Ebv = color_excess / curve_color_excess
#Ebv = color_excess / color_excess * 0.0940
color_excess_map = fits.PrimaryHDU(data=Ebv, header=balmer_dec_map[0].header)
color_excess_map.writeto("%s/total_color_excess.fits" % outdir, overwrite=True)
print("Color excess map saved to %s/total_color_excess.fits" % outdir)
print("Mean color excess: %.2f" % np.nanmean(Ebv))
for line in lines:
    wavemap = glob.glob('./camel_*/cleaned_images/%s_*[!e]wave_*%s.fits' % (args.rootname, line))
    if len(wavemap) == 0:
        print("Line map for %s line not found" % line)
        continue
    wavemap = wavemap[0]
    fluxmap = glob.glob('./camel_*/cleaned_images/%s_*[!e]flux_*%s.fits' % (args.rootname, line))[0]

    wavelenghts = fits.open(wavemap, extname="IMAGE")[0].data
    fluxes_fits = fits.open(fluxmap, extname="IMAGE")[0]
    fluxes = fluxes_fits.data
    header = fluxes_fits.header
    # deredden fluxes
    for x in range(wavelenghts.shape[0]):  # loop on the x range
        for y in range(wavelenghts.shape[1]):  # loop on the y range
            if np.isnan(Ebv[x, y]):
                fluxes[x, y] = np.nan
            else:
                fluxes[x, y] = fluxes[x, y] / extinction_function.extinguish(wavelenghts[x, y] * u.AA, Ebv=Ebv[x, y])
    dereddened_fits = fits.PrimaryHDU(data=fluxes, header=header)
    dereddened_fits.header['CURVE'] = "%s" % args.curve
    dereddened_fits.header['R_v'] = "%.1f" % args.ratio
    dereddened_fits.header['BALMER_MAP'] = "%s" % os.path.basename(args.decrement[0])
    outfile = fluxmap.replace(".fits", "deredden.fits")
    dereddened_fits.writeto(outfile, overwrite=True)
    python_argument = "%s -R %.1f -i %.1f -d %s" % (__file__, Rv, args.intrinsic, args.decrement[0])
    with open("%s/python_command.txt" % outdir, "w+") as f:
        f.write(python_argument)
    print("Dereddened flux map for line %s and stored it to %s" % (line, outfile))
