# @Author: Andrés Gúrpide <agurpide>
# @Date:   20-05-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 25-06-2021
# Script to compute extinction map
from astropy.io import fits
import argparse
import os
import numpy as np
import astropy.units as u
import glob


class C00:

    def __init__(self, Rv):
        self.Rv = Rv

    def evaluate(self, wavelength):
        """wavelength must be in wavelength units"""
        micron = wavelength.to(u.micron).value
        x = 1 / micron
        optical_indx = np.where(np.logical_and(0.63 <= micron, micron <= 2.20))
        ir_indx = np.where(np.logical_and(0.12 <= micron, micron <= 0.63))
        x = np.asarray(x)
        if x.ndim == 0:
            x = x[None]
        k = np.empty(len(x))
        k[optical_indx] = 2.659 * (-1.857 + 1.040 * x) + self.Rv
        k[ir_indx] = 2.659 * (-2.156 + 1.509 * x - 0.198 * x**2 + 0.011 * x**3) + self.Rv
        return k


ap = argparse.ArgumentParser(description='Calculate an extinction map using the extinction curve of Calzetti for a given value of R')
ap.add_argument("-r", "--rootname", nargs='?', help="Root name of the input line files", type=str, default="cleancamel")
ap.add_argument("-R", "--ratio", help="Ratio of total to selective extinction Av/E(B-V)", type=float, default=4.05, nargs="?")
ap.add_argument("-i", "--intrinsic", help="Intrinsic Balmer decrement ratio", type=float, default=2.86, nargs="?")
ap.add_argument("-b", "--hbeta", help="Path to HBETA flux maps", type=str, required=True, nargs=1)
ap.add_argument("-a", "--halpha", help="Path to HALPHA flux maps", type=str, required=True, nargs=1)
args = ap.parse_args()
halpha = 6563.0 * u.angstrom
hbeta = 4861.0 * u.angstrom

outdir = "deredden_momcheva"
lines = ['OII3727', 'OII3729', 'HBETA', 'NII6548', 'NII6583', 'SII6716', 'SII6731', 'OIII4959', 'OIII5007', 'HALPHA']


if not os.path.isdir(outdir):
    os.mkdir(outdir)
# Rv = E(B - V)/Av
Rv = args.ratio
curve = "calzetti"
extincton_curve = C00(Rv)

halpha_file = glob.glob('%s/*_*[!e]flux_*%s.fits' % (args.halpha[0], "HALPHA"))[0]
hbeta_file = glob.glob('%s/*_*[!e]flux_*%s.fits' % (args.hbeta[0], "HBETA"))[0]
print("HALPHA file found: %s" % halpha_file)
print("HBETA file found: %s" % hbeta_file)
halpha_map = fits.open(halpha_file, extname="IMAGE")
if halpha_map[0].header["WCSAXES"] == 3:
    halpha_map[0].header["WCSAXES"] = 2
hbeta_map = fits.open(hbeta_file, extname="IMAGE")


observed_ratios = halpha_map[0].data / hbeta_map[0].data
# get E(beta-halpha)
color_excess = 2.5 * np.log10(observed_ratios / args.intrinsic)
color_excess_map = fits.PrimaryHDU(data=color_excess, header=halpha_map[0].header)
color_excess_map.writeto("%s/decrement_color_excess.fits" % outdir, overwrite=True)

# E(alpha - beta) / A(V)
kbeta = extincton_curve.evaluate(hbeta)
kalpha = extincton_curve.evaluate(halpha)
curve_color_excess = (kbeta - kalpha)
print("Extinction curve at hbeta")
print(kbeta)
print("Extinction curve at Halpha")
print(kalpha)
print("k(beta) - k(halpha)")
print(curve_color_excess)
Ebv = color_excess / curve_color_excess
color_excess_map = fits.PrimaryHDU(data=Ebv, header=halpha_map[0].header)
color_excess_map.writeto("%s/color_excess_%s_Rv%.1f_i%.2f.fits" % (outdir, curve, Rv, args.intrinsic), overwrite=True)
print("Color excess map saved to %s/total_color_excess.fits" % outdir)
print("Mean color excess: %.2f" % np.nanmean(Ebv))
# compute uncertainties
halpha_efile = glob.glob('%s/*_*eflux_*%s.fits' % (args.halpha[0], "HALPHA"))[0]
hbeta_efile = glob.glob('%s/*_*eflux_*%s.fits' % (args.hbeta[0], "HBETA"))[0]
print("HALPHA error file found: %s" % halpha_efile)
print("HBETA error file found: %s" % hbeta_efile)
halpha_emap = fits.open(halpha_efile, extname="IMAGE")
hbeta_emap = fits.open(hbeta_efile, extname="IMAGE")
# derivative of log10(x) --> x'/x log10e (derivative of dHa/Hb/dHa = 1/Hb; d(Ha/Hb)/dHb = Ha/Hb^2
err_ebv = 2.5 * np.sqrt((np.log10(np.e) * (1 / halpha_map[0].data) * halpha_emap[0].data) ** 2 + (np.log10(np.e) * (1 / hbeta_map[0].data) * hbeta_emap[0].data)**2) / curve_color_excess
ecolor_excess_map = fits.PrimaryHDU(data=err_ebv, header=halpha_map[0].header)

ecolor_excess_map.writeto("%s/ecolor_excess_%s_Rv%.1f_i%.2f.fits" % (outdir, curve, Rv, args.intrinsic), overwrite=True)
for line in lines:
    print("\tDereddening %s line" % line)
    wavemap = glob.glob('./camel_*/cleaned_images/%s_*[!e]wave_*%s.fits' % (args.rootname, line))
    if len(wavemap) == 0:
        print("Line map for %s line not found" % line)
        continue
    wavemap = wavemap[0]
    fluxmap = glob.glob('./camel_*/cleaned_images/%s_*[!e]flux_*%s.fits' % (args.rootname, line))[0]
    print("Using wavelength map: %s" % wavemap)
    print("Using flux map: %s" % fluxmap)
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
                fluxes[x, y] = fluxes[x, y] * 10 ** (0.4 * extincton_curve.evaluate(wavelenghts[x, y] * u.angstrom) * Ebv[x, y])
    dereddened_fits = fits.PrimaryHDU(data=fluxes, header=header)
    dereddened_fits.header['CURVE'] = "%s" % curve
    dereddened_fits.header['R_v'] = "%.1f" % args.ratio
    dereddened_fits.header['COMMENT'] = "Balmer ratio: %s/%s" % (halpha_file, hbeta_file)
    outfile = fluxmap.replace(".fits", "deredden.fits")
    dereddened_fits.writeto(outfile, overwrite=True)
    python_argument = "%s -R %.1f -i %.1f -a %s -b %s" % (__file__, Rv, args.intrinsic, args.halpha[0], args.hbeta[0])
    with open("%s/python_command.txt" % outdir, "w+") as f:
        f.write(python_argument)
    print("Dereddened flux map for line %s and stored it to %s" % (line, outfile))
