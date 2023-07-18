# @Author: Andrés Gúrpide <agurpide>
# @Date:   20-05-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 10-02-2022
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
ap.add_argument("-r", "--rootname", nargs='?', help="Root name of the input line files", type=str, default="clean_camel")
ap.add_argument("-R", "--ratio", help="Ratio of total to selective extinction Av/E(B-V). Default Rv=4.05 from Calzetti", type=float, default=4.05, nargs="?")
ap.add_argument("-i", "--intrinsic", help="Intrinsic Balmer decrement ratio. Default 2.86", type=float, default=2.86, nargs="?")
ap.add_argument("-b", "--hbeta", help="Path to HBETA flux map directory", type=str, required=True, nargs=1)
ap.add_argument("-a", "--halpha", help="Path to HALPHA flux map directory", type=str, required=True, nargs=1)
args = ap.parse_args()
halpha = 6563.0 * u.angstrom
hbeta = 4861.0 * u.angstrom

outdir = "deredden_momcheva"
lines = ['OII3727', 'OII3729', 'HBETA', 'OI6300', 'NII6548', 'NII6583', 'SII6716', 'SII6731', 'OIII4959', 'OIII5007', 'HALPHA']
keyword = "flux" # could be "flux"

if not os.path.isdir(outdir):
    os.mkdir(outdir)
# Rv = E(B - V)/Av
Rv = args.ratio
intrinsic_ratio = args.intrinsic
curve = "calzetti"
extincton_curve = C00(Rv)
print('%s/%s_*[!e]%s_*%s.fits' % (args.halpha[0], args.rootname, keyword, "HALPHA"))
try:
    halpha_file = glob.glob('%s/%s_*[!e]%s_*%s.fits' % (args.halpha[0], args.rootname,  keyword, "HALPHA"))[0]
except IndexError:
    raise IndexError("Halpha %s file not found in %s!" % (keyword, args.halpha[0]))
print("Halpha file")
print(halpha_file)
hbeta_file = glob.glob('%s/%s_*[!e]%s_*%s.fits' % (args.hbeta[0], args.rootname, keyword, "HBETA"))[0]
print("Hbeta file")
print(hbeta_file)
halpha_map = fits.open(halpha_file, extname="IMAGE")
if halpha_map[0].header["WCSAXES"] == 3:
    halpha_map[0].header["WCSAXES"] = 2
hbeta_map = fits.open(hbeta_file, extname="IMAGE")
ext = 0
if hbeta_map[ext].data is None:
    ext = 1

observed_ratios = halpha_map[0].data / hbeta_map[ext].data
# get E(beta-halpha)
color_excess = 2.5 * np.log10(observed_ratios / intrinsic_ratio)
color_excess_map = fits.PrimaryHDU(data=color_excess, header=halpha_map[0].header)
color_excess_map.writeto("%s/decrement_color_excess_%.3f.fits" % (outdir, args.intrinsic), overwrite=True)

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
Ebv[np.where(Ebv < 0)] = 0
color_excess_map = fits.PrimaryHDU(data=Ebv, header=halpha_map[0].header)
out_color_excess = "color_excess_%s_Rv%.1f_i%.2f.fits" % (curve, Rv, args.intrinsic)
color_excess_map.writeto("%s/%s" % (outdir, out_color_excess), overwrite=True)
print("Color excess map saved to %s/%s" % (outdir, out_color_excess))
print("Mean color excess: %.2f+-%.2f" % (np.nanmean(Ebv), np.nanstd(Ebv)))

# compute uncertainties
halpha_efile = glob.glob('%s/%s_*e%s_*%s.fits' % (args.halpha[0], args.rootname, keyword, "HALPHA"))[0]
hbeta_efile = glob.glob('%s/%s_*e%s_*%s.fits' % (args.hbeta[0], args.rootname, keyword, "HBETA"))[0]
if len(halpha_efile) > 0:
    uncertainties = 1
    print("HALPHA error file found: %s" % halpha_efile)
    print("HBETA error file found: %s" % hbeta_efile)
    halpha_emap = fits.open(halpha_efile, extname="IMAGE")
    hbeta_emap = fits.open(hbeta_efile, extname="IMAGE")
    # derivative of log10(x) --> x'/x log10e (derivative of dHa/Hb/dHa = 1/Hb; d(Ha/Hb)/dHb = Ha/Hb^2

    err_ebv = 2.5 * np.sqrt((np.log10(np.e) * (halpha_emap[0].data / halpha_map[0].data)) ** 2 + (np.log10(np.e) * hbeta_emap[0].data / hbeta_map[ext].data)**2) / curve_color_excess # I believe there's a mistake here as it should be dlog10(x) = 1 / xln(10) instead of log10(e)/x which is ln(e)/xln(10_
    # I thought the above was incorrect, but comparing two different calculations I get exactly same results
    ecolor_excess_map = fits.PrimaryHDU(data=err_ebv, header=halpha_map[0].header)
    ecolor_excess_map.writeto("%s/ecolor_excess_%s_Rv%.1f_i%.2f.fits" % (outdir, curve, Rv, args.intrinsic), overwrite=True)
else:
    uncertainties = 0
    print("Warning; halpha error file not found. Uncertainties won't be computed")

for line in lines:
    print("\tDereddening %s line" % line)
    wavemap = glob.glob('./camel_*/cleaned_images/%s_*[!e]wave_*%s.fits' % (args.rootname, line))
    if len(wavemap) == 0:
        print("Line map for %s line not found" % line)
        continue
    wavemap = wavemap[0]

    fluxmap = glob.glob('./camel_*/cleaned_images/%s_*[!e]%s_*%s.fits' % (args.rootname, keyword, line))[0]
    efluxmap = glob.glob('./camel_*/cleaned_images/%s_*e%s_*%s.fits' % (args.rootname, keyword, line))[0]
    print("Using wavelength map: %s" % wavemap)
    print("Using flux map: %s" % fluxmap)
    print("Using eflux map: %s" % efluxmap)
    wavelenghts_fits = fits.open(wavemap, extname="IMAGE")
    ext = 0
    if wavelenghts_fits[ext].data is None:
        ext = 1
    wavelenghts = wavelenghts_fits[ext].data
    fluxes_fits = fits.open(fluxmap, extname="IMAGE")[ext]
    efluxes_fits = fits.open(efluxmap, extname="IMAGE")[0]
    fluxes = fluxes_fits.data
    efluxes = efluxes_fits.data
    # deredden fluxes
    for x in range(wavelenghts.shape[0]):  # loop on the x range
        for y in range(wavelenghts.shape[1]):  # loop on the y range
            if np.isnan(Ebv[x, y]):
                fluxes[x, y] = np.nan
                efluxes[x, y] = np.nan

            else:
                extinction_curve_value = extincton_curve.evaluate(wavelenghts[x, y] * u.angstrom)
                fluxes[x, y] = fluxes[x, y] * 10 ** (0.4 * extinction_curve_value * Ebv[x, y])
                if uncertainties:
                    efluxes[x, y] = np.sqrt((efluxes[x, y] * 10 ** (0.4 * extinction_curve_value * Ebv[x, y])) ** 2 + (err_ebv[x, y] * fluxes[x, y] * 10 ** (0.4 * extinction_curve_value * Ebv[x, y]) * np.log(10**(0.4 * extinction_curve_value)))**2)
    dereddened_fits = fits.PrimaryHDU(data=fluxes, header=fluxes_fits.header)
    dereddened_fits.header['CURVE'] = "%s" % curve
    dereddened_fits.header['R_v'] = "%.1f" % args.ratio
    dereddened_fits.header['Ratio'] = "%.1f" % args.intrinsic
    dereddened_fits.header['COMMENT'] = "Balmer ratio: %s/%s" % (halpha_file, hbeta_file)
    outfile = fluxmap.replace(".fits", "_deredden.fits")
    dereddened_fits.writeto(outfile, overwrite=True)
    edereddened_fits = fits.PrimaryHDU(data=efluxes, header=efluxes_fits.header)
    edereddened_fits.header['CURVE'] = "%s" % curve
    edereddened_fits.header['R_v'] = "%.1f" % args.ratio
    edereddened_fits.header['Ratio'] = "%.1f" % args.intrinsic
    edereddened_fits.header['COMMENT'] = "Balmer ratio: %s/%s" % (halpha_file, hbeta_file)
    if uncertainties:
        eoutfile = efluxmap.replace(".fits", "_deredden.fits")
        edereddened_fits.writeto(eoutfile, overwrite=True)
    python_argument = "%s -R %.1f -i %.3f -a %s -b %s" % (__file__, Rv, args.intrinsic, args.halpha[0], args.hbeta[0])
    with open("%s/python_command.txt" % outdir, "w+") as f:
        f.write(python_argument)
    print("Dereddened flux map for line %s and stored it to %s" % (line, outfile))
