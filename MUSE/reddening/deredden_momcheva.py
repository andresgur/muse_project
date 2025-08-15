# @Author: Andrés Gúrpide <agurpide>
# @Date:   20-05-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 10-02-2025
from astropy.io import fits
import argparse
import os
import numpy as np
import astropy.units as u
import glob
from deredden_utils import galactic_extinction, C00, run_get_ebv_script, parse_ebv_output


def parse_ebv_arg(EBV_arg, file):

    if EBV_arg is None:
        try:
            fitsfile = fits.open(file, extname="IMAGE")
        except:
            fitsfile = fits.open(file, extname=1)
        if "CRVAL1" in fitsfile[0].header:
            ra = fitsfile[ext].header["CRVAL1"]
            dec = fitsfile[ext].header["CRVAL2"]
        elif "CRVAL1" in fitsfile[1].header:
            ra = fitsfile[1].header["CRVAL1"]
            dec = fitsfile[1].header["CRVAL2"]
        else:
            raise ValueError("No RA and DEC found in the header of %s" % file)
        print(rf"Central coordinates (RA, DEC): %.4f, %.4f" % (ra, dec))
        output = run_get_ebv_script(ra, dec)
        print(output)
        EBV_gal, EBV_gal_err = parse_ebv_output(output)
    else:
        EBV_gal, EBV_gal_err = EBV_arg
        print(fr"Using EBV: %.3f$\pm$%.3f" % (EBV_gal, EBV_gal_err))
    return EBV_gal, EBV_gal_err

Rv_CARDELLI = 3.1

lines_path = "." 
# lines_path = "./camel_*/cleaned_images/"
#

ap = argparse.ArgumentParser(description='Applies extinction correction to a set of maps given an input HALPH and HBETA map. It uses the Cardelli extinction curve + the Calzetti extinction curve for extragalactic extinction')
ap.add_argument("-r", "--rootname", nargs='?', help="Root name of the input line files", type=str, default="clean_camel")
ap.add_argument("-R", "--ratio", help="Ratio of total to selective extinction Av/E(B-V). Default Rv=4.05 from Calzetti", type=float, default=4.05, nargs="?")
ap.add_argument("-i", "--intrinsic", help="Intrinsic Balmer decrement ratio. Default 2.86", type=float, default=2.86, nargs="?")
ap.add_argument("--ebv_gal", help="Line of sight (galactic) EBV to correct for foreground extinction. Uses Rv=3.1 and Cardelli extinction curve. Default gets the line of sight extinction from the cube centre coordinates",
               type=float, nargs=2)
ap.add_argument("-b", "--hbeta", help="Path to HBETA flux map directory", type=str, required=True, nargs=1)
ap.add_argument("-a", "--halpha", help="Path to HALPHA flux map directory", type=str, required=True, nargs=1)
args = ap.parse_args()
halpha = 6563.0 * u.angstrom
hbeta = 4861.0 * u.angstrom

outdir = "deredden_momcheva"
#'OII3727', 'OII3729',
lines = ['HBETA', 'OI6300', 'NII6548', 'NII6583', 'SII6716', 'SII6731', 'OIII4959', 'OIII5007', 'HALPHA', "HeII4686", "HeI5875", "HeI6678", "HeI7065", "ArIII7135", "SIII9069"]
fluxkeyword = "amplitude" # could be "flux"
wavelengthkeyword = "center" # wave

if not os.path.isdir(outdir):
    os.mkdir(outdir)
# Rv = E(B - V)/Av
Rv = args.ratio
intrinsic_ratio = args.intrinsic
EBV_gal = args.ebv_gal

curve = "calzetti"
extincton_curve = C00(Rv)
print('%s/%s_*[!e]%s_%s.fits' % (args.halpha[0], args.rootname, "HALPHA", fluxkeyword))
try:
    halpha_file = glob.glob(f'{args.halpha[0]}/{args.rootname}_HALPHA*[!e]{fluxkeyword}.fits')[0]
except IndexError:
    raise IndexError("Halpha %s file not found in %s!" % (fluxkeyword, args.halpha[0]))
print(f"Halpha file: {halpha_file}")
hbeta_file = glob.glob(f'{args.hbeta[0]}/{args.rootname}_HBETA*[!e]{fluxkeyword}.fits')[0]
print(f"Hbeta file: {hbeta_file}")
halpha_map = fits.open(halpha_file, extname="IMAGE")
ext = 0
if "WCSAXES" in  halpha_map[ext].header:
    ext = 0
        
elif "WCSAXES" in halpha_map[1].header:
    ext = 1
if halpha_map[ext].header["WCSAXES"] == 3:
    halpha_map[ext].header["WCSAXES"] = 2
print("Extension used for Halpha map: %d" % ext)
hbeta_map = fits.open(hbeta_file, extname="IMAGE")

# correct for galactic extinction
EBV_gal = parse_ebv_arg(args.ebv_gal, halpha_file) # [0] value and [1] error

if EBV_gal[0] > 0:
    print(rf"Applying Galactic extinction EBV = %.3f$\pm$%.5f" % (EBV_gal[0], EBV_gal[1]))
    Av_gal = Rv_CARDELLI * EBV_gal[0]
    zeros = np.zeros(halpha_map[ext].data.shape)
    ones = np.ones(halpha_map[ext].data.shape)
    # assume no errors for now
    halpha_map[ext].data, _ = galactic_extinction(halpha.value * ones, halpha_map[ext].data,
                                                 zeros, EBV_gal[0], 0, Rv_CARDELLI)
    hbeta_map[ext].data, _ = galactic_extinction(hbeta.value * ones, hbeta_map[ext].data,
                                                 zeros, EBV_gal[0], 0, Rv_CARDELLI)

observed_ratios = halpha_map[ext].data / hbeta_map[ext].data
# get E(beta-halpha)
color_excess = 2.5 * np.log10(observed_ratios / intrinsic_ratio)
color_excess_map = fits.PrimaryHDU(data=color_excess, header=halpha_map[ext].header)
color_excess_map.writeto("%s/decrement_color_excess_%.3f.fits" % (outdir, args.intrinsic),
                         overwrite=True)

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
# we assume if EBV is negative that there's no extinction
Ebv[Ebv < 0] = 0
color_excess_map = fits.PrimaryHDU(data=Ebv, header=halpha_map[0].header)
out_color_excess = "color_excess_%s_Rv%.1f_i%.2f.fits" % (curve, Rv, args.intrinsic)
color_excess_map.writeto(f"{outdir}/{out_color_excess}", overwrite=True)
print("Color excess map saved to %s/%s" % (outdir, out_color_excess))
print("Mean color excess: %.2f+-%.2f" % (np.nanmean(Ebv), np.nanstd(Ebv)))

# compute uncertainties
halpha_efile = glob.glob(f'{args.halpha[0]}/{args.rootname}_HALPHA_e{fluxkeyword}.fits')[0]
hbeta_efile = glob.glob(f'{args.hbeta[0]}/{args.rootname}_HBETA_e{fluxkeyword}.fits')[0]
if len(halpha_efile) > 0:
    uncertainties = 1
    print("HALPHA uncertainty file found: %s" % halpha_efile)
    print("HBETA uncertainty file found: %s" % hbeta_efile)
    halpha_emap = fits.open(halpha_efile, extname="IMAGE")
    hbeta_emap = fits.open(hbeta_efile, extname="IMAGE")

    halpha_map[ext].data, halpha_emap[ext].data = galactic_extinction(halpha.value * ones,
                                                                      halpha_map[ext].data, halpha_emap[ext].data,
                                                                      EBV_gal[0], EBV_gal[1], Rv_CARDELLI)
    hbeta_map[ext].data, hbeta_emap[ext].data = galactic_extinction(hbeta.value * ones,
                                                                      hbeta_map[ext].data, hbeta_emap[ext].data,
                                                                      EBV_gal[0], EBV_gal[1], Rv_CARDELLI)

    err_ebv = 2.5 * np.sqrt((np.log10(np.e) * (halpha_emap[ext].data / halpha_map[ext].data)) ** 2 + (np.log10(np.e) * hbeta_emap[ext].data / hbeta_map[ext].data)**2) / curve_color_excess
    # I thought the above was incorrect, but comparing two different calculations I get exactly same results
    ecolor_excess_map = fits.PrimaryHDU(data=err_ebv, header=halpha_map[ext].header)
    ecolor_excess_map.writeto("%s/ecolor_excess_%s_Rv%.1f_i%.2f.fits" % (outdir, curve, Rv, args.intrinsic), overwrite=True)
else:
    uncertainties = 0
    print("Warning; halpha error file not found. Uncertainties won't be computed")

for line in lines:
    print("\tDereddening %s line" % line)
    linepath  = f'{lines_path}/{args.rootname}_{line}*[!e]{wavelengthkeyword}.fits'
    wavemap = glob.glob(linepath)
    if len(wavemap) == 0:
        print(f"Line map for {line} line not found ({linepath})")
        continue
    wavemap = wavemap[0]

    fluxmap = glob.glob(f'{lines_path}/{args.rootname}_{line}*[!e]{fluxkeyword}.fits')[0]
    efluxmap = glob.glob(f'{lines_path}/{args.rootname}_{line}*e{fluxkeyword}.fits')[0]
    print("Using wavelength map: %s" % wavemap)
    print("Using flux map: %s" % fluxmap)
    print("Using eflux map: %s" % efluxmap)
    wavelenghts_fits = fits.open(wavemap, extname="IMAGE")
    ext = 0
    if wavelenghts_fits[ext].data is None:
        ext = 1
    wavelenghts = wavelenghts_fits[ext].data
    fluxes_fits = fits.open(fluxmap, extname="IMAGE")[ext]
    efluxes_fits = fits.open(efluxmap, extname="IMAGE")[ext]
    fluxes = fluxes_fits.data
    efluxes = efluxes_fits.data
    # deredden fluxes
    fluxes, efluxes = galactic_extinction(wavelenghts, fluxes, efluxes, EBV_gal[0],
                                          EBV_gal[1], Rv=Rv_CARDELLI)
    for x in range(wavelenghts.shape[0]):  # loop on the x range
        for y in range(wavelenghts.shape[1]):  # loop on the y range
            if np.isnan(Ebv[x, y]):
                # if EBV is nan we assume no extinction
                continue
                #fluxes[x, y] = np.nan
                #efluxes[x, y] = np.nan
            else:
                extinction_curve_value = extincton_curve.evaluate(wavelenghts[x, y] * u.angstrom)
                fluxes[x, y] = fluxes[x, y] * 10 ** (0.4 * extinction_curve_value * Ebv[x, y])
                if uncertainties:
                    efluxes[x, y] = np.sqrt((efluxes[x, y] * 10 ** (0.4 * extinction_curve_value * Ebv[x, y])) ** 2 \
                    + (err_ebv[x, y] * fluxes[x, y] * 10 ** (0.4 * extinction_curve_value * Ebv[x, y]) * np.log(10**(0.4 * extinction_curve_value)))**2)
    dereddened_fits = fits.HDUList(fits.PrimaryHDU(header=fluxes_fits.header))
    dereddened_fits.append(fits.PrimaryHDU(data=fluxes, header=fluxes_fits.header))
    dereddened_fits[0].header["EXTNAME"] = 'DATA'
    dereddened_fits[0].header["ERRDATA"] = 'STAT'
    dereddened_fits[0].header['CURVE'] = "%s" % curve
    dereddened_fits[0].header['R_v'] = "%.2f" % args.ratio
    dereddened_fits[0].header['R_v_gal'] = "%.2f" % Rv_CARDELLI
    dereddened_fits[0].header['EBV_gal'] = "%.3f$\pm$%.4f" % (EBV_gal[0], EBV_gal[1])
    dereddened_fits[0].header['COMMENT'] = "Balmer ratio files: %s/%s" % (halpha_file, hbeta_file)
    dereddened_fits[0].header['COMMENT'] = "The STAT contains the variance (if uncertainties were calculated)"
    # we save the variance
    dereddened_fits.append(fits.ImageHDU(data=efluxes**2, header=efluxes_fits.header,
                            name="STAT"))

    outfile = os.path.basename(fluxmap).replace(".fits", "_deredden.fits")
    dereddened_fits.writeto(f"{outdir}/{outfile}", overwrite=True)
    #edereddened_fits = fits.PrimaryHDU(data=efluxes, header=efluxes_fits.header)
    #edereddened_fits.header['CURVE'] = "%s" % curve
    #edereddened_fits.header['R_v'] = "%.2f" % args.ratio
    #edereddened_fits.header['Ratio'] = "%.2f" % args.intrinsic
    #edereddened_fits.header['COMMENT'] = "Balmer ratio: %s/%s" % (halpha_file, hbeta_file)
    #if uncertainties:
    #    eoutfile = efluxmap.replace(".fits", "_deredden.fits")
    #    edereddened_fits.writeto(eoutfile, overwrite=True)
    python_argument = "%s -R %.1f -i %.3f -a %s -b %s" % (__file__, Rv, args.intrinsic, args.halpha[0], args.hbeta[0])
    with open("%s/python_command.txt" % outdir, "w+") as f:
        f.write(python_argument)
    print("Dereddened flux map for line %s and stored it to %s" % (line, outfile))
