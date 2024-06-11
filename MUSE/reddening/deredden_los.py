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
import glob
import errno
from deredden_utils import galactic_extinction, run_get_ebv_script, parse_ebv_output


ap = argparse.ArgumentParser(description='Script to deredden the line maps from extinction along the line of sight. One fits file with "DATA" and "STAT" (variance) is created')
ap.add_argument("-r", "--rootname", nargs='?', help="Root name of the input line files", type=str, default="cleaned_images/clean_camel")
ap.add_argument("-R", "--ratio", help="Ratio of total to selective extinction Av/E(B-V). Default Rv=3.1 from Cardelli",
                type=float, default=3.1, nargs="?")
ap.add_argument("--ebv_gal", help="Line of sight (galactic) EBV to correct for foreground extinction. Uses Cardelli extinction curve." + \
                                "By default EBV (mean and std) is obtained from Schlafly & Finkbeiner 2011 from the cube central coordinates",
               type=float, required=False, nargs=2)
args = ap.parse_args()

#'OII3727', 'OII3729',
lines = ['HBETA', 'OI6300', 'NII6548', 'NII6583', 'SII6716', 'SII6731', 'OIII4959', 'OIII5007', 'HALPHA']
keyword = "flux" # could be "flux"
curve = "Cardelli"
# Rv = E(B - V)/Av
Rv = args.ratio

if args.ebv_gal is None:
    line = lines[0]
    print("Reading cube central coordinates from %s for EBV calculation:" % line)
    file = './camel_*/%s_*[!e]wave_*%s.fits' % (args.rootname, line)
    wavemap = glob.glob(file)
    if len(wavemap) == 0:
        print("Line map for %s line not found" % line)
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)
    else:
        try:
            fitsfile = fits.open(wavemap[0], extname="IMAGE")
            ext = 0
        except:
            ext = 1
            fitsfile = fits.open(wavemap, extname=ext)
        ra = fitsfile[0].header["CRVAL1"]
        dec = fitsfile[0].header["CRVAL2"]
        print("Cube central coordinates (RA, DEC): %.4f, %.4f" % (ra, dec))
        output = run_get_ebv_script(ra, dec)
        print(output)
        EBV_gal, EBV_gal_err = parse_ebv_output(output)
else:
    EBV_gal, EBV_gal_err = args.ebv_gal
    print("Using EBV: %.3f$\pm$%.3f" % (EBV_gal, EBV_gal_err))

for line in lines:
    print("\tDereddening %s line" % line)
    wavemap = glob.glob('./camel_*/%s_*[!e]wave_*%s.fits' % (args.rootname, line))
    if len(wavemap) == 0:
        print("Line map for %s line not found" % line)
        continue
    wavemap = wavemap[0]

    fluxmap = glob.glob('./camel_*/%s_*[!e]%s_*%s.fits' % (args.rootname, keyword, line))[0]
    efluxmap = glob.glob('./camel_*/%s_*e%s_*%s.fits' % (args.rootname, keyword, line))[0]
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
    fluxes, efluxes = galactic_extinction(wavelenghts, fluxes, efluxes, EBV_gal, EBV_gal_err, Rv=Rv)
    dereddened_fits = fits.HDUList(fits.PrimaryHDU(header=fluxes_fits.header))
    dereddened_fits[0].header["EXTNAME"] = 'DATA'
    dereddened_fits[0].header["ERRDATA"] = 'STAT'
    dereddened_fits[0].header['CURVE'] = "%s" % curve
    dereddened_fits[0].header['R_v'] = "%.2f" % args.ratio
    dereddened_fits[0].header['EBV'] = "%.2f$\pm$%.4f" % (EBV_gal, EBV_gal_err)
    dereddened_fits[0].header['COMMENT'] = "The STAT contains the variance"
    dereddened_fits.append(fits.ImageHDU(data=fluxes, header=fluxes_fits.header,
                           name="DATA"))
    # we save the variance
    dereddened_fits.append(fits.ImageHDU(data=efluxes**2, header=efluxes_fits.header,
                           name="STAT"))
    outfile = fluxmap.replace(".fits", "_deredden.fits")
    dereddened_fits.writeto(outfile, overwrite=True)

    print("Dereddened flux map for line %s and stored it to %s" % (line, outfile))
