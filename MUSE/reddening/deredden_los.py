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
from extinction import ccm89, remove
import subprocess
import errno
import re

PATH = "/home/andresgur/scripts/utils"

def run_get_ebv_script(ra, dec):
    # Path to the get_ebv.sh script
    script_path = "%s/get_ebv.sh" % PATH

    # Run the script with the provided RA and Dec
    result = subprocess.run([script_path, "%.5f" % ra, "%.5f" % dec], capture_output=True, text=True)

    # Check if the script ran successfully
    if result.returncode != 0:
        print("Error running the get_ebv.sh script:")
        print(result.stderr)
        return None

    # Capture the output
    output = result.stdout.strip()
    return output

def parse_ebv_output(output):
    """Parses the output of get_ebv.sh"""
    # Regular expressions to capture the mean and std values
    mean_pattern = re.compile(r"E\(B-V\) mean \(Schlafly & Finkbeiner 2011\) over \d+ arcminutes: ([0-9.]+) mag")
    std_pattern = re.compile(r"E\(B-V\) std \(Schlafly & Finkbeiner 2011\) over \d+ arcminutes: ([0-9.]+) mag")

    mean_match = mean_pattern.search(output)
    std_match = std_pattern.search(output)

    if mean_match and std_match:
        mean_value = float(mean_match.group(1))
        std_value = float(std_match.group(1))
        return mean_value, std_value
    else:
        print("Error: Could not parse the mean or std values from the output.")
        return None


def galactic_extinction(wavelengths, fluxes, efluxes, EBV_gal, EBV_gal_err=0, Rv=3.1):
    """Wavelengths in angstroms"""
    Av = EBV_gal * Rv
    Av_gal_err = EBV_gal_err * Rv
    wavs = np.array(wavelengths.flatten(), dtype="double")
    mag_ext = ccm89(wavs, Av, Rv, unit="aa")
    deredden_fluxes = remove(mag_ext, fluxes.flatten())
    ederedden_fluxes = np.sqrt(remove(mag_ext, efluxes.flatten())**2 + (Av_gal_err * fluxes.flatten() * np.log(10) * 0.4)**2)
    return deredden_fluxes.reshape(wavelengths.shape), ederedden_fluxes.reshape(wavelengths.shape)


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
