# @Author: Andrés Gúrpide <agurpide>
# @Date:   20-05-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 10-02-2025
from astropy.io import fits
import argparse
import numpy as np
from deredden_utils import galactic_extinction, run_get_ebv_script, parse_ebv_output
from mpdaf.obj import Cube



ap = argparse.ArgumentParser(description='Script to deredden a cube from extinction along the line of sight. Errors are propagated too.')
ap.add_argument("input_cube", help="Path to input cube")
ap.add_argument("-R", "--ratio", help="Ratio of total to selective extinction Av/E(B-V). Default Rv=3.1 from Cardelli",
                type=float, default=3.1, nargs="?")
ap.add_argument("--ebv_gal", help="Line of sight (galactic) EBV to correct for foreground extinction. Uses Cardelli extinction curve." + \
                                "By default EBV (mean and std) is obtained from Schlafly & Finkbeiner 2011 from the cube central coordinates",
               type=float, required=False, nargs=2)
args = ap.parse_args()

curve = "Cardelli"
# Rv = E(B - V)/Av
Rv = args.ratio
inputcube = args.input_cube

if args.ebv_gal is None:

    fitsfile = fits.open(inputcube)
    ra = fitsfile[1].header["CRVAL1"]
    dec = fitsfile[2].header["CRVAL2"]
    print("Cube central coordinates (RA, DEC): %.4f, %.4f" % (ra, dec))
    output = run_get_ebv_script(ra, dec)
    print(output)
    EBV_gal, EBV_gal_err = parse_ebv_output(output)
    fitsfile = None
else:
    EBV_gal, EBV_gal_err = args.ebv_gal
    print("Using EBV: %.3f$\pm$%.3f" % (EBV_gal, EBV_gal_err))

cube = Cube(inputcube)
wavelenghts = cube.wave.coord()
wavelengths_3d = cube.wave.coord()[:, np.newaxis, np.newaxis]
wavelengths_broadcasted = np.broadcast_to(wavelengths_3d, cube.data.shape)
print("Dereddening... this may take a while...")
fluxes, efluxes = galactic_extinction(wavelengths_broadcasted, cube.data,
                                          cube.var**0.5, EBV_gal, EBV_gal_err, Rv=Rv)

cube.data = fluxes
cube.var = efluxes**2.
cube.primary_header['CURVE'] = "%s" % curve
cube.primary_header['R_v'] = "%.2f" % args.ratio
cube.primary_header['EBV'] = "%.2f$\pm$%.4f" % (EBV_gal, EBV_gal_err)
cube.primary_header['COMMENT'] = "The STAT contains the variance"
outfile = inputcube.replace(".fits", "_deredden.fits")
cube.write(outfile)
print("Dereddened cube map stored to %s" % (outfile))
