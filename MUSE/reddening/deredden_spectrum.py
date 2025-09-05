    # @Author: Andrés Gúrpide <agurpide>
# @Date:   20-05-2025
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 10-02-2025
from astropy.io import fits
import argparse
import os
import astropy.units as u
from fitting.line import Lines
import sys
sys.path.append("/home/andresgur/scripts/pythonscripts/maxime_muse/fitting")
from fitutils import fit_spectrum, plot_fit
from deredden_utils import galactic_extinction, C00
from mpdaf.obj import Spectrum
from extinction import calzetti00, remove
from lmfit.model import save_model
from lmfit.printfuncs import ci_report
from math import log, log10

def division(a, b, a_err, b_err):
    """Compute the error on the division of two values or arrays (a/b)
    ----------
    Parameters
    a: value or np.array  of the numerator
    b: value or np.array of the denominator
    a_err: error on the numerator (or array)
    b_err: error on the denominator (or array)
    """

    return a/b, ((a_err / b) ** 2 + (a * b_err / b ** 2) ** 2)**0.5


def error_log10(value, error):
    """Compute the error on the log10 of a value
    
    """
    return error / (value * log(10))

CATALOG_LINES = Lines().lines


Rv_CARDELLI = 3.1

lines_path = "." 
# lines_path = "./camel_*/cleaned_images/"
#
if __name__ == "__main__":
    ap = argparse.ArgumentParser(description='Applies extinction correction to a spectrum based on the derived fluxes of Ha and Hb')
    ap.add_argument("input_spectrum", help="Path to input spectrum to be analysed", type=str)
    ap.add_argument("-Rv", "--Rv", help="Ratio of total to selective extinction Av/E(B-V). Default Rv=4.05 from Calzetti. Set to 0 to correct for galactic extinction only", 
                    type=float, default=4.05, nargs="?")
    ap.add_argument("-s", "--sigma", nargs='?', help="Initial guess for the width (sigma) of the lines in Angstroms. Default 1.4 Angstroms", default=1.4, type=float)
    ap.add_argument("-z", "--redshift", nargs='?', help="Initial guess for the redshift", default=0.0, type=float)
    ap.add_argument("-i", "--intrinsic", help="Intrinsic Balmer decrement ratio. Default 2.86", type=float, default=2.86, nargs="?")
    ap.add_argument("--EBV_gal", help="Line of sight (galactic) EBV to correct for foreground extinction. Uses Rv=3.1 and Cardelli extinction curve.",
                type=float, nargs=2, required=True)
    args = ap.parse_args()
    halpha = CATALOG_LINES["HALPHA"].wave
    hbeta = CATALOG_LINES["HBETA"].wave

    outdir = "deredden_spectrum"

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    # Rv = E(B - V)/Av
    Rv = args.Rv
    intrinsic_ratio = args.intrinsic
    EBV_gal = args.EBV_gal
    redshift = args.redshift
    sigma = args.sigma
    intrinsic_ratio = args.intrinsic

    curve = "calzetti"
    extincton_curve = C00(Rv)

    input_spectrum = args.input_spectrum

    spectrum = Spectrum(input_spectrum)
    datafluxes = spectrum.data
    wavelengths = spectrum.wave.coord()
    std = spectrum.var**0.5

    spectrumfits = fits.open(input_spectrum)

    dereddened_fits = fits.HDUList(fits.PrimaryHDU(header=spectrumfits[0].header))

    # correct for galactic extinction

    if EBV_gal[0] > 0:
        print(rf"Applying Galactic extinction EBV = %.3f$\pm$%.5f" % (EBV_gal[0], EBV_gal[1]))
        EBV_gal, EBV_gal_err = EBV_gal
        Av_gal = Rv_CARDELLI * EBV_gal
        galcorrectedfluxes, egalcorrectedfluxes = galactic_extinction(wavelengths, datafluxes, std, EBV_gal, EBV_gal_err, Rv_CARDELLI)
        dereddened_fits[0].header['R_v_gal'] = r"%.2f" % Rv_CARDELLI
        dereddened_fits[0].header['EBV_gal'] = r"%.3f" % EBV_gal
        dereddened_fits[0].header['EBV_galerr'] = r"%.4f" % EBV_gal_err

        correctedspectrum = spectrum.copy()
        correctedspectrum.data = galcorrectedfluxes
        correctedspectrum.var = egalcorrectedfluxes**2.
    else:
        correctedspectrum = spectrum

    if args.Rv==0:
        print("Correcting for galactic extinction only, not applying any other extinction correction")
        exit(1)

    margin = 10
    step = spectrum.wave.get_step()

    line_fluxes = {}

    for line in "HALPHA", "HBETA":
        print(f"Fitting {line}")
        linegroups = [[line]]
        fit_lines = {}
        minwav = CATALOG_LINES[line].wave * (1 + redshift) - margin
        maxwav = CATALOG_LINES[line].wave * (1 + redshift) + margin
        fit_lines[line] = CATALOG_LINES[line]

        data_spectrum = correctedspectrum.copy()
        
        # cut the cube over the wavelength range we will not needed
        data_spectrum.mask_region(minwav - step *1.01, maxwav + step * 1.01, inside=False)
        data_spectrum.crop()
        linewavelengths = data_spectrum.wave.coord()
        result, conf = fit_spectrum(data_spectrum, fit_lines, redshift=redshift, sigma=sigma, wavelengths=linewavelengths, degree=1, uncertainties=True)
        save_model(result, f"{outdir}/{line}_model.sav")
        print(ci_report(conf, ndigits=3))
        fig = plot_fit(linewavelengths[~data_spectrum.mask], result)

        fig.savefig(f"{outdir}/{line}_fit.png", dpi=200)

        flux = result.best_values[line + "_amplitude"]
        eflux = abs(flux - conf["%s_amplitude" % line][0][1])
        line_fluxes["%s" % line] = (flux, eflux)

    observed_ratio, eobserved_ratio = division(line_fluxes["HALPHA"][0], line_fluxes["HBETA"][0], line_fluxes["HALPHA"][1], line_fluxes["HBETA"][1])
    ratio = observed_ratio / intrinsic_ratio
    color_excess, ecolor_excess = 2.5 * log10(ratio), 2.5 * error_log10(ratio, eobserved_ratio / intrinsic_ratio)
    extinction_curve = C00(Rv)
    kbeta = extinction_curve.evaluate(CATALOG_LINES["HBETA"].wave * u.AA)
    kalpha = extinction_curve.evaluate(CATALOG_LINES["HALPHA"].wave * u.AA)
    curve_color_excess = (kbeta - kalpha)
    EBV, EBV_err = division(color_excess, curve_color_excess, ecolor_excess, 0)
    print(r"Observed ratio: %.2f/%.2f (%.2f\pm%.2f)" % (line_fluxes["HALPHA"][0], line_fluxes["HBETA"][0], observed_ratio, eobserved_ratio))

    if EBV - EBV_err < 0:
        print("EBV is negative, avoiding extra extinction correction")
        Av = 0    
        Av_err = 0

    else:
        Av = Rv * EBV
        Av_err = Rv * EBV_err
        print(r"Derived color excess E(B-V) = %.3f$\pm$%.3f" % (EBV, EBV_err))
        print(r"Applying extinction correction with Rv=%.2f, Av=%.3f$\pm$%.3f" % (Rv, Av, Av_err))

    wavelengths = correctedspectrum.wave.coord()
    mag_ext = calzetti00(wavelengths, Av, Rv, unit="aa") # wavelength in angstroms
    fluxes = correctedspectrum.data
    dereddened_fluxes = remove(mag_ext, correctedspectrum.data)
    vardereddened_fluxes = remove(mag_ext, correctedspectrum.var**0.5)**2 + (Av_err * fluxes * log(10) * 0.4)**2
    dereddened_fits.append(fits.ImageHDU(data=dereddened_fluxes.data, header=spectrumfits[1].header, name="DATA"))
    dereddened_fits[0].header['CURVE'] = r"%s" % curve
    dereddened_fits[0].header['R_v'] = "%.2f" % args.Rv
    dereddened_fits[0].header['EBV'] = "%.3f" % EBV
    dereddened_fits[0].header['EBV_err'] = "%.3f" % EBV_err
    dereddened_fits[0].header['COMMENT'] = "Observed ratio: %.2f/%.2f (%.2f)" % (line_fluxes["HALPHA"][0], line_fluxes["HBETA"][0], observed_ratio)
    dereddened_fits[0].header['COMMENT'] = "The STAT contains the variance (if uncertainties were calculated)"
    # we save the variance
    dereddened_fits.append(fits.ImageHDU(data=vardereddened_fluxes.data, header=spectrumfits[2].header,
                            name="STAT"))
    spectrumfits.close()

    outfile = os.path.basename(input_spectrum).replace(".fits", "_deredden.fits")
    dereddened_fits.writeto(f"{outdir}/{outfile}", overwrite=True)
    print("Dereddened spectrum stored to %s" % (outfile))
