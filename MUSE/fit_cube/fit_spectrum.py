#!/usr/bin/env python
# -*- coding: utf-8 -*
import numpy as np
import argparse, os, logging
from line import Lines
from lineutils import get_instrument_FWHM, peculiar_velocity, sigma_to_fwhm, ckms
from fitutils import fit_spectrum, DofError, SpectrumMaskedError, plot_fit
from mpdaf.obj import Spectrum
import numpy.ma as ma
from lmfit.model import save_model
from lmfit.printfuncs import ci_report
from math import pi


logger = logging.getLogger('fit_spectrum')

def setup_logging(logfile, loglevel=logging.DEBUG):
    handlers = [logging.StreamHandler()]
    if loglevel==logging.DEBUG:
        handlers.append(logging.FileHandler(logfile))
        logger.debug(f"Started logging to {logfile}")
    """Set up logging configuration"""
    logging.basicConfig(
        level=loglevel,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )

logging.getLogger('matplotlib').setLevel(logging.WARNING)
CATALOG_LINES = Lines().lines
sqrt2pi = (2. * np.pi) ** 0.5
CM_TO_PC = 3.24077929e-19

def check_input_lines(linegroups):
    """Check if the input lines are valid and exist in the catalog, and are in ascending wavelength order
    
    """
    for group in linegroups:
        lines_in_group = group.split(",")
        
        # Check if lines exist in catalog
        for line in lines_in_group:
            if line not in CATALOG_LINES:
                raise ValueError(f"Line {line} not found in catalog. Please check the line name or add it to the catalog. Available lines: {CATALOG_LINES.keys()}")
        
        # Check if lines are in ascending wavelength order
        if len(lines_in_group) > 1:
            wavelengths = [CATALOG_LINES[line].wave for line in lines_in_group]
            if wavelengths != sorted(wavelengths):
                raise ValueError(f"Lines in group '{group}' are not in ascending wavelength order. "
                               f"Current order: {', '.join([f'{line}({CATALOG_LINES[line].wave:.1f}Å)' for line in lines_in_group])}. "
                               f"Please reorder them by wavelength.")
    
    return True


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Fit data cube using voronoi binning')
    parser.add_argument("-d", "--degree", nargs='?', help="Polynomial degree for the continuum", default=1, type=int)
    parser.add_argument("-z", "--redshift", nargs='?', help="Initial guess for the redshift", default=0.0, type=float)
    parser.add_argument("-D", "--distance", nargs='?', help="Optional distance in Mpc to estimate line luminosities", type=float)
    parser.add_argument("-s", "--sigma", nargs='?', help="Initial guess for the width (sigma) of the lines in Angstroms. Default 1.4 Angstroms", default=1.4, type=float)
    parser.add_argument("-o", "--outdir", nargs='?', help="Name of the output directory", default="fit_spectrum", type=str)
    parser.add_argument("input_spectrum", help="Path to input spectrum to be analysed", type=str)
    parser.add_argument("-w", "--window", help="+- Angstroms to cut the spectrum around each line centroid. Default 25 Angstroms. Reduce for lines close to telluric lines (e.g. [OI]6300)", type=float, default=25)
    parser.add_argument("-l", "--linegroups", nargs='+', default="HeII4686", help="Space-separated groups of comma-separated lines to fit. (e.g. if you want two groups with same continuum: HBETA OIII4959,OIII5007)", type=str) # 
    args = parser.parse_args()

    outpath = args.outdir
    sigma = args.sigma
    degree = args.degree
    if args.distance is not None:
        D = args.distance * 1e6  / CM_TO_PC

    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    check_input_lines(args.linegroups)

    # Create log file name based on line groups
    linegroupout = "_".join([line for group in args.linegroups for line in group.split(",")])
    logfile = f"{outpath}/fit_spectrum_{linegroupout}.log"
    setup_logging(logfile, logging.DEBUG)
    
    inputgroups = args.linegroups#[["HeII4686"]]#[["HBETA", "OIII4959", "OIII5007"], ["OI6300"],
            #["NII6548","HALPHA", "NII6583"], ["SII6716", "SII6731"]]
            # ["NII5755"], ["HeI5875.64"], 
    logger.info("Line groups: %s",inputgroups)

    redshift = args.redshift
    # Margins to cut the spectrum around the blue and red lines of each group
    margin = args.window # Angstroms

    # create waveleneght cuts
    wav_cuts = np.zeros((len(inputgroups), 2), dtype=float)
    # we need these variables for later
    fit_lines = {}
    linegroups = []

    for index, group in enumerate(inputgroups):
        grouplines = group.split(",")
        linegroups.append(grouplines)
        # make wavelength range based on the first and last line of the group
        
        minwav = CATALOG_LINES[grouplines[0]].wave * (1 + redshift) - margin
        maxwav = CATALOG_LINES[grouplines[-1]].wave * (1 + redshift) + margin
        wav_cuts[index] = (minwav, maxwav)
        for line in grouplines:
            if line in fit_lines.keys():
                logger.info(f"Adding broad component for {line}")
                fit_lines[f"{line}_broad"] =  CATALOG_LINES[line]
            else:
                fit_lines[line] = CATALOG_LINES[line]


    # get the first line of the first group and the last line of the last group
    waveinit, waveend = list(fit_lines.values())[0].wave * (1 + redshift) - margin, list(fit_lines.values())[-1].wave * (1 + redshift) + margin
    specname = args.input_spectrum
    data_spectrum = Spectrum(specname)
    step = data_spectrum.wave.get_step()
    # cut the cube over the wavelength range we will not needed
    data_spectrum.mask_region(waveinit - step *1.01, waveend + step * 1.01, inside=False)
    data_spectrum.crop()

    # mask the data spectrum outside the wavelength range of interest, important to copy here
    original_mask = data_spectrum.data.mask.copy()  
    # mask the entire range by default
    data_spectrum[:] = ma.masked
    # in the wavelength range of interest, keep the original mask
    for index, wav_range in enumerate(wav_cuts):
        lmin, lmax = wav_range
        logger.info(f"Masking wavelength range {lmin} - {lmax} Angstroms")
        minindex = data_spectrum.wave.pixel(lmin, nearest=True)
        maxindex = data_spectrum.wave.pixel(lmax, nearest=True)
        # Unmask this wavelength range, but preserve original mask
        data_spectrum.data.mask[minindex:maxindex] = original_mask[minindex:maxindex]

    wavelengths = data_spectrum.wave.coord()

    logger.info(f"Fitting spectrum {specname} with {len(fit_lines)} lines,  {len(linegroups)} line groups")
    
    try:
        result, conf = fit_spectrum(data_spectrum, fit_lines, redshift=redshift, sigma=sigma, wavelengths=wavelengths, degree=degree, uncertainties=True)
    except (SpectrumMaskedError, DofError) as e:
        logger.error(f"Fitting failed: {e}")
        exit(1)
    save_model(result, f"{outpath}/{linegroupout}_model.sav")
    print(result.fit_report())
    print(ci_report(conf, ndigits=3))
    fig = plot_fit(wavelengths[~data_spectrum.mask], result, lref=fit_lines[line].wave, annotate=False, z_sys=redshift, normalize=False)

    fig.savefig(f"{outpath}/{linegroupout}_fit")

    cont_prefix = 'cont'
    noise = np.std(result.data - result.best_fit)
    redchi = result.redchi
    bic = result.bic
    rsquared = result.rsquared
    ndata = data_spectrum.data.size
    nvarys = result.nvarys
    # Collect all parameter values and errors
    param_values = []
    param_names = []
    
    # Add basic fit statistics
    param_names.extend(['redchi',  'bic', 'rsquared', 'ndof', 'ndata', 'ier'])
    param_values.extend([redchi, bic, rsquared, 
                        ndata - nvarys, ndata, result.ier])
    # Process each line
    previouslinename = None
    for linename in fit_lines.keys():

        if previouslinename ==linename:
            linename += "_broad"
        line = fit_lines[linename]
        best_values = result.params
        
        # Calculate SNR
        signal = best_values["%s_height" % linename].value
        snr = signal / noise if noise > 0 else np.nan
        
        param_names.append(f"{linename}_snr")
        param_values.append(snr)
        
        # Store basic line parameters
        linepars = ["amplitude", "fwhm", "center"]
        for par in linepars:
            param_name = "%s_%s" % (linename, par)
            if param_name in best_values:
                param = best_values[param_name]
                
                # Parameter value
                param_names.append(param_name)
                param_values.append(param.value)
                
                # Parameter uncertainty
                error_value = np.nan
                if conf is not None:
                    # Handle FWHM specially (stored as sigma)
                    if par == "fwhm":
                        parname = "%s_sigma" % linename
                        if parname in conf:
                            confpar = conf[parname]
                            error_value = (confpar[1][1] - confpar[0][1]) * sigma_to_fwhm
                        # if it was a broad line
                        elif "broad" in linename:
                            confpar = conf["%s_sigma_factor" %linename]
                            # multiply sigma of the narrow line with the factor
                            print(confpar[0][1], confpar[1][1])
                            error_value = (confpar[1][1] - confpar[0][1]) * best_values["%s_fwhm" % linename.replace("_broad", "")]
                        
                        # store for later
                        fwhm_value = param.value
                        efwhm_value = error_value
                    elif param_name in conf:
                        confpar = conf[param_name]
                        error_value = confpar[1][1] - confpar[0][1]

                elif param.stderr is not None:
                    error_value = param.stderr
                if par=="center":
                    # store for later
                    wavelength = param.value
                    ewavelength = error_value
                param_names.append(f"{param_name}_err")
                param_values.append(error_value)

                if par=="amplitude" and args.distance is not None:
                    param_names.append("%s_lum" % linename)
                    param_values.append(param.value * 4 * pi * D**2)
                    param_names.append("%s_lum_err" % linename)
                    param_values.append(error_value * 4 * pi * D**2)
        
        # FWHM corrected and FWHM velocity
        FWHM_inst, eFWHM_inst = get_instrument_FWHM(wavelength, ewavelength)
        if fwhm_value > FWHM_inst:
            fwhm_corrected = np.sqrt(fwhm_value**2 - FWHM_inst**2)
            efwhm_corrected = np.sqrt((fwhm_value * efwhm_value)**2 + (FWHM_inst * eFWHM_inst)**2) / fwhm_corrected
            fwhm_vel = fwhm_corrected / wavelength * ckms
            efwhm_vel = np.sqrt(efwhm_corrected**2 + (fwhm_corrected * ewavelength / wavelength)**2) * ckms / wavelength
        # unresolved line
        else:
            fwhm_corrected = 0
            efwhm_corrected = np.nan
            fwhm_vel = 0
            efwhm_vel = np.nan


        param_names.extend([f"{linename}_fwhm_corrected", f"{linename}_fwhm_corrected_err", f"{linename}_fwhm_vel", f"{linename}_fwhm_vel_err"])
        param_values.extend([fwhm_corrected, efwhm_corrected, fwhm_vel, efwhm_vel])

        # Velocity
        if "sky" in linename:
            z_sys = 0.0
        else:
            z_sys = redshift
        velocity, evelocity = peculiar_velocity(line.wave, wavelength, 
                                              ewavelength=ewavelength, z_sys=z_sys)
        
        param_names.extend([f"{linename}_vel", f"{linename}_vel_err"])
        param_values.extend([velocity, evelocity])

        previouslinename = linename
    
    # Add continuum parameters
    for i in range(degree + 1):
        param = best_values[f"{cont_prefix}_c{i}"]
        param_names.append(f"cont_c{i}")
        param_values.append(param.value)
        
        # Continuum parameter uncertainty
        error_value = np.nan
        if conf is not None and f"{cont_prefix}_c{i}" in conf:
            confpar = conf[f"{cont_prefix}_c{i}"]
            error_value = np.abs(confpar[0][1] - confpar[1][1])
        elif param.stderr is not None:
            error_value = param.stderr
        
        param_names.append(f"cont_c{i}_err")
        param_values.append(error_value)
    
    # Write results to text file
    header = "\t".join(param_names)
    np.savetxt(f"{outpath}/{linegroupout}_fit_results.dat",
               np.array(param_values).reshape(1, -1),
               header=header,
               fmt="%.6g",
               delimiter="\t")
    logger.info("Output stored in %s" % outpath)