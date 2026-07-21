#!/usr/bin/env python
# -*- coding: utf-8 -*
import numpy as np
import argparse, os, logging
from line import Lines
from lineutils import get_instrument_FWHM, compute_shift, ckms, correct_FWHM
from fitutils import fit_spectrum, DofError, SpectrumMaskedError, plot_fit, get_error
from mpdaf.obj import Spectrum
import numpy.ma as ma
from lmfit.model import save_model
from lmfit.printfuncs import ci_report
from math import pi
import matplotlib.pyplot as plt

logger = logging.getLogger("fit_spectrum")


def setup_logging(logfile, loglevel=logging.DEBUG):
    handlers = [logging.StreamHandler()]
    if loglevel == logging.DEBUG:
        handlers.append(logging.FileHandler(logfile))
        logger.debug(f"Started logging to {logfile}")
    """Set up logging configuration"""
    logging.basicConfig(
        level=loglevel,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=handlers,
    )


logging.getLogger("matplotlib").setLevel(logging.WARNING)
CATALOG_LINES = Lines().lines
CM_TO_PC = 3.24077929e-19


def check_input_lines(linegroups):
    """Check if the input lines are valid and exist in the catalog, and are in ascending wavelength order"""
    for group in linegroups:
        lines_in_group = group.split(",")

        # Check if lines exist in catalog
        for line in lines_in_group:
            if line.strip("_") not in CATALOG_LINES:
                raise ValueError(
                    f"Line {line} not found in catalog. Please check the line name or add it to the catalog. Available lines: {CATALOG_LINES.keys()}"
                )

        # Check if lines are in ascending wavelength order
        if len(lines_in_group) > 1:
            wavelengths = [CATALOG_LINES[line.strip("_")].wave for line in lines_in_group]
            if wavelengths != sorted(wavelengths):
                raise ValueError(
                    f"Lines in group '{group}' are not in ascending wavelength order. "
                    f"Current order: {', '.join([f'{line}({CATALOG_LINES[line].wave:.1f}Å)' for line in lines_in_group])}. "
                    f"Please reorder them by wavelength."
                )

    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Fit data cube using voronoi binning")
    parser.add_argument(
        "-d",
        "--degree",
        nargs="?",
        help="Polynomial degree for the continuum",
        default=1,
        type=int,
    )
    parser.add_argument(
        "-z",
        "--redshift",
        nargs="?",
        help="Initial guess for the redshift",
        default=0.0,
        type=float,
    )
    parser.add_argument(
        "-D",
        "--distance",
        nargs="?",
        help="Optional distance in Mpc to estimate line luminosities",
        type=float,
    )
    
    parser.add_argument(
        "--FWHM",
        nargs="?",
        help="Initial guess for the FWHM width of the lines in numbers of instrumental resolution",
        default=1,
        type=float,
    )
    parser.add_argument(
        "-s",
        "--sigma",
        nargs="?",
        help="Number of sigmas for the uncertainties",
        default=1,
        type=float,
    )
    parser.add_argument("--vmargin", default=300, nargs="?", type=float, help="Margin in velocity space in km/s")
    parser.add_argument(
        "-o",
        "--outdir",
        nargs="?",
        help="Name of the output directory",
        default="fit_spectrum",
        type=str,
    )
    parser.add_argument(
        "input_spectrum", help="Path to input spectrum to be analysed", type=str
    )
    parser.add_argument(
        "-w",
        "--window",
        help="+- Angstroms to cut the spectrum around each line centroid. Default 25 Angstroms. Reduce for lines close to telluric lines (e.g. [OI]6300)",
        type=float,
        default=25,
    )
    parser.add_argument(
        "-l",
        "--linegroups",
        nargs="+",
        default="HeII4686",
        help="Space-separated groups of comma-separated lines to fit. (e.g. if you want two groups with same continuum: HBETA OIII4959,OIII5007)",
        type=str,
    )  #
    parser.add_argument(
        "--FWHM_factor",
        nargs="?",
        default=4,
        help="Maximum FWHM allowed for the narrow line. This is times the instrumental resolution. Default 4",
        type=float,
    )  #
    args = parser.parse_args()

    outpath = args.outdir
    sigma = args.sigma
    degree = args.degree
    if args.distance is not None:
        D = args.distance * 1e6 / CM_TO_PC

    if not os.path.isdir(outpath):
        os.mkdir(outpath)

    check_input_lines(args.linegroups)

    # Create log file name based on line groups
    linegroupout = "_".join(
        [line for group in args.linegroups for line in group.split(",")]
    )
    logfile = f"{outpath}/fit_spectrum_{linegroupout}.log"
    setup_logging(logfile, logging.DEBUG)

    inputgroups = (
        args.linegroups
    )  # [["HeII4686"]]#[["HBETA", "OIII4959", "OIII5007"], ["OI6300"],
    # ["NII6548","HALPHA", "NII6583"], ["SII6716", "SII6731"]]
    # ["NII5755"], ["HeI5875.64"],
    logger.info("Line groups: %s", inputgroups)

    logger.info(f"Degree of continuum polynomial: {degree}")

    redshift = args.redshift
    # Margins to cut the spectrum around the blue and red lines of each group
    margin = args.window  # Angstroms

    FWHM_max_factor = args.FWHM_factor

    # create waveleneght cuts
    wav_cuts = np.zeros((len(inputgroups), 2), dtype=float)
    # we need these variables for later
    fit_lines = {}
    linegroups = []

    for index, group in enumerate(inputgroups):
        grouplines = group.split(",")
        linegroups.append(grouplines)
        # make wavelength range based on the first and last line of the group

        minwav = CATALOG_LINES[grouplines[0].strip("_")].wave * (1 + redshift) - margin
        maxwav = CATALOG_LINES[grouplines[-1].strip("_")].wave * (1 + redshift) + margin
        wav_cuts[index] = (minwav, maxwav)
        for line in grouplines:
            if line.strip("_") in fit_lines.keys():
                logger.info(f"Adding broad component for {line}")
                if line.endswith("_"):
                    fit_lines[f"{line.strip("_")}_broad_"] = CATALOG_LINES[line.strip("_")]
                else:    
                    fit_lines[f"{line}_broad"] = CATALOG_LINES[line]
            if line.endswith("_"):
                logger.info(f"Adding {line.strip("_")} with frozen velocity")
                fit_lines[line] = CATALOG_LINES[line.strip("_")]
            else:
                fit_lines[line] = CATALOG_LINES[line]

    # get the first line of the first group and the last line of the last group
    waveinit, waveend = (
        list(fit_lines.values())[0].wave * (1 + redshift) - margin,
        list(fit_lines.values())[-1].wave * (1 + redshift) + margin,
    )
    specname = args.input_spectrum
    data_spectrum = Spectrum(specname)
    step = data_spectrum.wave.get_step()
    # cut the cube over the wavelength range we will not needed
    data_spectrum.mask_region(
        waveinit - step * 2.01, waveend + step * 2.01, inside=False
    )
    data_spectrum.crop()
    if np.all(data_spectrum.mask):
        raise SpectrumMaskedError(
            "Spectrum is completely masked! Check the input spectrum, try different lines or wavelenght cuts"
        )

    # mask the data spectrum outside the wavelength range of interest, important to copy here

    # here is symmetric
    original_mask = data_spectrum.data.mask.copy()
    # data_spectrum.plot()
    # mask the entire range by default
    data_spectrum[:] = ma.masked
    # in the wavelength range of interest, keep the original mask
    for index, wav_range in enumerate(wav_cuts):
        lmin, lmax = wav_range
        logger.info(f"Masking wavelength range {lmin} - {lmax} Angstroms")
        minindex = data_spectrum.wave.pixel(lmin, nearest=True)
        maxindex = data_spectrum.wave.pixel(lmax + 2, nearest=True)
        # Unmask this wavelength range, but preserve original mask
        data_spectrum.data.mask[minindex:maxindex] = original_mask[minindex:maxindex]

    wavelengths = data_spectrum.wave.coord()
    # data_spectrum.plot()
    
    logger.info(
        f"Fitting spectrum {specname} with {len(fit_lines)} lines,  {len(linegroups)} line groups"
    )
    vmargin=args.vmargin
    try:
        result, conf = fit_spectrum(
            data_spectrum,
            fit_lines,
            redshift=redshift,
            sigma=sigma,
            wavelengths=wavelengths,
            degree=degree,
            uncertainties=True,
            FWHM_max_factor=FWHM_max_factor,
            vmargin=vmargin,
        )
    except (SpectrumMaskedError, DofError) as e:
        logger.error(f"Fitting failed: {e}")
        exit(1)
    save_model(result, f"{outpath}/{linegroupout}_model.sav")
    print(result.fit_report())
    if conf is not None:
        print(ci_report(conf, ndigits=3))
    if len(fit_lines) == 1:
        fig, axes = plot_fit(
            result,
            lref=fit_lines[line].wave,
            z_sys=redshift,
            annotate=False,
            normalize=False,
        )
    else:
        fig, axes = plot_fit(
            result,
            # lref=fit_lines[line].wave,
            z_sys=redshift,
            annotate=False,
            normalize=True,
        )
    #axes[0].text(
    #        0.15, 0.8, "HeII" + r"$\lambda$4686", fontsize=24, transform=axes[0].transAxes
    #)
    #plt.show()
    fig.savefig(f"{outpath}/{linegroupout}_fit.png")

    cont_prefix = "cont"
    best_fit = result.best_fit
    data = result.data
    res = result.data - result.best_fit
    x = result.userkws["x"]
    xnew = np.linspace(x.min(), x.max(), 500)
    comps = result.eval_components(x=xnew)
    # store outputs
    # --- Save data file (observed wavelength grid) ---
    data_out = np.column_stack([
        x,
        data,
        1 / result.weights,
        best_fit,
        res
    ])
    np.savetxt(
        f"{outpath}/{linegroupout}_data.tsv",
        data_out,
        delimiter='\t',
        header='wave\tflux\terror\tbest_fit\tresidual'
    )

    # --- Save model file (fine wavelength grid) ---
    xnew = np.linspace(x.min(), x.max(), len(x) * 500)
    comps = result.eval_components(x=xnew)

    cont_cols   = [comps[k] for k in comps if k.startswith(cont_prefix)]
    line_cols   = [comps[k] for k in comps if not k.startswith(cont_prefix)]
    line_names  = [k.rstrip('_') for k in comps if not k.startswith(cont_prefix)]
    cont_names  = [k.rstrip('_') for k in comps if k.startswith(cont_prefix)]

    continuum = np.sum(cont_cols, axis=0) if cont_cols else np.zeros_like(xnew)
    total     = result.eval(x=xnew)

    model_out = np.column_stack([xnew, total, continuum] + line_cols)
    header    = '\t'.join(['wave', 'total', 'continuum'] + line_names)
    np.savetxt(
        f"{outpath}/{linegroupout}_model.tsv",
        model_out,
        delimiter='\t',
        header=header
    )
    noise = np.std(res)
    redchi = result.redchi
    bic = result.bic
    rsquared = result.rsquared
    ndata = data_spectrum.data.size
    nvarys = result.nvarys
    # Collect all parameter values and errors
    param_values = []
    param_names = []

    # Add basic fit statistics
    param_names.extend(["redchi", "bic", "rsquared", "ndof", "ndata", "ier"])
    param_values.extend([redchi, bic, rsquared, ndata - nvarys, ndata, result.ier])
    # Process each line
    previouslinename = None
    for linename in fit_lines.keys():

        if previouslinename == linename:
            linename += "_broad"
        
        line = fit_lines[linename]
        best_values = result.params

        # Calculate SNR
        signal = best_values["%s_height" % linename].value
        snr = signal / noise if noise > 0 else np.nan

        param_names.append(f"{linename}_snr")
        param_values.append(snr)

        # Store basic line parameters
        linepars = ["amplitude", "vel", "fwhm_vel"]
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
                    if line.ref in fit_lines:
                        # replace only the actual line name, not the broad bit
                        refparam = param_name.replace(linename.replace("broad", "_"), line.ref)
                        if (par=="vel" or par=="fwhm_vel"):
                            confpar = conf[param_name.replace(linename, line.ref)]
                            error_value = get_error(confpar)
                            if np.isinf(error_value):
                                error_value = best_values[param_name.replace(linename, line.ref)].stderr
                        elif par=="amplitude":
                            print(f"ref line {refparam} for {linename}. Param name {param_name}")
                            refamplitude = best_values[refparam].value
                            amplitudefactor = "%s_factor" % linename
                            confpar = conf[amplitudefactor]
                            # the error is the amplitude error times the factor
                            error_value =  get_error(confpar) * refamplitude
                            if np.isinf(error_value):
                                error_value = best_values[amplitudefactor].stderr * refamplitude
                    # handle the broad factor
                    elif par == "fwhm_vel" and "broad" in linename:
                        confpar = conf["%s_fwhm_factor" % linename]
                        # multiply sigma of the narrow line with the factor
                        error_value = (
                            get_error(confpar)
                            * best_values[
                                "%s_%s" % (linename.replace("_broad", ""), par)
                            ]
                        )
                        
                    # normal line
                    elif param_name in conf:
                        confpar = conf[param_name]
                        error_value = get_error(confpar)
                        if np.isinf(error_value):
                            error_value = param.stderr

                elif param.stderr is not None:
                    # get the error from the line or from it's ref if tied
                    if line.ref in fit_lines:
                        refparam = param_name.replace(linename, line.ref)
                        if (par=="vel" or par=="fwhm_vel"):
                            error_value = best_values[refparam].stderr
                        elif par=="amplitude":
                            amplitudefactor = "%s_factor" % linename
                            error_value = best_values[refparam].value * best_values[amplitudefactor].stderr
                    else:
                        error_value = param.stderr

                param_names.append(f"{param_name}_err")
                param_values.append(error_value)

                if par == "vel":
                    # store for later
                    velocity = param.value
                    evelocity = error_value
                    wavelength, ewavelength = compute_shift(
                        line.wave,
                        z_sys=redshift,
                        velocity=velocity,
                        evelocity=evelocity,
                    )
                    param_names.append(f"{linename}_center")
                    param_values.append(wavelength)

                    param_names.append(f"{linename}_center_err")
                    param_values.append(ewavelength)

                if par == "amplitude" and args.distance is not None:
                    param_names.append("%s_lum" % linename)
                    param_values.append(param.value * 4 * pi * D**2)
                    param_names.append("%s_lum_err" % linename)
                    param_values.append(error_value * 4 * pi * D**2)

                if par == "fwhm_vel":
                    # FWHM corrected and FWHM velocity
                    FWHM_inst, eFWHM_inst = get_instrument_FWHM(wavelength, ewavelength)
                    FWHM_instvel = FWHM_inst / wavelength * ckms
                    if param.value > FWHM_instvel:
                        eFWHM_instvel = (
                            (eFWHM_inst / wavelength) ** 2.0
                            + (FWHM_inst * ewavelength / wavelength**2)
                        ) ** 0.5 * ckms
                        fwhm_vel_corrected, efwhm_vel_corrected = correct_FWHM(
                            param.value, FWHM_instvel, error_value, eFWHM_instvel
                        )
                        print(f"FWHM inst angstroms {FWHM_inst:.2f} {eFWHM_inst:.3f} {eFWHM_inst/FWHM_inst} {eFWHM_instvel/FWHM_instvel}")
                        print(f"param value {param.value:.3f} {error_value:.3f}, FWHM_inst {FWHM_instvel:.2f} {eFWHM_instvel:.3f}, FWHM_VEL corrected {fwhm_vel_corrected:.2f}, {efwhm_vel_corrected:.1f}, ewavelength {ewavelength:.2f}")

                        # f = A**2 - B **2
                        # df/dA = 2A
                        # df/dB = -2B
                        # df = sqrt( (df/da dA)**2 + (-2B dB)**2)
                        # df = sqrt( (2A dA)**2 + (2B dB)**2) = 2 sqrt(A dA + B dB)
                        # F = sqrt(f) = sqrt(A**2 - B**2)
                        # d sqrt(f) = 1/(2 * sqrt(f))  * df = 1 / (2 F) * df =   sqrt( (A dA)**2 + (B dB)**2) / F
                    else:
                        fwhm_vel_corrected = 0
                        efwhm_vel_corrected = np.nan

                    param_names.append(f"{linename}_fwhm_vel_corr")
                    param_values.append(fwhm_vel_corrected)
                    param_names.append(f"{linename}_fwhm_vel_corr_err")
                    param_values.append(efwhm_vel_corrected)

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
            error_value = get_error(confpar)
        elif param.stderr is not None:
            error_value = param.stderr

        param_names.append(f"cont_c{i}_err")
        param_values.append(error_value)

    # Write results to text file
    header = "\t".join(param_names)
    np.savetxt(
        f"{outpath}/{linegroupout}_fit_results.tsv",
        np.array(param_values).reshape(1, -1),
        header=header,
        fmt="%.6g",
        delimiter="\t",
    )
    logger.info("Output stored in %s" % outpath)
