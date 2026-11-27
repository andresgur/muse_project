#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Collect fit results from batch jobs and create FITS maps.
"""

import numpy as np
import argparse
import os
import json
import pickle
import glob
import re
import logging
from mpdaf.obj import Cube
from lineutils import correct_FWHM, get_instrument_FWHM_vspace, compute_shift
from fit_cube import create_out_maps, save_image, get_error
from line import Lines

CATALOG_LINES = Lines().lines
ckms = 299792.458

# Configure logging for SLURM stderr
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def calculate_line_noise(center, sigma, residuals, wavelengths, window=7):
    """Calculate the noise underlying a particular line
    This is compute as the rms of the residuals over a 2sigma window, encompassing approximately 95% of the line profile
    """
    left, right = center - window * sigma, center + window * sigma
    mask = (wavelengths > left) & (wavelengths < right)
    noise = np.std(residuals[mask])
    return noise


def collect_results(spectra_dir, results_pattern, noisewindow):
    """Collect all result files and create maps"""
    outdir = spectra_dir + "/" + "results"
    os.makedirs(outdir, exist_ok=True)

    # Find all result files
    line_result_files = glob.glob(spectra_dir + "/*" + results_pattern + "*.pkl")
    if not line_result_files:
        logger.error(f"No result files found matching pattern {results_pattern} in {spectra_dir}")
        return
    print(f"Found {len(line_result_files)} result files")
    results_stem = re.escape(results_pattern)
    # note we already insert the preceding "_" if not given in the stem
    results_regex = re.compile(rf"^(?P<linegroups>.+?){results_stem}_\d+_\d+\.pkl$")

    matches = [results_regex.match(os.path.basename(f)) for f in line_result_files]
    linegroups = np.unique([m.group("linegroups") for m in matches if m is not None])

    # Load all results
    all_results = []
    
    try:
        # each linegroup may contain several files if the batch jobs were split, so we loop over linegroups and then over files for each linegroup
        for linegroup in linegroups:

            # Load manifest for this group of lines
            with open(os.path.join(spectra_dir, f"manifest_{linegroup}.json")) as f:
                manifestlines = json.load(f)
            print(manifestlines)

            data_cube = Cube(manifestlines["input_subcube"])
            redshift = data_cube.primary_header["REDSHIFT"]
            margin = data_cube.primary_header["WINDOW"]
            degree = " ".join(f"{i}" for i in manifestlines["degrees"])
            try:
                FWHM_min = manifestlines["FWHM_min_factor"]
            except KeyError:
                FWHM_min = None
            try: 
                FWHM_max = manifestlines["FWHM_max_factor"]
            except KeyError:
                FWHM_max = None

            saveimagedict = {'redshift': redshift, 'FWHM_min': FWHM_min, 
                                             'margin': margin, 'degree': degree, 
                                             'noisewindow': noisewindow}

            # dummy skeleton for the output cube
            templatefile = data_cube[0, :, :]

            outshape = templatefile.shape

            print(f"Creating output maps...")

            # unpack lines, note we cannot simply split by "_" because some lines have "_" in their names (e.g. OI_sky)
            fit_lines = [line for group in manifestlines["linegroups"] for line in group]
            outmaps = create_out_maps(
                fit_lines, outshape, degree=max(manifestlines["degrees"])
            )
            # grab all result files -- important, do not add "*" otherwise we may get extra files for a given group line
            # e.g. if HeII4686_results*.pkl and HeII4686_ArIV4711_ArIV4740_results*.pkl exist, we want to make sure we only grab the files for the current linegroup
            line_results_file = glob.glob(
                spectra_dir + f"/*{linegroup}" + results_pattern + "*.pkl"
            )
            print(f"Found {len(line_results_file)} result files for linegroup {linegroup}")
            
            for rfile in line_results_file:

                print(f"Processing result file {rfile}")

                with open(rfile, "rb") as f:
                    fits_data = pickle.load(f)
                # skip failed fits
                validfits = [f for f in fits_data if f["success"]]

                print(f"Processing {len(validfits)} valid fits...")

                for fit_data in validfits:
                    try:
                        x, y = fit_data["x"], fit_data["y"]
                        outmaps["redchi"][y, x] = fit_data["redchi"]
                        outmaps["BIC"][y, x] = fit_data["bic"]
                        outmaps["rsquared"][y, x] = fit_data["rsquared"]
                        best_values = fit_data["params"]
                        conf = fit_data["conf"]

                        for linename in fit_lines:
                        
                            noise = calculate_line_noise(
                                best_values["%s_center" % (linename)]["value"],
                                best_values["%s_sigma" % (linename)]["value"],
                                fit_data["model"] - fit_data["data"],
                                fit_data["wave"],
                                window=noisewindow
                            )

                            # SNR calculation
                            line = CATALOG_LINES[linename]
                            #signal_param_name = f"{linename}_amplitude"
                            signal_param_name = f"{linename}_height"
                            #if line.ref is not None: --> NO! We always get the height from its line, even if it had a ref, otherwise it will be smaller/higher
                            # get the amplitude for the ref line, the error is already reduced by fititng the two lines together, so they technically ahve the same SNR
                                #signal_param_name = f"{line.ref}_amplitude"
                            #    signal_param_name = f"{line.ref}_height"
                        
                            signalparam = best_values[signal_param_name]
                        
                            error = signalparam["stderr"] if signalparam["stderr"] is not None else np.nan
                            #if conf is not None:
                            #       error = get_error(conf[signal_param_name])
                            # noise = fit_data["noise"]
                            signal = signalparam["value"]
                            #noise = error
                            #outmaps["%s_snr" % (linename)][y, x] = (
                            #    signal / noise if noise > 0 else np.nan
                            #)
                            outmaps["%s_snr" % (linename)][y, x] = signal / noise 

                            for par in ["amplitude", "vel", "fwhm_vel"]:
                                param_name = "%s_%s" % (linename, par)
                                if param_name in best_values:
                                    param = best_values[param_name]
                                    paramstderr = (
                                        param["stderr"]
                                        if param["stderr"] is not None
                                        else np.nan
                                    )
                                    outmaps[param_name][y, x] = param["value"]

                                    # Store parameter uncertainties
                                    # get the uncertainty from the confidence intervals if available
                                    if conf is not None:
                                        # handle FWHM specially, as it is stored as sigma
                                        # handle the broad factor
                                        if par == "fwhm_vel" and "broad" in linename:
                                            confpar = conf["%s_fwhm_factor" % linename]
                                            # multiply sigma of the narrow line with the factor
                                            errorbroadfactor = get_error(confpar)
                                            errorbroadfactor = (
                                                errorbroadfactor
                                                if np.isfinite(errorbroadfactor)
                                                else paramstderr
                                            )
                                            fwhmvelnarrowline = best_values[
                                                "%s_%s" % (linename.replace("_broad", ""), par)
                                            ]
                                            outmaps[f"{linename}_e{par}"][y, x] = (
                                                errorbroadfactor * fwhmvelnarrowline["value"]
                                            )
                                        # if we have conf but some params are missing they are tied, so we just skip them if they are not in the conf interval
                                        elif param_name in conf:

                                            confpar = conf[param_name]
                                            error = get_error(confpar)
                                            outmaps["%s_e%s" % (linename, par)][y, x] = (
                                                error if np.isfinite(error) else paramstderr
                                            )
                                        # keep the stderr for tied parameters if available
                                        # elif param.stderr is not None:
                                        #   outmaps["%s_e%s" % (line.name, par)][y, x] = param.stderr

                                    # otherwise, use the standard error
                                    elif paramstderr is not None:
                                        outmaps[f"{linename}_e{par}"][y, x] = paramstderr

                        # continuum parameters now
                        for i in range(2):
                            contparam = f"cont_c{i}"
                            if contparam in best_values:
                                param = best_values[contparam]
                                outmaps[contparam][y, x] = param["value"]

                                if conf is not None:
                                    confpar = conf[contparam]
                                    error = get_error(confpar)
                                    outmaps["cont_ec%d" % i][y, x] = (
                                        error if np.isfinite(error) else param["stderr"]
                                    )
                                elif paramstderr is not None:
                                    outmaps["cont_ec%d" % i][y, x] = param["stderr"]
                    except KeyError as exc:
                        params = fit_data.get("params")
                        conf = fit_data.get("conf")
                        param_keys = sorted(params) if isinstance(params, dict) else params
                        conf_keys = sorted(conf) if isinstance(conf, dict) else conf
                        raise KeyError(
                            f"Missing key {exc.args[0]!r} while processing line group {linegroup}, "
                            f"file {rfile}, pixel x={fit_data.get('x', 'unknown')}, "
                            f"y={fit_data.get('y', 'unknown')}. "
                            f"Available params: {param_keys}; confidence interval keys: {conf_keys}"
                        ) from exc

            linepars = ["amplitude", "fwhm_vel", "vel", "snr"]

            for linename in fit_lines:
                line = CATALOG_LINES[linename]
                # store line param values
                for par in linepars:
                    # parameter value
                    save_image(
                        data=outmaps["%s_%s" % (linename, par)],
                        templatefile=templatefile,
                        **saveimagedict,
                        filename="%s/%s_%s.fits" % (outdir, linename, par),
                    )
                    # parameter uncertainty
                    if par == "snr":
                        # do not store the uncertainty for SNR
                        continue
                    save_image(
                        data=outmaps["%s_e%s" % (linename, par)],
                        templatefile=templatefile,
                        **saveimagedict,
                        filename="%s/%s_e%s.fits" % (outdir, linename, par),
                    )

                # store wavelenght
                wavelength, ewavelength = compute_shift(
                    line.wave,
                    z_sys=redshift,
                    velocity=outmaps["%s_%s" % (linename, "vel")],
                    evelocity=outmaps["%s_e%s" % (linename, "vel")],
                )

                par = "center"
                save_image(
                    data=wavelength,
                    templatefile=templatefile,
                    **saveimagedict,
                    filename=f"{outdir}/{linename}_{par}.fits",
                    bunit="Angstrom",
                )

                # ecenter
                save_image(
                    data=ewavelength,
                    templatefile=templatefile,
                    **saveimagedict,
                    filename=f"{outdir}/{linename}_e{par}.fits",
                    bunit="Angstrom",
                )

                # FWHM
                FWHM_o = outmaps["%s_fwhm_vel" % (linename)]
                eFWHM_o = outmaps["%s_efwhm_vel" % (linename)]
                # Instrumental FWHM
                FWHM_instvel, eFWHM_instvel = get_instrument_FWHM_vspace(wavelength, ewavelength)

                # get pixels where FWHM_obs > FWHM_inst
                valid_selector = FWHM_o > FWHM_instvel

                par = "fwhm_vel_corrected"
                FWHM_corrected, eFWHM_corrected = correct_FWHM(
                    FWHM_o[valid_selector],
                    FWHM_instvel[valid_selector],
                    eFWHM_o[valid_selector],
                    eFWHM_instvel[valid_selector],
                )
                outmaps[f"{linename}_{par}"][valid_selector] = FWHM_corrected
                outmaps[f"{linename}_{par}"][~valid_selector] = 0
                # store FWHM corrected maps
                save_image(
                    data=outmaps[f"{linename}_{par}"],
                    templatefile=templatefile,
                    **saveimagedict,
                    filename="%s/%s_%s.fits" % (outdir, linename, par),
                    bunit="km/s",
                )

                outmaps[f"{linename}_e{par}"][
                    valid_selector
                ] = eFWHM_corrected  # the non valid valus are simply nan
                save_image(
                    data=outmaps[f"{linename}_e{par}"],
                    templatefile=templatefile,
                    **saveimagedict,
                    filename="%s/%s_e%s.fits" % (outdir, linename, par),
                    bunit="km/s",
                )

            for i in range(2):
                bunit = templatefile.unit if i == 0 else f"{templatefile.unit}/Angstrom^{i}"
                save_image(
                    outmaps["cont_c%d" % i],
                    templatefile=templatefile,
                    **saveimagedict,
                    filename="%s/cont_c%d_%s.fits" % (outdir, i, linegroup),
                    bunit=bunit,
                )
                save_image(
                    outmaps["cont_ec%d" % i],
                    templatefile=templatefile,
                    **saveimagedict,
                    filename="%s/cont_ec%d_%s.fits" % (outdir, i, linegroup),
                    bunit=bunit,
                )

            # save stat files
            for fitstat in ["redchi", "rsquared", "BIC"]:
                # write the reduced chi square and r squared maps
                save_image(
                    data=outmaps[fitstat],
                    templatefile=templatefile,
                    **saveimagedict,
                    filename="%s/%s_%s.fits" % (outdir, fitstat, linegroup),
                )

        print(f"Processed lines {linegroup}")
    except Exception:
        logger.error(f"Failed to process line group: {linegroup}")
        raise

    return all_results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("spectra_dir", help="Directory with split spectra")
    parser.add_argument(
        "results_pattern",
        help="Pattern for result files. The script will look for files named '*results_pattern.pkl' in the spectra_dir.",
    )
    parser.add_argument('-n', '--noisewindow', help="Specify the noise window for the signal-to-noise ratio calculation. Normally don't change, but you can experiment. " \
    "Values below 5 are not recommended. Values above 10 will cause other lines to get into the residual error calculation. Default is 7.", type=float, default=7)
    args = parser.parse_args()
    if args.results_pattern.endswith(".pkl"):
        print("Error: results_pattern cannot end with .pkl")
        exit(1)
    elif not args.results_pattern.startswith("_"):
        results_pattern = "_" + args.results_pattern

    collect_results(args.spectra_dir, results_pattern, args.noisewindow)
