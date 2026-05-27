#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Fit spectra directly from one subcube file.
Metadata (linegroup/redshift/window) is read from subcube FITS headers.
"""

import argparse
import json
import os
import pickle
import time
from itertools import islice

import numpy as np
from astropy.io import fits
from mpdaf.obj import Cube

from fit_cube import parse_linegroups
from fitutils import DofError, SpectrumMaskedError, fit_spectrum
from line import Lines
from tqdm import tqdm

CATALOG_LINES = Lines().lines


def load_subcube_metadata(subcube_file):
    with fits.open(subcube_file) as hdul:
        header = hdul[0].header
        linegroup_raw = header.get("LINEGRP")
        redshift = header.get("REDSHIFT")
        window = header.get("WINDOW")

    if linegroup_raw is None:
        raise ValueError(
            "LINEGRP keyword not found in FITS header. Recreate subcube with split_cube_subcube.py"
        )
    if redshift is None:
        raise ValueError(
            "REDSHIFT keyword not found in FITS header. Recreate subcube with split_cube_subcube.py"
        )
    if window is None:
        window = 25.0

    linegroups = [group for group in str(linegroup_raw).split() if group]
    if len(linegroups) == 0:
        raise ValueError("LINEGRP header is empty. Cannot determine lines to fit")

    return linegroups, float(redshift), float(window)


def crop_and_mask_spectrum_to_windows(spectrum, waveinit, waveend, wav_cuts):
    wave = spectrum.wave.coord()
    step = wave[1] - wave[0]
    selector = (wave > (waveinit - step * 1.01)) & (wave < (waveend + step * 1.01))

    wave_data = wave[selector]
    data_data = spectrum.data.data[selector]
    var_data = spectrum.var.data[selector] if spectrum.var is not None else None
    original_mask = spectrum.mask[selector]
    final_mask = np.ones_like(original_mask, dtype=bool)

    for lmin, lmax in wav_cuts:
        inrange = (wave_data >= lmin) & (wave_data <= lmax)
        final_mask[inrange] = original_mask[inrange]

    return {
        "wave": wave_data,
        "data": data_data,
        "mask": final_mask,
        "var": var_data,
    }


class SpectrumObj:
    def __init__(self, wave, data, var, mask=None):
        self.wave = wave
        self.data = np.ma.MaskedArray(data, mask=mask)
        self.var = var
        self.mask = mask


def fit_single_spectrum(
    spectrum,
    fit_lines,
    redshift,
    wavelengths,
    degrees=[1],
    sigma=1.4,
    uncertainties=False,
    x=None,
    y=None,
    deltaBIC_threshold=6,
):
    best_result = None
    conf = None
    best_degree = -1
    best_bic = np.inf
    # we calculate uncertainties only for the best model later on
    for degree in degrees:
        try:
            # we compute uncertainties later on --> note this messes up lmfit internally, so best to do it for every fit as it seems some variables are stored internally
            result, conf = fit_spectrum(
                spectrum,
                fit_lines,
                redshift,
                sigma,
                wavelengths,
                degree,
                uncertainties,
            )
            # set the best result to degree 0 by default
            if best_result is None:
                best_degree = degree
                best_bic = result.bic
                best_result = result
            elif result.success and result.bic + deltaBIC_threshold < best_bic:
                best_bic = result.bic
                best_result = result

        except SpectrumMaskedError as err:

            fit_result = {
                "x": int(x),
                "y": int(y),
                "success": False,
                "error": str(err),
            }
            return fit_result, -1
        except DofError as err:
            # if the degree = 0 worked, then we just break as the others will all give dof problems
            if best_result is not None:
                break
            else:
                # if it's the first time we just return the fit_result with the error
                fit_result = {
                    "x": int(x),
                    "y": int(y),
                    "success": False,
                    "error": str(err),
                }
                return fit_result, -1

    fit_result = {
        "x": int(x),
        "y": int(y),
        "success": True,
        "params": {
            k: {"value": v.value, "stderr": v.stderr}
            for k, v in best_result.params.items()
        },
        "redchi": best_result.redchi,
        "bic": best_result.bic,
        "rsquared": best_result.rsquared,
        "ndata": best_result.ndata,
        "nvarys": best_result.nvarys,
        "ier": best_result.ier,
        "errorbars": best_result.errorbars,
        "model": best_result.best_fit,
        "data": best_result.data,
        "wave": result.userkws["x"],
        "conf": conf,
    }

    return fit_result, best_degree


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_subcube", help="Path to subcube FITS file")
    parser.add_argument("-s", "--start-idx", type=int, default=0)
    parser.add_argument("-e", "--end-idx", type=int, default=None)
    parser.add_argument(
        "-d",
        "--degrees",
        type=int,
        default=[1],
        nargs="+",
        help="Polynomial degree(s) for continuum fit to try",
    )
    parser.add_argument("--sigma", type=float, default=1.4)
    parser.add_argument("--uncertainties", action="store_true")
    parser.add_argument(
        "-o", "--output", help="Output root for pickle results", required=True
    )

    args = parser.parse_args()
    if not args.input_subcube.endswith(".fits"):
        raise ValueError("Input subcube must be a FITS file")
    linegroupsstr, redshift, window = load_subcube_metadata(args.input_subcube)

    fit_lines, linegroups, _, _, _ = parse_linegroups(linegroupsstr, redshift, window)

    linegroup_header = " ".join(linegroupsstr)

    print(
        f"Using subcube metadata: linegroups={linegroup_header}, redshift={redshift:.6f}, window={window:.2f}"
    )

    data_cube = Cube(args.input_subcube)
    total_spectra = data_cube.shape[1] * data_cube.shape[2]
    start_idx = max(0, args.start_idx)
    end_idx = (
        total_spectra if args.end_idx is None else min(args.end_idx, total_spectra)
    )
    if start_idx >= end_idx:
        raise ValueError(
            f"Invalid index range [{start_idx}, {end_idx}) for {total_spectra} spectra"
        )
    print(f"Using spectral index range: [{start_idx}, {end_idx})")
    totalspec = end_idx - start_idx
    print(f"Total spectra to fit: {totalspec}")

    results = []

    starttime = time.time()

    for y, x in tqdm(
        islice(np.ndindex(*data_cube.shape[1:]), start_idx, end_idx),
        total=totalspec,
        mininterval=1,
        unit="spectrum",
    ):
        spectrum = data_cube[:, y, x]
        if np.all(spectrum.mask):
            continue

        best_result, _ = fit_single_spectrum(
            spectrum=spectrum,
            fit_lines=fit_lines,
            redshift=redshift,
            wavelengths=spectrum.wave.coord(),
            degrees=args.degrees,
            sigma=args.sigma,
            uncertainties=args.uncertainties,
            x=x,
            y=y,
        )

        results.append(best_result)

    endtime = time.time()
    ellapsed_time = endtime - starttime

    fitted_count = len(results)
    if fitted_count > 0:
        print(
            f"Fitted {fitted_count}/{totalspec} spectra in {ellapsed_time:.2f} seconds ({ellapsed_time/3600:.2f} hours) ({ellapsed_time/fitted_count:.2f} seconds/spectrum)"
        )
    else:
        print("No spectra were fitted in requested index range (likely all masked)")

    linegroup_out = "_".join(
        [line for group in linegroupsstr for line in group.split(",")]
    )

    if args.output.endswith(".pkl"):
        outroot = args.output.replace(".pkl", "")
    else:
        outroot = args.output

    output_dir = os.path.dirname(args.input_subcube) or "."

    outfile = os.path.join(
        output_dir,
        f"{linegroup_out}_{outroot}_{start_idx}_{end_idx}.pkl",
    )
    with open(outfile, "wb") as fobj:
        pickle.dump(results, fobj)

    manifestfit = {
        "input_subcube": args.input_subcube,
        "linegroups": linegroups,
        "degrees": args.degrees,
        "sigma": args.sigma,
        "uncertainties": args.uncertainties,
        "start_idx": start_idx,
        "end_idx": end_idx,
        "fitted_count": fitted_count,
        "ellapsed_time": f"{ellapsed_time /3600:.2f}",
    }

    manifest_name = f"manifest_{linegroup_out}.json"
    manifest_path = os.path.join(output_dir, manifest_name)
    with open(manifest_path, "w") as fobj:
        json.dump(manifestfit, fobj, indent=2)

    print(f"Saved {len(results)} results to {outfile}")
    print(f"Saved fit manifest to {manifest_path}")


if __name__ == "__main__":
    main()
