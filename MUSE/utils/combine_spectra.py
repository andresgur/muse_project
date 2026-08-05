# combine_spectra.py
"""Combine N MUSE/mpdaf spectra via sum or mean."""

import argparse
import sys
from mpdaf.obj import Spectrum


def combine_spectra(input_files, mode):
    spectra = [Spectrum(f) for f in input_files]

    result = spectra[0].copy()
    for sp in spectra[1:]:
        if sp.get_step()==spectra[0].get_step() and sp.get_start()==spectra[0].get_start() and sp.get_end()==spectra[0].get_end():
            result = result + sp
        elif sp.get_step()==spectra[0].get_step() and (sp.get_start()-spectra[0].get_start()) < 0.1  and (sp.get_end()-spectra[0].get_end()) < 0.1:
            sp.set_wcs(wave=spectra[0].wave)
            result = result + sp
        else:
            raise ValueError(f"Spectra {spectra[0].filename} and {sp.filename} have different wavelength grids.")

    if mode == "mean":
        result = result / len(spectra)

    return result


def main():
    ap = argparse.ArgumentParser(description="Combine N mpdaf spectra.")
    ap.add_argument("inputs", nargs="+", help="Input spectrum FITS files")
    ap.add_argument("-m", "--mode", choices=["sum", "mean"], default="mean",
                    help="Combination mode. Default: mean")
    ap.add_argument("-o", "--output", default="combined.fits",
                    help="Output FITS file. Default: combined.fits")
    args = ap.parse_args()

    if len(args.inputs) < 2:
        print("Error: provide at least 2 input spectra.", file=sys.stderr)
        sys.exit(1)

    print(f"Combining {len(args.inputs)} spectra (mode={args.mode})")
    print(f"Input files: {args.inputs}")
    combined = combine_spectra(args.inputs, args.mode)
    # add keywords to header
    combined.primary_header['COMBMODE'] = (args.mode, 'Combination mode: sum or mean')
    combined.primary_header['NCOMBINE'] = (len(args.inputs), 'Number of combined spectra')

    outfile = f"{args.mode}_{args.output}"
    combined.write(outfile)
    print(f"Written to {outfile}")


if __name__ == "__main__":
    main()