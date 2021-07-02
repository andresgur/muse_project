# @Author: Andrés Gúrpide <agurpide>
# @Date:   29-06-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 29-06-2021

import numpy as np
from mpdaf.obj import Cube
import argparse

ap = argparse.ArgumentParser(description='Compute numer of nan and negative pixels in cube.')
ap.add_argument("input_cubes", nargs='+', help='List of input cubes', type=str)
args = ap.parse_args()

for input_cube in args.input_cubes:
    cube = Cube(input_cube)
    nspaxels = cube.shape[0] * cube.shape[1] * cube.shape[2]
    print("\t### Stats for Cube %s ###" % input_cube)
    print("\t-------------------------------------------")
    print("\tNumber of spaxels: %d" % (nspaxels))
    nnans = np.count_nonzero(np.isnan(cube.data))
    per_nans = nnans / nspaxels * 100
    print("\tNumber of NaNs: %d (%.1f %%)" % (nnans, per_nans))
    negative = np.count_nonzero(cube.data < 0)
    per_neg = negative / nspaxels * 100
    print("\tNumber of negative values: %d (%.1f %%)" % (negative, per_neg))
