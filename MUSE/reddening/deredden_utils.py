from extinction import ccm89, remove
import re
import subprocess
import numpy as np
from math import log
import astropy.units as u
import os

HOME = os.getenv("HOME")
PATH = "%s/scripts/utils" % HOME

def galactic_extinction(wavelengths, fluxes, efluxes, EBV_gal, EBV_gal_err=0, Rv=3.1):
    """Wavelengths in angstroms
    wavelength: 2D array
    fluxes: 2D array
        Same shape as wavelength
    efluxes: 2D array,
        Same shape as wavelength and fluxes
    EBV_gal: float,
    EBV_gal_err: float
        Error on EBV_gal
    Rv: float,
        3.1 by Default

    Returns
    -------
    deredden_fluxes: array-like
        Same shape as fluxes
    varredden_fluxes: array-like
        The propagated error of the dereddened fluxes (i.e. sqrt(variance)
    """
    Av = EBV_gal * Rv
    Av_gal_err = EBV_gal_err * Rv
    wavs = np.array(wavelengths.flatten(), dtype="double")
    mag_ext = ccm89(wavs, Av, Rv, unit="aa")
    flattened_fluxes = fluxes.flatten()
    deredden_fluxes = remove(mag_ext, flattened_fluxes)
    ederedden_fluxes = (remove(mag_ext, efluxes.flatten())**2 + (Av_gal_err * flattened_fluxes * log(10) * 0.4)**2) **0.5
    return deredden_fluxes.reshape(wavelengths.shape), ederedden_fluxes.reshape(wavelengths.shape)


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
    mean_pattern = re.compile(r"E\(B-V\) mean \(Schlafly & Finkbeiner 2011\) over \d+ degrees: ([0-9.]+) mag")
    std_pattern = re.compile(r"E\(B-V\) std \(Schlafly & Finkbeiner 2011\) over \d+ degrees: ([0-9.]+) mag")

    mean_match = mean_pattern.search(output)
    std_match = std_pattern.search(output)

    if mean_match and std_match:
        mean_value = float(mean_match.group(1))
        std_value = float(std_match.group(1))
        return mean_value, std_value
    else:
        print("Error: Could not parse the mean or std values from the output.")
        return None



class C00:
    """Calzetti extinction curve"""
    def __init__(self, Rv=4.05):
        self.Rv = Rv

    def evaluate(self, wavelength):
        """wavelength must be in wavelength units"""
        micron = wavelength.to(u.micron).value
        micron = np.asarray(micron)
        if micron.ndim == 0:
            micron = micron[None]
        x = 1 / micron
        optical_indx = np.where(np.logical_and(0.63 <= micron, micron <= 2.20))
        ir_indx = np.where(np.logical_and(0.12 <= micron, micron <= 0.63))
        k = np.empty(len(micron))
        k[optical_indx] = 2.659 * (-1.857 + 1.040 * x[optical_indx]) + self.Rv
        k[ir_indx] = 2.659 * (-2.156 + 1.509 * x[ir_indx] - 0.198 * x[ir_indx]**2 + 0.011 * x[ir_indx]**3) + self.Rv
        
        # If input was scalar, return scalar
        if wavelength.to(u.micron).value.ndim == 0:
            return k[0]
        return k
