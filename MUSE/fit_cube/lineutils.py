#!/usr/bin/env python
# -*- coding: utf-8 -*
import numpy as np
ckms = 299792.458
sigma_to_fwhm = 2.355  # FWHM = 2.355 * sigma



def peculiar_velocity(reference, wavelength, ewavelength=0, z_sys=0):
    """Calculate the peculiar velocity of a line given its reference wavelength and observed wavelength. See https://academic.oup.com/mnras/article/442/2/1117/983284 Eq 12.
    Works for relativistic velocities too. All wavelengths in same units.
    
    Parameters
    ----------
    reference: float
        Reference wavelength
    wavelength: float
        Observed wavelength
    ewavelength: float
        Uncertainty in observed wavelength
    z_sys: float
        Systemic redshift

    Returns
    -------
    float or array-like
        Peculiar velocity in km/s
    float or array-like
        Uncertainty in peculiar velocity in km/s
    """
    z = wavelength / reference - 1
    ez = ewavelength / reference
    zp = (z - z_sys) / (1 + z_sys)
    ezp = ez / (1 + z_sys)
    factor = (1 + zp)**2
    efactor = 2. * (1 + zp) * ezp
    vp = ckms * (factor - 1) / (factor + 1)
    #df/factor = 1 x (factor + 1) - (factor -1)X1 / (factor + 1)^2 ==> 2 / (Factor + 1)^2
    evp = ckms * 2. * efactor / (factor + 1)**2.
    return vp, evp


def get_instrument_FWHM(wavelength, ewavelength=0):
    """Get MUSE FWHM LSF based on Bacon et al. 2017 Equation 8 (see also Benoit et al. 2018) 1
    Parameters
    ----------
    wavelength: float or array-like
        Wavelength in Angstroms

    Returns
    ------- 
    fwhm: float or array-like
        FWHM in Angstroms
    efwhm: float or array-like
        Uncertainty in FWHM in Angstroms
    """
    A = 5.866 * 1e-8
    B = 9.187 * 1e-4
    return A * wavelength ** 2. - B * wavelength + 6.04, np.abs(2 * A * wavelength - B * ewavelength)
