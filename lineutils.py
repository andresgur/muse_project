#!/usr/bin/env python
# -*- coding: utf-8 -*
import numpy as np
from scipy.stats import norm
ckms = 299792.458
sigma_to_fwhm = 2.355  # FWHM = 2.355 * sigma
# LSF based on Bacon et al. 2017 Equation 8 (see also Benoit et al. 2018)
LSF_A = 5.866e-8
LSF_B = 9.187e-4
LSF_C = 6.04

def compute_shift(reference, z_sys, velocity, evelocity=0):
    """Compute the observed wavelength of a line given its reference wavelength, systemic redshift and peculiar velocity. 
    See https://academic.oup.com/mnras/article/442/2/1117/983284 Eq 8 and 9.
    Works for relativistic velocities too.
    
    Parameters
    ----------
    reference: float
        Reference wavelength
    z_sys: float
        Systemic redshift
    velocity: float or array-like
        Peculiar velocity in km/s

    Returns
    -------
    float or array-like
        Observed wavelength in same units as reference
    """
    v_c = velocity / ckms
    evc = evelocity / ckms
    zp_plus_1 = ((1. + v_c) / (1. - v_c)) **0.5 # this is 1 + zp
    # d(1 +v/(1-v)) = (1/(1-v) + (1+v)/(1-v)^2) = 1/(1-v) [1 + (1+v)/(1-v)] = 1/(1-v) * (1-v +1 +v)/(1-v) = 2/(1-v)^2
    #edivision = 2 * evc / (1 - v_c)**2
    # 0.5 * 2 cancel out
    #ezp_plus_1 = 0.5 * edivision / zp_plus_1
    # so simplyfing:
    ezp_plus_1 = evc / (1. - v_c)**2. / zp_plus_1
    # 1 + z = (1 + zp) * (1 + z_sys) Eq 9 from https://doi.org/10.1093/mnras/stu920
    sysshifted = reference * (1. + z_sys)
    values, errors = sysshifted * zp_plus_1, sysshifted * ezp_plus_1
    return values, errors



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
    """Get MUSE FWHM LSF based on Bacon et al. 2017 Equation 8 (see also Benoit et al. 2018)

 
    FWHM_inst = LSF_A * wavelength ** 2. - LSF_B * wavelength + LSF_C
    eFWHM_inst = np.abs(2 * LSF_A * wavelength - LSF_B ) * ewavelength
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
    return LSF_A * wavelength ** 2. - LSF_B * wavelength + LSF_C, np.abs(2 * LSF_A * wavelength - LSF_B ) * ewavelength


def get_instrument_FWHM_vspace(wavelength, ewavelength=0):
    """Get MUSE FWHM LSF based on Bacon et al. 2017 Equation 8 (see also Benoit et al. 2018) in velocity space
    This is just the FWHM divided by the wavelength, converted to velocity units (km/s) (FWHM/lambda x c).


    Parameters
    ----------
    wavelength: float or array-like
        Wavelength in Angstroms

    Returns
    ------- 
    fwhm: float or array-like
        FWHM in km/s
    efwhm: float or array-like
        Uncertainty in FWHM in km/s
    """
    fwhm, _ = get_instrument_FWHM(wavelength)
    fwhm_vel = ckms * fwhm / wavelength
    # this is the derivative of c FWHM/lambda = c x (LSF_A * wavelength - LSF_B + LSF_C / wavelength)
    efwhm_vel = ckms * np.abs(LSF_A - LSF_C / wavelength**2) * ewavelength
    return fwhm_vel, efwhm_vel



def correct_FWHM(FWHM_obs, FWHM_inst, eFWHM_obs, eFWHM_inst):
    """Correct observed FWHM by instrumental FWHM. All in same units.

    #d(A^2 - B^2) = sqrt((2AdA)^2 - (2B dB)^2) --> 2 x sqrt(Ada^2 + Bdb^2)
    #eFWHM f = sqrt(A^2-B^2); df/dA = 2A/2f = A/f; df/dB = -2B/2f = -B/f; eFWHM = sqrt((eA/f)^2 + (eB/f)^2) = sqrt(eA^2 + eB^2)/f
    
    
    """
    
    FWHM_vel_corrected = (FWHM_obs**2. - FWHM_inst**2.)**0.5
    eFWHM_vel_corrected = ((FWHM_obs * eFWHM_obs)**2. + (FWHM_inst * eFWHM_inst)**2.)**0.5 / FWHM_vel_corrected
    return FWHM_vel_corrected, eFWHM_vel_corrected



def get_SFR(LHalpha):
    """"Returns the SFR in M_sun/yr from the Halpha luminosity in erg/s (Equation 2 from Kennicutt 1998, ApJ, 498, 541).

    Parameters
    ----------
    LHalpha: float or array-like
        The (extinction-corrected) Halpha luminosity (can also be its uncertainty, in which case the uncertainty is what is calculated)
    """
    return 7.9 * 1e-42 * LHalpha

def get_EW(Fline, linewave, result, eFline=0, cont_prefix="cont"):
    """Compute the equivalent width (and its propagated uncertainty, accounting for flux and continuum uncertainties) of a line given its flux and continuum level. The continuum is estimated from the fit result at the line wavelength.
    
    Parameters
    ----------
    Fline: float or array-like
        The line flux
    linewave: float or array-like
        The line wavelength
    result: lmfit.model.ModelResult
        The fit result from which to evaluate the continuum level and its uncertainty
    eFline: float or array-like, optional
        The uncertainty in the line flux (default is 0)
    cont_prefix: str, optional
        The prefix of the continuum component in the fit result (default is "cont")
    """
    cont = result.eval_components(x=linewave)[f"{cont_prefix}"]
    if result.errorbars:
        result.eval_uncertainty(sigma=1, x=linewave)
        econt = result.dely_comps[f"{cont_prefix}"]
    else:
        print(f"Warning: No error bars available for continuum at {linewave}. Setting econt=0.")
        econt = 0
    EW = Fline / cont
    eEW = ((eFline / cont)**2. + (Fline * econt / cont**2.)**2.) ** 0.5
    return EW, eEW


