#!/usr/bin/env python
# -*- coding: utf-8 -*
import numpy as np
import logging, os
import sys
sys.path.append("/home/andresgur/scripts/pythonscripts/maxime_muse/fitting")
from lineutils import get_instrument_FWHM, sigma_to_fwhm, peculiar_velocity, ckms, compute_shift
from astropy.stats import sigma_clipped_stats
from lmfit.models import PolynomialModel, GaussianModel
import warnings
import matplotlib.pyplot as plt


if os.path.isfile("/home/andresgur/.config/matplotlib/stylelib/paper.mplstyle"):
    plt.style.use("/home/andresgur/.config/matplotlib/stylelib/paper.mplstyle")
from math import pi


logger = logging.getLogger('fit_utils')


class SpectrumMaskedError(Exception):
    """Raised when the spectrum is completely masked and cannot be fitted."""
    pass


class DofError(Exception):
    """Raised when there are insufficient degrees of freedom for fitting (too few data points for the number of parameters)."""
    pass

def setup_logging(logfile, loglevel=logging.DEBUG):
    handlers = [logging.StreamHandler()]
    if loglevel==logging.DEBUG:
        handlers.append(logging.FileHandler(logfile))
    """Set up logging configuration"""
    logging.basicConfig(
        level=loglevel,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=handlers
    )

sqrt2pi = (2. * pi) ** 0.5



def get_error(confpar):
    """Get the error from confidenece interval output from lmfit. Return the average if both intervals are valid, otherwise return the valid one.
    Parameters
    ----------
    conf: tuple
        The tuple containing the confidence interval information from lmfit and the best fit value
        """
    bestpar = confpar[1][1]
    negative_error =  (bestpar - confpar[0][1])
    positive_error =  (confpar[2][1] - bestpar)
    if np.isinf(negative_error) and np.isinf(positive_error):
        return np.inf
    elif np.isinf(negative_error):
        return positive_error
    elif np.isinf(positive_error):
        return negative_error
    else:
        return (negative_error + positive_error) / 2.


def get_initial_fit(model, x):
    """
    Parameters
    ----------
    model: lmfit.Model
        The lmfit model to evaluate
    x: array-like
        The x values to evaluate the model at

    Returns
    -------
    array-like
        The evaluated model values
    """
    params = model.make_params()
    return model.eval(params, x=x)


def create_line_model(lines, redshift=0.0, sigma=1.4, degree=1, vmargin=400., cont_prefix="cont"):
    """Create lmfit model with all the lines to be fitted. Lines sigmas and position are tight to the reference values of other lines.
    Parameters
    ---------
    lines: dict,
        Lines to be fitted with names as keys
    redshift: float
        Initial guess for the redshift
    sigma: float
        Initial guess for the line sigma in Angstroms
    vmargin: float
        +-Margin in km/s (i.e. in velocity space) to considered around the redshifted line centroid for the min and maximum allowed peculiar velocity values

    Returns
    -------
    line_model: lmfit.Model
        A lmfit model with all the lines to be fitted. The lines are defined by their name, wavelength, and reference line if the ratio has to be constrained.
        Note this model does not contain the continuum
    """
     # initialize model with the continuum first
    line_model = PolynomialModel(degree=degree, prefix=f"{cont_prefix}_")

    for linename in lines.keys():
        line = lines[linename]
        # assign parameters, redshift the line unless it's a sky line
        velexpr = f"((1 + {linename}_vel/ {ckms:.3f} ) / (1 - {linename}_vel/ {ckms:.3f})) **0.5 * (1 + {redshift:.8f}) * {line.wave:.4f}" if not "sky" in linename else f"((1 + {linename}_vel/ {ckms:.3f}) / (1 - {linename}_vel/ {ckms:.3f})) **0.5 * {line.wave:.4f}"
        redshifted = (1 + redshift) * line.wave if not "sky" in linename else line.wave
            
        line_model += GaussianModel(prefix="%s_" % linename)

        # add a broadened version of the same line
        if "broad" in linename:
            logger.info(f"Adding broad component for {linename}")
            line_model.set_param_hint("%s_vel" % linename, value=0., min=-vmargin, max=vmargin, vary=True)

            line_model.set_param_hint("%s_center" % linename, expr=velexpr)
            # add param to force the width of the broad line always larger than the narrow component
            # add the name of the line in case there are TWO broad lines
            line_model.set_param_hint("%s_fwhm_factor" % linename, value=2, min=1., max=7.5, vary=True)
            originallinename = linename.replace("_broad", "")
            # the FWHM is a factor broader than the narrow line
            line_model.set_param_hint("%s_fwhm_vel" % linename, expr=f"{linename}_fwhm_factor * {originallinename}_fwhm_vel")
            line_model.set_param_hint("%s_sigma" % linename, 
                                      expr=f"{linename}_fwhm_vel / 2.355 * {linename}_center / {ckms:.3f}", value=sigma * 2.0)
            # everything set, we are done here
            continue

        # new narrow line
        else:
            if line.ref is None:
                ##line_model.set_param_hint("%s_%s_tied" % (linename, previousline.name), min=0, max=100, value=line.wave - previousline.wave)
                if not "sky" in linename:
                    line_model.set_param_hint("%s_vel" % linename, value=0., min=-vmargin, max=vmargin, vary=True)
                else:
                    # sky lines should be centered on the ref line, so allow less margin for velocity shifts
                    line_model.set_param_hint("%s_vel" % linename, value=0., min=-100, max=100, vary=True)
                line_model.set_param_hint("%s_center" % linename, expr=velexpr,)

                FWHMinst = get_instrument_FWHM(redshifted)[0] / redshifted * ckms  # convert to velocity space
                minFWHM = FWHMinst * 0.7
                maxFWHM = FWHMinst * 4 if not "sky" in linename else FWHMinst * 2
                FWHMvalue = sigma * sigma_to_fwhm / redshifted * ckms  # convert to sigma
                line_model.set_param_hint("%s_fwhm_vel" % linename, value=FWHMvalue, min=minFWHM, max=maxFWHM)
                line_model.set_param_hint("%s_sigma" % linename, expr=f"{linename}_fwhm_vel / 2.355 * {linename}_center / {ckms:.3f}",
                                          value=sigma)
            else:
    
                # tight lines center and sigma to the reference lines
                refline = lines[line.ref]
                # link velocity
                line_model.set_param_hint("%s_vel" % linename, expr="%s_vel" % (refline.name))
                line_model.set_param_hint("%s_center" % linename, expr=velexpr)
                # link sigma to the other line, and infer sigma for this line based on fwhm_vel
                line_model.set_param_hint("%s_fwhm_vel" % linename, expr="%s_fwhm_vel" % refline.name)
                line_model.set_param_hint("%s_sigma" % linename, 
                                          expr="%s_fwhm_vel / 2.355 / %.3f * %s_center" % (linename, ckms, linename))

                # link the flux
                line_model.set_param_hint("%s_factor" % linename, value=line.th, 
                                            min=line.low, max=line.up)
                line_model.set_param_hint(name="%s_amplitude" % linename, expr='%s_amplitude * %s_factor' % (refline, linename))

    return line_model


def plot_fit(x, result, initial_model=None, cont_prefix="cont", normalize=True, annotate=True, lref=None, z_sys=0):
    xnew = np.linspace(x.min(), x.max(), 500)
    comps = result.eval_components(x=xnew)

    fig, axes = plt.subplots(2, 1,
                             height_ratios=(2, 1), sharex=True, gridspec_kw={'hspace': 0.00})
    ax = axes[0]
    if normalize:
        norm = result.data.max()
        ax.set_ylabel("Normalized Flux")
    else:
        norm = 1
        ax.set_ylabel(rf"Flux ($10^{-20}$ erg cm$^{-2}$ s${-1}$ $\AA^{-1}$)")

    if lref:
        axes[-1].set_xlabel(fr'Velocity (km/s)')
        xplotnew = peculiar_velocity(lref, xnew, z_sys=z_sys)[0]
        xplot = peculiar_velocity(lref, x, z_sys=z_sys)[0]
    else:
        axes[-1].set_xlabel(r'Wavelength ($\AA$)')
        xplotnew = xnew
        xplot = x


    dely = result.eval_uncertainty(sigma=1, x=xnew)
    bestfit = result.eval(x=xnew)

    lines = ax.plot(xplotnew, bestfit / norm, label=r'Best Fit', color="C0")
    ax.fill_between(xplotnew, (bestfit - dely) / norm, 
                (bestfit + dely) / norm, color=lines[0].get_color(), alpha=0.2, zorder=-10)
    if annotate:
        ax.text(0.7, 0.9, rf"$\chi_r^2 = {result.redchi:.2f}$", fontsize=18, transform=ax.transAxes)
        ax.text(0.2, 0.9, rf"$R^2= {result.rsquared:.2f}$", fontsize=18, transform=ax.transAxes)
    #plt.text(0.2, 0.8, rf"SNR= {fit_data["params"]["%s_height" % (linename)].value / fit_data["noise"]:.2f}", fontsize=18, transform=ax.transAxes)
    ax.errorbar(xplot, result.data / norm, yerr=1 / result.weights / norm, label=r'Data', color="black")
    if initial_model is not None:
        ax.plot(xplot, get_initial_fit(initial_model, x) / norm, label=r'Initial Fit', ls="--")
    for i,comp in enumerate(comps.keys()):
        if cont_prefix not in comp:
            ax.plot(xplotnew, (comps[comp] + comps[f"{cont_prefix}_"]) / norm, ls="--", zorder=-10, color=f"C{i+1}")
    res = result.best_fit - result.data
    ax = axes[1]
    ax.errorbar(xplot, res  * result.weights, yerr=1, color="black")
    ax.axhline(0, ls="--", color="black")
    axes[-1].set_ylabel(fr'Residuals')
    axes[0].legend()
    return fig



def fit_spectrum(spectrum, fit_lines, redshift, sigma, wavelengths, degree=1, uncertainties=False):
    """Fit a spectrum to the lines defined in fit_lines. The lines are fitted with a Gaussian model and the continuum is fitted with a polynomial model.
    Parameters
    ----------
    spectrum: mpdaf.obj.Spectrum,

    redshift: float,
        Redshift of the object
    degree:int,
        Degree of the polynomial to fit the continuum
    """
    cont_prefix = 'cont'
    spectrumfluxes = spectrum.data
    spectrummask = spectrum.mask
    if np.all(spectrummask):
        raise SpectrumMaskedError("Spectrum is completely masked and cannot be fitted")
    else:
        usefulfluxes = spectrumfluxes[~spectrummask].data
        usefulwavelengths = wavelengths[~spectrummask]
        usefulstd = spectrum.var[~spectrummask] **0.5
    # create lmfit model
        line_model = create_line_model(fit_lines, redshift, sigma=sigma, degree=degree, cont_prefix=cont_prefix)
        nparams = len(line_model.param_names)
    # check if we have enough fluxes to fit the model
        if nparams > len(usefulfluxes):
            raise DofError(f"Insufficient degrees of freedom: {len(usefulfluxes)} data points but {nparams} parameters")
   
    _, median, stddev = sigma_clipped_stats(usefulfluxes, sigma=3)
    # min should be 0, but data is noisy so let's set it to the min of the data. Add 0.1 to avoid same min max
    # Convert to standard Python float for consistency
    c1 = (usefulfluxes[-1] - usefulfluxes[0]) / (usefulwavelengths[-1] - usefulwavelengths[0])
    c0 = usefulfluxes[0] - c1 * usefulwavelengths[0]
    line_model.set_param_hint("%s_c0" % cont_prefix, value=c0, min=-np.abs(c0) * 4., max = np.abs(c0) * 4.)#) min = usefulfluxes.min(), max=usefulfluxes.max())# min=usefulfluxes.min() - 5. *stddev, max=median + stddev * 5.0)
    if degree>=1:
        line_model.set_param_hint("%s_c1" % cont_prefix, value=c1, min = -np.abs(c1) * 2., max = np.abs(c1) * 3.)#min=-10000., max=10000.)
    for linename in fit_lines.keys():
        # assign sigma and amplitude base on data
        refline = fit_lines[linename].ref
        # no reference line
        if refline is None:
            centroid = compute_shift(fit_lines[linename].wave, redshift, line_model.param_hints[f"{linename}_vel"]["value"])[0]
            central_wavelength = np.argmin(np.abs(usefulwavelengths - centroid))
            min_ind = central_wavelength - 4 if central_wavelength - 4 > 0 else 0
            max_ind = central_wavelength + 4 if central_wavelength + 4 < len(usefulwavelengths) else len(usefulwavelengths) - 1
            peak_index = np.argmax(usefulfluxes[min_ind:max_ind]) + min_ind
            sigma = line_model.param_hints[f"{linename}_sigma"]["value"]
            init_amplitude = (usefulfluxes[peak_index] - median) * sqrt2pi * sigma 
            init_amplitude = 0. if init_amplitude < 0 else init_amplitude
            # don't set min to 0 in case init_amplitude is 0
            maxheight = np.max(usefulfluxes) - (median - stddev) + 5 * usefulstd.max()
            # set 10 times the initial amplitude, but if too large restrict it to the max height of the data  + the 5 * std
            maxamplitude = init_amplitude * 10 if init_amplitude * 10 < maxheight * sqrt2pi * (2. * sigma) else maxheight * sqrt2pi * (2. *sigma)
            # Convert all values to standard Python float for consistency
            line_model.set_param_hint("%s_amplitude" % linename, value=init_amplitude, min=-0.1, max=maxamplitude, vary=True)
        # reference line already linked

    logger.debug(f"Fitting spectrum with {len(usefulfluxes)} datapoints with {nparams} parameters")
    logger.debug("Initial model")
    logger.debug(line_model.param_hints)
    result = line_model.fit(usefulfluxes, x=usefulwavelengths, weights=1.0 / usefulstd, nan_policy='raise')
    
    conf = None
    if result.errorbars and result.rsquared > 0.0 and uncertainties:
        try:
            conf = result.conf_interval(sigmas=[1])
        except ValueError as e:
            return result, -1

    #if result.redchi>100:
      # 

    return result, conf

