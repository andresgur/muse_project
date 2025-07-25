# -*- coding: utf-8 -*
import numpy as np
import argparse, os, time, logging
from line import Line, Lines
from mpdaf.obj import Cube, iter_spe
from tqdm.contrib.concurrent import process_map
from functools import partial
import warnings
from astropy.stats import sigma_clipped_stats
from lmfit.models import PolynomialModel, GaussianModel
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.utils.exceptions import AstropyWarning
from astropy.constants import c


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('fit_cube')
warnings.filterwarnings('ignore', category=AstropyWarning, module='astropy')
warnings.filterwarnings('ignore', message='.*UFloat objects with std_dev==0.*')
CATALOG_LINES = Lines()
ckms = 299792.458


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
    efactor = 2 * (1 + zp) * ezp
    vp = ckms * (factor - 1) / (factor + 1)
    #df/factor = 1 x (factor + 1) - (factor -1)X1 / (factor + 1)^2 ==> 2 / (Factor + 1)^2
    evp = ckms * 2 * efactor / (factor + 1)**2
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
    return A * wavelength ** 2 - B * wavelength + 6.04, A * 2 * wavelength * ewavelength - B * ewavelength

def get_initial_fit(model, x):
    params = model.make_params()
    return model.eval(params, x=x)

def create_line_model(lines, redshift=0.0, sigma=1.4, margin=10):
    """Create lmfit model with all the lines to be fitted. Lines sigmas and position are tight to the reference values of other lines.
    Parameters
    ---------
    lines: list or array of Line objects,
        Lines to be fitted
    redshift: float
        Initial guess for the redshift
    sigma: float
        Initial guess for the line sigma in Angstroms
    margin: float
        +-Margin in Angstroms to considered around the redshifted line centroid for the min and maximum allowed values

    Returns
    -------
    line_model: lmfit.Model
        A lmfit model with all the lines to be fitted. The lines are defined by their name, wavelength, and reference line if the ratio has to be constrained.
        Note this model does not contain the continuum
    """
    line_model = None
     # initialize model with first line
    for line in lines:
        linename = line.name
        # create the model for the first line
        if line_model is None:
            line_model = GaussianModel(prefix="%s_" % linename)
        else:
            line_model += GaussianModel(prefix="%s_" % linename)
    # assign parameters
        redshifted = (1 + redshift) * line.wave

        line_model.set_param_hint("%s_center" % line.name, value=redshifted, vary=True, max=redshifted + margin,
                                      min=redshifted - margin)
        if line.ref is None:
            line_model.set_param_hint("%s_sigma" % line.name, value=sigma, vary=True, max=4, min=0.8)
        # tight lines center and sigma to the reference lines
        else:
            refline = CATALOG_LINES.lines[line.ref]
            line_model.set_param_hint("%s_center" % line.name, expr="%s_center + %.1f" % (refline.name, (line.wave - refline.wave)))
            line_model.set_param_hint("%s_sigma" % line.name, expr="%s_sigma" % refline.name)

    return line_model


def create_out_maps(lines, shape):
    """Create output maps
    
    lines: list of Line objects
        List of lines that have been fitted to the data
    shape: tuple
        Shape of the output maps, usually the 2D shape of the cube (y, x)

    Returns
    -------
    outmaps: dict
        Dictionary with the output maps for each line and parameter. The keys are the line name and parameter name, e.g. "HBETA_amplitude", "OIII5007_fwhm
    """
    outmaps = { }
    outpars = ["amplitude", "fwhm", "center", "snr"]
    for line in lines:
        for par in outpars:
            outmaps["%s_%s" % (line, par)] = np.full(shape, np.nan, dtype=float)
            # uncertainty map
            outmaps["%s_e%s" % (line, par)] = np.full(shape, np.nan, dtype=float)

    outmaps["redchi"] = np.full(shape, np.nan, dtype=float)
    outmaps["residuals"] = np.full(shape, np.nan, dtype=float)
    return outmaps


def fit_spectrum(spectrum_pixels, fit_lines, redshift, sigma, wavelengths, wav_cuts, linegroups, degree=1):
    """Fit a spectrum to the lines defined in fit_lines. The lines are fitted with a Gaussian model and the continuum is fitted with a polynomial model.
    Parameters
    ----------
    redshift: float,
        Redshift of the object
    degree:int,
        Degree of the polynomial to fit the continuum
    """
    spectrum, (y, x) = spectrum_pixels
    spectrumfluxes = spectrum.data
    if spectrumfluxes.all() is np.ma.masked:
        logger.debug("Spectrum at X:%d Y%d completely masked, skipping" % (x, y))
        return None
    else:
        usefulfluxes = spectrumfluxes[~spectrumfluxes.mask].data
        usefulwavelengths = wavelengths[~spectrumfluxes.mask]
        if spectrum.var is not None:
            usefulvar = spectrum.var[~spectrumfluxes.mask]
        else:
            usefulvar = np.full(usefulfluxes.shape, np.var(usefulfluxes), dtype=float)
    # create lmfit model
        line_model = create_line_model(fit_lines, redshift, sigma=args.sigma)
        cont_prefix = 'cont'
        cont_model = PolynomialModel(degree=degree, prefix="%s_" %cont_prefix)
        line_model += cont_model     
        nparams = len(line_model.param_names)
    # check if we have enough fluxes to fit the model
        if nparams > len(usefulfluxes):
            logger.debug("Spectrum at X:%d Y%d has only %d fluxes but %d parameters, skipping" % (x, y, len(usefulfluxes), nparams))
            return None
   
    _, median, _ = sigma_clipped_stats(usefulfluxes, sigma=3)
    # min should be 0, but data is noisy so let's set it to the min of the data. Add 0.1 to avoid same min max
    line_model.set_param_hint("%s_c0" % cont_prefix, value=median, vary=True, min=np.min(usefulfluxes), max=np.max(usefulfluxes) + 0.1)
    line_model.set_param_hint("%s_c1" % cont_prefix, value=0, vary=True)
    
    for group in linegroups:
        # mask the fluxes and return the new fluxes with more masked values
        # mask outside the wavelength range of interest
        # these fluxes are just to compute the initial values of the lines
        for linename in group:
            # assign sigma and amplitude base on data
            refline = CATALOG_LINES.lines[linename].ref
            if refline is None:
                centroid = line_model.param_hints[f"{linename}_center"]["value"]
                central_wavelength = np.argmin(np.abs(usefulwavelengths - centroid))
                min_ind = central_wavelength - 4 if central_wavelength - 4 > 0 else 0
                max_ind = central_wavelength + 4 if central_wavelength + 4 < len(usefulwavelengths) else len(usefulwavelengths) - 1
                peak_index = np.argmax(usefulfluxes[min_ind:max_ind]) + min_ind

                sigma = line_model.param_hints[f"{linename}_sigma"]["value"]
                init_amplitude = (usefulfluxes[peak_index] - median) * np.sqrt(2 * np.pi * (sigma ** 2))
                line_model.set_param_hint("%s_amplitude" % linename, value=init_amplitude, vary=True, min=0, max=init_amplitude * 50 + 1)
            else:
                line_model.set_param_hint("%s_factor" % linename, value=CATALOG_LINES.lines[linename].th, 
                                          min=CATALOG_LINES.lines[linename].low, max=CATALOG_LINES.lines[linename].up)
                line_model.set_param_hint(name="%s_amplitude" % linename, expr='%s_amplitude * %s_factor' % (refline, linename))
    # mask the range we want to use
    # do this in a two step process
    #fluxes = ma.masked_where( (usefulwavelengths < wav_cuts[0][0]) | (usefulwavelengths > wav_cuts[0][1]), usefulfluxes)

    # the for loop won't be entered if wav_cuts has only one element
    #for wav_range in wav_cuts[1:]:
     #   fluxes = ma.masked_where((usefulwavelengths < wav_range[0]) | (usefulwavelengths > wav_range[1]), fluxes)
    
    #if fluxes.all() is np.ma.masked:
     #   logger.debug("Spectrum at X:%d Y%d completely masked after filtering, skipping" % (x, y))
      #  return None
    #if nparams > len(fluxes[~fluxes.mask]):
     #   logger.debug("Spectrum at X:%d Y%d has only %d fluxes but %d parameters, skipping" % (x, y, len(fluxes[~fluxes.mask]), nparams))
      #  return None
        
    logger.debug(f"Full model")
    logger.debug(line_model.param_hints)
    std_dev = np.sqrt(usefulvar)
    logger.debug(f"Fitting spectrum at X:{x} Y:{y} with {len(usefulfluxes)} datapoints with {len(line_model.param_hints)} parameters")
    result = line_model.fit(usefulfluxes, x=usefulwavelengths, weights=1.0 / std_dev, nan_policy='raise')

    fit_data = {
        'x': x,
        'y': y,
        'ndata': result.ndata,
        'nvarys': result.nvarys,
        'redchi': result.redchi,
        'residual': result.residual.copy(),
        'noise': np.std(usefulfluxes - result.best_fit),
        'best_fit': result.best_fit.copy(),
        'errorbars': result.errorbars if hasattr(result, 'errorbars') else False,
        #'covar': result.covar.copy() if result.covar is not None else None,
        'params': result.params
    }

    #plt.errorbar(usefulwavelengths, usefulfluxes, yerr=std_dev, label='Data')
    #plt.plot(usefulwavelengths, get_initial_fit(line_model, usefulwavelengths), label=r'Initial Fit', ls="--")
    #plt.plot(usefulwavelengths, result.best_fit, label='Best Fit')
    #plt.show()
    return fit_data

def check_input_lines(linegroups):
    """Check if the input lines are valid and exist in the catalog
    
    """
    for group in linegroups:
        for line in group.split(","):
            if line not in CATALOG_LINES.lines:
                raise ValueError(f"Line {line} not found in catalog. Please check the line name or add it to the catalog. Available lines: {CATALOG_LINES.lines.keys()}")
    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Fit data cube using voronoi binning')
    parser.add_argument("-z", "--redshift", nargs='?', help="Initial guess for the redshift", default=0.0, type=float)
    parser.add_argument("-s", "--sigma", nargs='?', help="Initial guess for the width (sigma) of the lines in Angstroms. Default 1.4 Angstroms", default=1.4, type=float)
    parser.add_argument("-o", "--outdir", nargs='?', help="Name of the output directory", default="voronoi", type=str)
    parser.add_argument("input_cube", help="Path to input cube to be analysed", type=str)
    parser.add_argument("-w", "--window", help="+- Angstroms to cut the spectrum around each line centroid. Default 25 Angstroms. Reduce for lines close to telluric lines (e.g. [OI]6300)", type=float, default=25)
    parser.add_argument("-c", "--cpus", nargs='?', help="Number of CPUs to use for parallelization. By default it uses all -1 available cpus", type=int)
    parser.add_argument("-l", "--linegroups", nargs='+', default="HeII4686", help="Space-separated groups of comma-separated lines to fit. (e.g. if you want two groups with same continuum: HBETA OIII4959,OIII5007)", type=str) # 
    args = parser.parse_args()

    outpath = args.outdir
    sigma = args.sigma
    z = args.redshift

    if not os.path.isdir(outpath):
        os.mkdir(outpath)


    if args.cpus is None:
        cpus = os.cpu_count() - 1
    else:
        cpus = args.cpus

    check_input_lines(args.linegroups)

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
    fit_lines = []
    linegroups = []

    for index, group in enumerate(inputgroups):
        grouplines = group.split(",")
        linegroups.append(grouplines)
        # make wavelength range based on the first and last line of the group
        minwav = CATALOG_LINES.lines[grouplines[0]].wave * (1 + redshift) - margin
        maxwav = CATALOG_LINES.lines[grouplines[-1]].wave * (1 + redshift) + margin
        wav_cuts[index] = (minwav, maxwav)
        for line in grouplines:
            fit_lines.append(CATALOG_LINES.lines[line])
            
    
    # get the first line of the first group and the last line of the last group
    waveinit, waveend = fit_lines[0].wave * (1 + redshift) - margin, fit_lines[-1].wave * (1 + redshift) + margin
    cubename = args.input_cube
    data_cube = Cube(cubename)
    step = data_cube.wave.get_step()
    # cut the cube over the wavelength range we will not needed
    data_cube = data_cube.select_lambda(waveinit - step *1.01, waveend + step * 1.01)

    # mask the data cube outside the wavelength range of interest, important to copy here
    original_mask = data_cube.data.mask.copy()  
    # mask the entire range by default
    data_cube[:,:,:] = ma.masked
    # in the wavelength range of interest, keep the original mask
    for index, wav_range in enumerate(wav_cuts):
        lmin, lmax = wav_range
        print(f"Masking wavelength range {lmin} - {lmax} Angstroms")
        minindex = data_cube.wave.pixel(lmin, nearest=True)
        maxindex = data_cube.wave.pixel(lmax, nearest=True)
        # Unmask this wavelength range, but preserve original mask
        data_cube.data.mask[minindex:maxindex, :, :] = original_mask[minindex:maxindex, :, :]

    if data_cube.data.all() is np.ma.masked:
        logger.error("All data in the cube is masked. Please check the wavelength range and the input cube.")
        exit(1)
    wavelengths = data_cube.wave.coord()

    logger.info(f"Fitting cube {cubename} with {len(fit_lines)} lines,  {len(linegroups)} line groups, and {cpus} CPUs")
    total_spectra = data_cube.shape[1] * data_cube.shape[2]  # y * x dimensions
    logger.info(f"Total spectra to fit: {total_spectra}")
    start = time.time()
    if cpus>1:
        results = process_map(partial(fit_spectrum, fit_lines=fit_lines, redshift=z, sigma=sigma, wavelengths=wavelengths, wav_cuts=wav_cuts, linegroups=linegroups), 
                iter_spe(data_cube, index=True), max_workers=cpus, chunksize=1,
                unit=" spectra", desc="Fitting spectra", total=total_spectra)
    else:
        results = []
        for spectrum in iter_spe(data_cube, index=True):
            result = fit_spectrum(spectrum, fit_lines, redshift=z, sigma=sigma, wavelengths=wavelengths, wav_cuts=wav_cuts, linegroups=linegroups)
            results.append(result)
    
    end = time.time()
    logger.info(f"Fitting completed in {end - start:.2f} seconds in {cpus} CPUs")
    logger.info("Storing results in output maps")
    # the 2D coordinates shape of the cube (index 0 is wavelength axis)
    outshape = data_cube.shape[1:]
    outmaps = create_out_maps(fit_lines, outshape)
    # cube for residuals
    #data_cube.crop()
    #residual_cube = data_cube.clone()
    #residual_cube.data = np.zeros(data_cube.shape)

    redchis = []
    ndofs = []
    ndata = []
    xs = []
    ys = []
    
    for fit_data in results:

        if fit_data is None:
            continue
        else:
            x = fit_data['x']
            y = fit_data['y']
            
            # Store reduced chi square
            outmaps["redchi"][y, x] = fit_data['redchi']

            xs.append(x)
            ys.append(y)
            redchis.append(fit_data['redchi'])
            ndofs.append(fit_data['ndata'] - fit_data['nvarys'])    
            ndata.append(fit_data['ndata'])

            best_values = fit_data['params']
            for line in fit_lines:
                # Calculate SNR
                signal = best_values["%s_height" % (line.name)].value
                noise = fit_data["noise"]
                outmaps["%s_snr" % (line.name)][y, x] = signal / noise if noise > 0 else np.nan
                
                # Store parameter values
                for par in ["amplitude", "fwhm", "center"]:
                    param_name = "%s_%s" % (line.name, par)
                    if param_name in best_values:
                        param = best_values[param_name]
                        outmaps[param_name][y, x] = param.value

                        if param.stderr is not None:
                            outmaps["%s_e%s" % (line.name, par)][y, x] = param.stderr
            #warnings.warn("Warning: noise computation might be incorrect: see meaning of residuals here https://lmfit.github.io/lmfit-py/model.html#modelresult-attributes")
            #residual_cube[:, y, x] = fit_data['residual'] --> we need to fix not all spectra having same number of pixels


    linepars = ["amplitude", "fwhm", "center"]
    # dummy skeleton for the output cube
    outfile = data_cube[0 ,: , :]
    for line in fit_lines:
        # FWHM corrected map
        FWHM_inst, eFWHM_inst = get_instrument_FWHM(outmaps["%s_center" % (line.name)], outmaps["%s_ecenter" % (line.name)])
        # FWHM
        selector = FWHM_inst > outmaps["%s_fwhm" % (line.name)]
        FWHM_corrected = np.sqrt(outmaps["%s_fwhm" % (line.name)] ** 2 - FWHM_inst ** 2)
        FWHM_corrected[selector] = 0 # these are upper limits
        outfile.data = FWHM_corrected
        par = "fwhm_corrected"
        outfile.write("%s/%s_%s.fits" % (outpath, line, par))
        #eFWHM f = sqrt(A^2-B^2); df/dA = 2A/2f = A/f; df/dB = -2B/2f = -B/f; eFWHM = sqrt((eA/f)^2 + (eB/f)^2) = sqrt(eA^2 + eB^2)/f

        Ae = outmaps["%s_efwhm" % (line.name)] * outmaps["%s_fwhm" % (line.name)]
        Be = FWHM_inst * eFWHM_inst
        eFWHM_corrected = np.sqrt(Ae ** 2 + Be ** 2) / FWHM_corrected
        eFWHM_corrected[selector] = np.nan # these are upper limits
        outfile.data = eFWHM_corrected
        par = "efwhm_corrected"
        outfile.write("%s/%s_%s.fits" % (outpath, line, par))

        # velocity map
        velocity, evelocity = peculiar_velocity(line.wave, outmaps["%s_center" % (line.name)], 
                                    ewavelength=outmaps["%s_ecenter" % (line.name)], z_sys=redshift)
        par = "vel"
        outfile.data = velocity
        outfile.write("%s/%s_%s.fits" % (outpath, line, par))

        par = "evel"
        outfile.data = evelocity
        outfile.write("%s/%s_%s.fits" % (outpath, line, par))

        for par in linepars:
            # parameter value
            outfile.data = outmaps["%s_%s" % (line, par)]
            outfile.write("%s/%s_%s.fits" % (outpath, line, par))
            # parameter uncertainty
            outfile.data = outmaps["%s_e%s" % (line, par)]
            outfile.write("%s/%s_e%s.fits" % (outpath, line, par))
        
        outfile.data = outmaps["%s_snr" % (line)]
        outfile.write("%s/%s_snr.fits" % (outpath, line))
        
    outfile.data = outmaps["redchi"]
    linegroupout = "_".join([line for group in linegroups for line in group])
    outfile.write(f"{outpath}/redchi_{linegroupout}.fits")

    # residual_cube.write("%s/residual_cube.fits" % outpath)
    np.savetxt(f"{outpath}/fit_results_{linegroupout}.txt",
               np.array([xs, ys, redchis, ndofs, ndata]).T,
               header="x\ty\tredchi\tndof\tndata",
               fmt="%d\t%d\t%.4f\t%d\t%d")
    print("Output stored in %s" % outpath)