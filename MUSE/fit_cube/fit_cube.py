#!/usr/bin/env python
# -*- coding: utf-8 -*
import numpy as np
import argparse, os, time, logging
from line import Lines
from lineutils import get_instrument_FWHM, compute_shift, correct_FWHM
from fitutils import fit_spectrum, SpectrumMaskedError, DofError, get_error
from mpdaf.obj import Cube, Image, iter_spe
import warnings
import numpy.ma as ma
import matplotlib.pyplot as plt
from astropy.utils.exceptions import AstropyWarning


logger = logging.getLogger('fit_cube')

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

warnings.filterwarnings('ignore', category=AstropyWarning, module='astropy')
warnings.filterwarnings('ignore', category=UserWarning, module='lmfit')
warnings.filterwarnings('ignore', message='.*UFloat objects with std_dev==0.*')
logging.getLogger('matplotlib').setLevel(logging.WARNING)
CATALOG_LINES = Lines()
ckms = 299792.458
sigma_to_fwhm = 2.355  # FWHM = 2.355 * sigma
sqrt2pi = (2. * np.pi) ** 0.5

def save_image(data, templatefile, dtype=float, degree=1, redshift=0.0, margin=25, filename="output.fits", bunit=""):
    outfile = Image(data=data, wcs=templatefile.wcs, dtype=dtype, unit=bunit)
    outfile.primary_header = templatefile.primary_header
    outfile.primary_header["DEGREE"] = degree
    outfile.primary_header["REDSHIFT"] = redshift
    outfile.primary_header["WINDOW"] = margin
    outfile.data_header["OBJECT"] = templatefile.data_header.get("OBJECT", "Unknown")
    outfile.write(filename)

def create_out_maps(lines, shape, degree=1):
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
    outpars = ["amplitude", "vel", "fwhm_vel", "fwhm_vel_corrected", "snr"]
    for line in lines:
        for par in outpars:
            outmaps["%s_%s" % (line, par)] = np.full(shape, np.nan, dtype=float)
            # uncertainty map
            outmaps["%s_e%s" % (line, par)] = np.full(shape, np.nan, dtype=float)

    for i in range(degree + 1):
        outmaps["cont_c%d" % i] = np.full(shape, np.nan, dtype=float)
        outmaps["cont_ec%d" % i] = np.full(shape, np.nan, dtype=float)
    outmaps["redchi"] = np.full(shape, np.nan, dtype=float)
    outmaps["BIC"] = np.full(shape, np.nan, dtype=float)
    outmaps["rsquared"] = np.full(shape, np.nan, dtype=float)
    outmaps["residuals"] = np.full(shape, np.nan, dtype=float)
    return outmaps


def fit_spectrum_par(spectrum, fit_lines, redshift, sigma, wavelengths, degree=1, out="outfile", uncertainties=False):
    """Fit a spectrum to the lines defined in fit_lines. The lines are fitted with a Gaussian model and the continuum is fitted with a polynomial model.
    Parameters
    ----------
    redshift: float,
        Redshift of the object
    degree:int,
        Degree of the polynomial to fit the continuum
    """
    try:
        result, conf = fit_spectrum(spectrum, fit_lines, redshift, sigma, wavelengths, degree, uncertainties)
    except SpectrumMaskedError:
        return None
    except DofError as dof:
        return None
    #fig = plot_fit(spectrum.wave.coord()[~spectrum.data.mask], result)
    #plt.savefig(f"{out}_fit_spectrum_X{x}_Y{y}.png")
    #plt.close(fig)
        #outstr = f"Conf interval problem with for X:{x} Y:{y}\n"
        #outstr += f"Fitting info:\n"
        #outstr += result.fit_report()
        #logger.info(outstr)
        

    fit_data = {
        'ndata': result.ndata,
        'nvarys': result.nvarys,
        'redchi': result.redchi,
        'BIC': result.bic,
        'residual': result.residual,
        'noise': np.std(result.data - result.best_fit),
        'best_fit': result.best_fit,
        'errorbars': result.errorbars,
        'rsquared': result.rsquared,
        'conf': conf,
        'ier': result.ier,
        #'covar': result.covar.copy() if result.covar is not None else None,
        'params': result.params,
        'report': result.fit_report()
    }


    #print(f"Fitting result for X:{x} Y:{y}")
    #print(result.fit_report())
    #if conf is not None:
     #   print(report_ci(conf))
    #fig = plot_fit(usefulwavelengths, usefulfluxes, result, std_dev, line_model)
    #plt.savefig(f"{outpath}/fit_spectrum_X{x}_Y{y}.png")
    #plt.close(fig)
    
    return fit_data

def check_input_lines(linegroups):
    """Check if the input lines are valid and exist in the catalog, and are in ascending wavelength order
    
    """
    for group in linegroups:
        lines_in_group = group.split(",")
        
        # Check if lines exist in catalog
        for line in lines_in_group:
            if line not in CATALOG_LINES.lines:
                raise ValueError(f"Line {line} not found in catalog. Please check the line name or add it to the catalog. Available lines: {CATALOG_LINES.lines.keys()}")
        
        # Check if lines are in ascending wavelength order
        if len(lines_in_group) > 1:
            wavelengths = [CATALOG_LINES.lines[line].wave for line in lines_in_group]
            if wavelengths != sorted(wavelengths):
                raise ValueError(f"Lines in group '{group}' are not in ascending wavelength order. "
                               f"Current order: {', '.join([f'{line}({CATALOG_LINES.lines[line].wave:.1f}Å)' for line in lines_in_group])}. "
                               f"Please reorder them by wavelength.")
    
    return True


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Fit data cube using voronoi binning')
    parser.add_argument("-d", "--degree", nargs='?', help="Polynomial degree for the continuum", default=1, type=int)
    parser.add_argument("-z", "--redshift", nargs='?', help="Initial guess for the redshift", default=0.0, type=float)
    parser.add_argument("-s", "--sigma", nargs='?', help="Initial guess for the width (sigma) of the lines in Angstroms. Default 1.4 Angstroms", default=1.4, type=float)
    parser.add_argument("-o", "--outdir", nargs='?', help="Name of the output directory", default="fit_cube", type=str)
    parser.add_argument("input_cube", help="Path to input cube to be analysed", type=str)
    parser.add_argument("-w", "--window", help="+- Angstroms to cut the spectrum around each line centroid. Default 25 Angstroms. Reduce for lines close to telluric lines (e.g. [OI]6300)", type=float, default=25)
    parser.add_argument("-c", "--cpus", nargs='?', help="Number of CPUs to use for parallelization. By default it uses all -1 available cpus", type=int)
    parser.add_argument("-l", "--linegroups", nargs='+', default="HeII4686", help="Space-separated groups of comma-separated lines to fit. (e.g. if you want two groups with same continuum: HBETA OIII4959,OIII5007)", type=str) # 
    parser.add_argument("-u", "--uncertainties", action="store_true", help="Compute confidence intervals for fit parameters. This will slow down the fitting process.")
    args = parser.parse_args()

    outpath = args.outdir
    sigma = args.sigma
    degree = args.degree

    if not os.path.isdir(outpath):
        os.mkdir(outpath)


    if args.cpus is None:
        cpus = os.cpu_count() - 1
    else:
        cpus = args.cpus

    check_input_lines(args.linegroups)

    # Create log file name based on line groups
    linegroupout = "_".join([line for group in args.linegroups for line in group.split(",")])
    logfile = f"{outpath}/fit_cube_{linegroupout}.log"
    setup_logging(logfile, logging.INFO)
    logger.info(f"Started logging to {logfile}")

    inputgroups = args.linegroups#[["HeII4686"]]#[["HBETA", "OIII4959", "OIII5007"], ["OI6300"],
            #["NII6548","HALPHA", "NII6583"], ["SII6716", "SII6731"]]
            # ["NII5755"], ["HeI5875.64"], 
    logger.info("Line groups: %s",inputgroups)

    logger.info(f"Degree of continuum polynomial: {degree}")

    redshift = args.redshift
    # Margins to cut the spectrum around the blue and red lines of each group
    margin = args.window # Angstroms

    # create waveleneght cuts
    wav_cuts = np.zeros((len(inputgroups), 2), dtype=float)
    # we need these variables for later
    fit_lines = {}
    linegroups = []

    for index, group in enumerate(inputgroups):
        grouplines = group.split(",")
        linegroups.append(grouplines)
        # make wavelength range based on the first and last line of the group
        
        minwav = CATALOG_LINES.lines[grouplines[0]].wave * (1 + redshift) - margin
        maxwav = CATALOG_LINES.lines[grouplines[-1]].wave * (1 + redshift) + margin
        wav_cuts[index] = (minwav, maxwav)
        for line in grouplines:
            fit_lines[line] = CATALOG_LINES.lines[line]


    # get the first line of the first group and the last line of the last group
    waveinit, waveend = list(fit_lines.values())[0].wave * (1 + redshift) - margin, list(fit_lines.values())[-1].wave * (1 + redshift) + margin
    cubename = args.input_cube
    data_cube = Cube(cubename)
    step = data_cube.wave.get_step()
    # cut the cube over the wavelength range we will not needed
    data_cube = data_cube.select_lambda(waveinit - step *1.01, waveend + step * 1.01)

    # mask the data cube outside the wavelength range of interest, important to copy here
    original_mask = data_cube.data.mask.copy()  
    # mask the entire range by default
    data_cube.data[:,:,:] = ma.masked
    # in the wavelength range of interest, keep the original mask
    for index, wav_range in enumerate(wav_cuts):
        lmin, lmax = wav_range
        logger.info(f"Masking wavelength range {lmin} - {lmax} Angstroms")
        minindex = data_cube.wave.pixel(lmin, nearest=True)
        maxindex = data_cube.wave.pixel(lmax, nearest=True)
        # Unmask this wavelength range, but preserve original mask
        data_cube.mask[minindex:maxindex, :, :] = original_mask[minindex:maxindex, :, :]

    if np.all(data_cube.mask):
        logger.error("All data in the cube is masked. Please check the wavelength range and the input cube.")
        exit(1)
    wavelengths = data_cube.wave.coord()

    logger.info(f"Fitting cube {cubename} with {len(fit_lines)} lines,  {len(linegroups)} line groups, and {cpus} CPUs")
    total_spectra = data_cube.shape[1] * data_cube.shape[2]  # y * x dimensions
    logger.info(f"Total spectra to fit: {total_spectra}")
    start = time.time()
    if cpus>1:
        results = data_cube.loop_spe_multiprocessing(f=fit_spectrum_par, cpu=cpus, verbose=True, fit_lines=fit_lines, redshift=redshift, sigma=sigma, 
                                      wavelengths=wavelengths, degree=degree, 
                                      out=f"{outpath}/{linegroupout}", 
                                      uncertainties=args.uncertainties,)
    else:
        results = []
        for spectrum in iter_spe(data_cube, index=True):
            result = fit_spectrum_par(spectrum, fit_lines, redshift=redshift, sigma=sigma, wavelengths=wavelengths, 
                                      degree=degree, out=f"{outpath}/{linegroupout}", uncertainties=args.uncertainties)
            results.append(result)
    
    end = time.time()
    timetaken = end - start
    logger.info(f"Fitting completed in {timetaken:.2f} seconds in {cpus} CPUs")
    logger.info("Storing results in output maps")
    # the 2D coordinates shape of the cube (index 0 is wavelength axis)
    outshape = data_cube.shape[1:]
    outmaps = create_out_maps(fit_lines.keys(), outshape, degree=degree)
    # cube for residuals
    #data_cube.crop()
    #residual_cube = data_cube.clone()
    #residual_cube.data = np.zeros(data_cube.shape)

    redchis = []
    BICs = []
    rsquareds = []
    iers = []
    ndofs = []
    ndata = []
    xs = []
    ys = []
    cont_prefix = 'cont'
    for x in range(data_cube.shape[2]):
        for y in range(data_cube.shape[1]):
            fit_data = results[y, x]

            if fit_data is None:
                continue

            #fig, ax = plt.subplots()
            #plt.plot(wavelengths, data_cube.data[:, y, x], label='Data', color='black')
            #xnew = np.linspace(wavelengths.min(), wavelengths.max(), len(fit_data['best_fit']))
            #plt.plot(xnew, fit_data['best_fit'], label='Best Fit', color='red')
            #ax.text(0.7, 0.9, rf"$\chi_r^2 = {fit_data['redchi']:.2f}$", fontsize=18, transform=ax.transAxes)
            #ax.text(0.2, 0.9, rf"$R^2= {fit_data['rsquared']:.2f}$", fontsize=18, transform=ax.transAxes)

            # Store reduced chi square
            outmaps["redchi"][y, x] = fit_data['redchi']
            outmaps["BIC"][y, x] = fit_data['BIC']
            outmaps["rsquared"][y, x] = fit_data['rsquared']

            xs.append(x)
            ys.append(y)
            redchis.append(fit_data['redchi'])
            BICs.append(fit_data['BIC'])
            rsquareds.append(fit_data['rsquared'])
            ndofs.append(fit_data['ndata'] - fit_data['nvarys'])    
            ndata.append(fit_data['ndata'])
            iers.append(fit_data['ier'])

            best_values = fit_data['params']
            conf = fit_data["conf"]
            for linename in fit_lines.keys():
                # Calculate SNR
                signal = best_values["%s_height" % (linename)].value
                noise = fit_data["noise"]
                outmaps["%s_snr" % (linename)][y, x] = signal / noise if noise > 0 else np.nan
                #ax.text(0.2, 0.6, rf"$SNR = {outmaps['%s_snr' % (linename)][y, x]:.1f}$", fontsize=18, 
                 #       transform=ax.transAxes)

                # Store parameter values
                for par in ["amplitude", "vel", "fwhm_vel"]:
                    param_name = "%s_%s" % (linename, par)
                    if param_name in best_values:
                        param = best_values[param_name]
                        outmaps[param_name][y, x] = param.value
                        # Store parameter uncertainties
                        # get the uncertainty from the confidence intervals if available
                        if conf is not None:
                            # handle FWHM specially, as it is stored as sigma
                                                # handle the broad factor
                            if par == "fwhm_vel" and "broad" in linename:
                                confpar = conf["%s_fwhm_factor" %linename]
                                # multiply sigma of the narrow line with the factor
                                errorbroadfactor = get_error(confpar)
                                errorbroadfactor = errorbroadfactor if np.isfinite(errorbroadfactor) else param.stderr
                                fwhmvelnarrowline = best_values["%s_%s" % (linename.replace("_broad", ""), par)]
                                outmaps["%s_e%s" % (linename, par)][y, x] = errorbroadfactor * fwhmvelnarrowline
                            # if we have conf but some params are missing they are tied, so we just skip them if they are not in the conf interval
                            elif param_name in conf:
                                confpar = conf[param_name]
                                error = get_error(confpar)
                                outmaps["%s_e%s" % (linename, par)][y, x] = error if np.isfinite(error) else param.stderr
                            # keep the stderr for tied parameters if available
                            #elif param.stderr is not None:
                             #   outmaps["%s_e%s" % (line.name, par)][y, x] = param.stderr

                        # otherwise, use the standard error
                        elif param.stderr is not None:
                            outmaps["%s_e%s" % (linename, par)][y, x] = param.stderr

            #plt.savefig("%s/fit_spectrum_%s_X%d_Y%d.png" % (outpath, linegroupout, x, y))

            # Store continuum parameters
            for i in range(degree + 1):
                param = best_values[f"{cont_prefix}_c{i}"]
                outmaps["cont_c%d" % i][y, x] = param.value
                if conf is not None:
                        confpar = conf[f"{cont_prefix}_c{i}"]
                        outmaps["cont_ec%d" % i][y, x] = get_error(confpar)
                elif param.stderr is not None:
                    outmaps["cont_ec%d" % i][y, x] = param.stderr
            #warnings.warn("Warning: noise computation might be incorrect: see meaning of residuals here https://lmfit.github.io/lmfit-py/model.html#modelresult-attributes")
            #residual_cube[:, y, x] = fit_data['residual'] --> we need to fix not all spectra having same number of pixels


    linepars = ["amplitude", "fwhm_vel", "vel", "snr"]
    # dummy skeleton for the output cube
    templatefile = data_cube[0 ,: , :]
    for line in fit_lines.values():
        
        linename = line.name
        # store line param values
        for par in linepars:
            # parameter value
            save_image(data=outmaps["%s_%s" % (linename, par)], templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename="%s/%s_%s.fits" % (outpath, linename, par))
            # parameter uncertainty
            if par=="snr":
                # do not store the uncertainty for SNR
                continue
            save_image(data=outmaps["%s_e%s" % (linename, par)], templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename="%s/%s_e%s.fits" % (outpath, linename, par))

        # store wavelenght
        wavelength, ewavelength = compute_shift(line.wave, z_sys=redshift, velocity=outmaps["%s_%s" % (linename, "vel")], 
                                                                    evelocity=outmaps["%s_e%s" % (linename, "vel")])

        par = "center"
        save_image(data=wavelength, templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename=f"{outpath}/{linename}_{par}.fits", bunit="Angstrom")

        #ecenter
        save_image(data=ewavelength, templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename=f"{outpath}/{linename}_e{par}.fits", bunit="Angstrom")

        # FWHM
        FWHM_o = outmaps["%s_fwhm_vel" % (linename)]
        eFWHM_o = outmaps["%s_efwhm_vel" % (linename)]
        # Instrumental FWHM
        FWHM_inst, eFWHM_inst = get_instrument_FWHM(wavelength, ewavelength)
        FWHM_instvel = FWHM_inst / wavelength * ckms
        eFWHM_instvel =  ((eFWHM_inst / wavelength)**2. + (FWHM_inst * ewavelength / wavelength**2.))**0.5 * ckms

        # get pixels where FWHM_obs > FWHM_inst
        valid_selector = FWHM_o > FWHM_instvel
        
        par = "fwhm_vel_corrected"
        FWHM_corrected, eFWHM_corrected = correct_FWHM(FWHM_o[valid_selector], FWHM_instvel[valid_selector], eFWHM_o[valid_selector], eFWHM_instvel[valid_selector])
        outmaps[f"{linename}_{par}"][valid_selector] = FWHM_corrected
        outmaps[f"{linename}_{par}"][~valid_selector] = 0
        # store FWHM corrected maps
        save_image(data=outmaps[f"{linename}_{par}"], templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename="%s/%s_%s.fits" % (outpath, linename, par), bunit="km/s")

        outmaps[f"{linename}_e{par}"][valid_selector] = eFWHM_corrected # the non valid valus are simply nan
        save_image(data=outmaps[f"{linename}_e{par}"], templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename="%s/%s_e%s.fits" % (outpath, linename, par), bunit="km/s")

    for fitstat in ["redchi", "rsquared", "BIC"]:
        # write the reduced chi square and r squared maps
        save_image(data=outmaps[fitstat], templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename="%s/%s_%s.fits" % (outpath, fitstat, linegroupout))
        

    for i in range(degree + 1):
        bunit = templatefile.unit if i ==0 else f"{templatefile.unit}/Angstrom^{i}"
        save_image(outmaps["cont_c%d" % i], templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename="%s/cont_c%d_%s.fits" % (outpath, i, linegroupout), bunit=bunit)
        save_image(outmaps["cont_ec%d" % i], templatefile=templatefile, degree=degree, redshift=redshift, margin=margin,
                         filename="%s/cont_ec%d_%s.fits" % (outpath, i, linegroupout), bunit=bunit)


    # residual_cube.write("%s/residual_cube.fits" % outpath)
    np.savetxt(f"{outpath}/fit_results_{linegroupout}.txt",
               np.array([xs, ys, iers, redchis, rsquareds, ndofs, ndata]).T,
               header="x\ty\tier\tredchi\trsquared\tndof\tndata",
               fmt="%d\t%d\t%d\t%.3f\t%.3f\t%d\t%d")
    logger.info("Output stored in %s" % outpath)