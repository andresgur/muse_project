# @Author: Andrés Gúrpide <agurpide>
# @Date:   22-12-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 16-03-2022
import numpy as np
from voronoi import voronoi
import numpy as np
from mpdaf.obj import Cube, Spectrum
import astropy.units as u
import os
import time
from IPython.core.display import display, HTML
from lmfit.models import PolynomialModel, GaussianModel
from astropy.stats import sigma_clipped_stats
import numpy.ma as ma
import argparse


def fit_bin(bin_num, vor_output=None, data_cube=None, line_model=None, wav_ranges=None, els=None, line_groups=None, sigma=1.4, plot=False, outpath=None):
    """Performs fit to one Voronoi bin"""
    idx = vor_output[:,4] == bin_num
    x, y, x_mean, y_mean, _bn = vor_output[idx].T
    xint, yint= x.astype("int"), y.astype("int")
    xrav = np.ravel(np.asarray(xint))
    yrav = np.ravel(np.asarray(yint))
    bin_spectrum = get_weighted_spectrum(data_cube, xy[0], xy[1])
    mean, median, stddev = sigma_clipped_stats(bin_spectrum[:, 1], sigma=3)
    line_model.set_param_hint("%s_c0" % cont_prefix, value=median, vary=True, min=0)
    line_model.set_param_hint("%s_c1" % cont_prefix, value=0, vary=True)

    for index, wav_range in enumerate(wav_ranges):
        fluxes = ma.masked_where((bin_spectrum[:, 0] < wav_range[0]) & (bin_spectrum[:,0] > wav_range[1]), bin_spectrum[:,1])
        for linename in line_groups[index]:
            # assign sigma and amplitude base on data
            refline = els.lines[linename].ref
            if refline is None:
                centroid = els.lines[linename].wave * (1 + redshift)
                central_wavelength = np.argmin(np.abs(bin_spectrum[:,0] - centroid))
                min_ind = central_wavelength - 4
                max_ind = central_wavelength + 4
                peak_index = np.argmax(bin_spectrum[min_ind:max_ind,1]) + min_ind
                peak_wavelength = bin_spectrum[peak_index, 0]
                init_amplitude = (bin_spectrum[peak_index, 1] - median) * np.sqrt(2 * np.pi * (sigma ** 2))
                line_model.set_param_hint("%s_amplitude" % linename, value=init_amplitude, vary=True, min=0)
            else:
                line_model.set_param_hint("%s_factor" % linename, value=els.lines[linename].th, min=els.lines[linename].low, max=els.lines[linename].up)
                line_model.set_param_hint(name="%s_amplitude" % linename, expr='%s_amplitude * %s_factor' % (refline, linename))
    # mask the range we want to use
    fluxes = ma.masked_where((bin_spectrum[:,0] > wav_ranges[0][0]) & (bin_spectrum[:,0] < wav_ranges[0][1]), bin_spectrum[:,1])
    for wav_range in wav_ranges[1:]:
        fluxes = ma.masked_where((bin_spectrum[:,0] > wav_range[0]) & (bin_spectrum[:,0] < wav_range[1]), fluxes)
    #lmfit does not work well with masked values
    wavelengths = bin_spectrum[:, 0][np.where(fluxes.mask)]
    std_dev = bin_spectrum[:, 2][np.where(fluxes.mask)]
    weights = 1.0 / std_dev
    fluxes = ma.getdata(fluxes)[np.where(fluxes.mask)]

    result = line_model.fit(fluxes, x=wavelengths, weights=weights, nan_policy='omit')
    if plot:
        fig, axes = plt.subplots(1, len(wav_ranges), sharey=True)
        fig.suptitle("SN = %.1f, (X, Y) = (%.2f, %.2f)" % (bin_sn[bin_num], x_bar[bin_num], y_bar[bin_num]))

        for ax, wave_range in zip(axes, wav_ranges):
            ax.errorbar(wavelengths, fluxes)
            ax.plot(wavelengths, get_inital_fit(line_model, wavelengths), color='orange', ls='--')
            ax.plot(wavelengths, result.best_fit, color="green")
            ax.set_xlim(wave_range[0], wave_range[1])  # outliers only
        for ax in axes[1:]:
            ax.spines['left'].set_visible(False)
        for ax in axes[:-1]:
            ax.spines['right'].set_visible(False)
        #ax.tick_params(labeltop=False)  # don't put tick labels at the top
        fig.savefig("%s/bin_%d.png" % (outpath, bin_num))
        plt.close()


class line:
    """This class is the basic for defining a line. A line is defined by its name, its wavelength and the reference line to which it is attached if the ratio has to be constrained.

    """
    def __init__(self, name, wave, ref=None, fit=False, low=0, up=None, th=None, index=None):
        self.index = index
        self.name = name
        self.wave = wave
        self.fit = fit
        self.ref = ref
        if ref == name:
            self.low = 1
            self.up = 1
            self.th = 1
        else:
            self.low = low
            self.up = up
            self.th = th

class lines:
    """This class enables to deal with lines.  A dictionary stored in lines will contain informations on each lines.

    """

    def append(self, line):
        self.lines[line.name] = line
        self.lines[line.name].index = self.index
        self.index += 1

    def __init__(self):
        self.index = 0
        self.lines = {}
        # ref always to the higher wavelength line
        self.append(line('HBETA', 4859))
        self.append(line('HALPHA', 6563, ref='HBETA', low=2.65, th=2.85, up=6))
        self.append(line('NII6548', 6548.05))
        self.append(line('NII6583', 6583.45, ref='NII6548', low=2.5, up=3.3, th=3.))
        # there ratios are the inverse of 0.4 and 1.5 respectively
        self.append(line('SII6716', 6716.44))
        self.append(line('SII6731', 6730.82, ref='SII6716', low=0.667, up=2.5, th=1.))
        self.append(line('OIII4959', 4958.911))
        self.append(line('OIII5007', 5006.843, ref='OIII4959', low=2.5, up=3.3, th=3.))
        self.append(line('OI6300', 6300.3))
        #self.append(line("NII5755", 5754.8))
        #self.append(line("HeII",   4685.682))


def get_inital_fit(model, x):
    params = model.make_params()
    return model.eval(params, x=x)

def create_line_model(els, redshift=0.001568, sigma=1.4 * u.AA, margin=10 * u.AA):
    """Create lmfit model with all the lines to be fitted. Lines sigmas and position are tight to the reference values of other lines.
    Parameters
    ---------
    els: lines
    redshift: float
        Initial guess for the redshift
    sigma: astropy.Quantity
        Initial guess for the line sigma
    margin: astropy.Quantity
        +-Margin to considered around the redshifted line centroid for the min and maximum allowed values
    """
    line_model = None
     # initialize model with first line
    for i in els.lines:
        linename = els.lines[i].name
        if line_model is None:
            line_model = GaussianModel(prefix="%s_" % linename)
        else:
            line_model += GaussianModel(prefix="%s_" % linename)
    # assign parameters
    argsorted = np.array([els.lines[i].index for i in els.lines])
    linenames = np.array([els.lines[i].name for i in els.lines])[argsorted]
    margin = margin.to(u.AA).value
    sigma = sigma.to(u.AA).value

    for i in linenames:
        line = els.lines[i]
        redshifted = (1 + redshift) * line.wave

        line_model.set_param_hint("%s_center" % line.name, value=redshifted, vary=True, max=redshifted + margin,
                                      min=redshifted - margin)
        if line.ref is None:

            line_model.set_param_hint("%s_sigma" % line.name, value=sigma, vary=True, max=4, min=0.8)
        # tight lines center and sigma to the reference lines
        else:
            refline = els.lines[line.ref]
            line_model.set_param_hint("%s_center" % line.name, expr="%s_center + %.1f" % (refline.name, (line.wave - refline.wave)))
            line_model.set_param_hint("%s_sigma" % line.name, expr="%s_sigma" % refline.name)

    return line_model

# method take from ifuanal but modified (the output spectrum was multiplied by the pixel size for ome odd reason)
def get_weighted_spectrum(data_cube, x, y):
    """
    Return the weighted mean spectrum for spaxels at locations ``x``,
    ``y``.

    Similar to ``_get_single_spectrum`` except ``x`` and ``y`` are arrays.
    The single spectra given by these locations are combined using the
    weighted mean of the fluxes. Returns array of same form as
    ``_get_single_spectrum``.

    """
    x = np.ravel(np.asarray(x))
    y = np.ravel(np.asarray(y))
    if x.size != y.size:
        raise AttributeError("``x`` and ``y`` should be same size")
    # Use weighted arthimetic mean and variance of weighted mean
    # arrays to hold all the flux and flux_stddev values in spaxels
    spaxels_flux = data_cube.data[:,y,x] # numpy axes switch
    spaxels_stddev = data_cube.var[:,y,x]

    # Find any wavelengths where we have >75% spaxels as nans
    # and flag that wavelength as bad for starlight
    bad_idx = np.isnan(spaxels_flux) | np.isnan(spaxels_stddev)
    num_bad = np.sum(bad_idx, axis=1)
    bad_lamb = num_bad > 0.75 * x.size

    # Array to hold final weighted-mean spectrum - same format
    # as _get_single_spectrum()
    spec = np.empty((data_cube.wave.coord().size, 4))
    spec[:,0] = data_cube.wave.coord()
    # Use masked arrays to cover the nans while preserving shape
    spaxels_flux_ma = np.ma.masked_array(spaxels_flux, bad_idx)
    spaxels_stddev_ma = np.ma.masked_array(spaxels_stddev, bad_idx)
    # Calculate the weighted mean and uncertainty
    w = 1/spaxels_stddev_ma**2
    spec[:,1] = np.ma.average(spaxels_flux_ma, weights=w, axis=1)

    #* x.size --> why multiply by the size of the pixel?
    spec[:,2] = 1/np.sum(w, axis=1)**0.5
    #* x.size
    # STARLIGHT ignores flags >=2
    spec[:,3] = bad_lamb.astype("int") * 2

    return spec

def create_out_maps(lines, shape):
    """Create output maps"""
    outmaps = { }
    outpars = ["amplitude", "fwhm", "center", "snr"]
    for line in lines:
        for par in outpars:
            outmaps["%s_%s" % (line, par)] = np.ones(shape) * np.nan
            outmaps["%s_e%s" % (line, par)] = np.ones(shape) * np.nan
    return outmaps


parser = argparse.ArgumentParser(description='Fit data cube using voronoi binning')
parser.add_argument("-z", "--redshift", nargs='?', help="Initial guess for the redshift", default=0.001568, type=float)
parser.add_argument("-s", "--sigma", nargs='?', help="Initial guess for the width (sigma) of the lines in Angstroms. Default 1.4 Angstroms", default=1.4, type=float)
parser.add_argument("-t", "--target_sn", nargs='?', help="Target S/N for the Voronoi binning", default=5, type=float)
parser.add_argument("-c", "--cpu", nargs='?', help="Number of CPUs to use for the Voronoi binning", default=6, type=float)
parser.add_argument("--min_sn", nargs='?', help="Minimum S/N to reject pixels for the Voronoi binning", default=0.1, type=float)
parser.add_argument("input_cube", nargs=1, help="Path to input cube to be analysed", type=str)
parser.add_argument("--plot", help="Plot resulting fits of each inidivual Voronoi pixel. Default false", action='store_true')
parser.add_argument("-o", "--outdir", nargs='?', help="Name of the output directory", default="voronoi", type=str)

args = parser.parse_args()

outpath = args.outdir

if not os.path.isdir(outpath):
    os.mkdir(outpath)



# hb [OIII], [OI]6300, [NII] Ha and [SII]
groups = [["HBETA", "OIII4959", "OIII5007"], ["OI6300"],
          ["NII6548","HALPHA", "NII6583"], ["SII6716", "SII6731"]]
          # ["NII5755"],
          # "HeII",
els = lines()
redshift = args.redshift
# Margins to cut the spectrum around the blue and red lines of each group
margin = 25 * u.AA
wav_ranges = np.zeros((len(groups), 2))
for index, group in enumerate(groups):
    minwav = els.lines[group[0]].wave * (1 + redshift) - margin.value
    maxwav = els.lines[group[-1]].wave * (1 + redshift) + margin.value
    wav_ranges[index] = (minwav, maxwav)

min_sn = args.min_sn
target_sn = args.target_sn
# read cube
cube = args.input_cube[0]
data_cube = Cube(cube)

# signal
lambda_low =  6550 * u.AA
lambda_high = 6600 * u.AA
#lambda_low = 5008 * u.AA
#lambda_high = 5020* u.AA
#print(lambda_low)
wavelengths = data_cube.wave.coord()
idx_low = np.abs(wavelengths - lambda_low.to(u.AA).value).argmin()
idx_upp = np.abs(wavelengths - lambda_high.to(u.AA).value).argmin() + 1
# catching the warning when spaxels only have nans
signal = np.nanmean(data_cube.data[idx_low:idx_upp, :, :],
                    axis=0)
noise = np.nanmean(np.sqrt(data_cube.var[idx_low:idx_upp, :, :]), axis=0)

# continuum
#lambda_low =  4835 * u.AA
#lambda_high = 4855 * u.AA
#lambda_low =  6450 * u.AA
#lambda_high = 6520 * u.AA
#cont_idx_low = np.abs(wavelengths - lambda_low.to(u.AA).value).argmin()
##cont_idx_upp = np.abs(wavelengths - lambda_high.to(u.AA).value).argmin() + 1
#cont = np.nanmean(data_cube.data[cont_idx_low:cont_idx_upp, :, :], axis=0)
#signal = signal - cont
num_nans = np.sum(np.isnan(data_cube.data[idx_low:idx_upp, :, :]),
                  axis=0)
perc = 0.2
idx_nans = np.where(num_nans > perc * (idx_upp - idx_low))
print("Ignoring %d pixels with %d%% of nan pixels" % (len(idx_nans), perc * 100))
signal[idx_nans] = np.nan
noise[idx_nans] = np.nan
xx, yy = np.meshgrid(np.arange(signal.shape[1]),
                             np.arange(signal.shape[0]))
# Need to clean the input of the bad spectra:
#  remove those with signal == nan
vor_input = np.column_stack((xx.ravel(), yy.ravel(),
                                     signal.ravel(), noise.ravel()))
vor_input = vor_input[~np.isnan(vor_input[:,3])]
#  remove any with negative or zero noise
vor_input = vor_input[vor_input[:,3] > 0]
#  also reduce to only spaxels with S/N >= min_sn
vor_input = vor_input[(vor_input[:,2]/vor_input[:,3]) >= min_sn]
x, y, sig, noi = vor_input[:, [0,1,2,3]].T
vor_plot = "%s/%s" %(outpath, cube.replace(".fits","_%d_bins_voronoi.pdf" % target_sn))
start = time.time()
res  = voronoi.voronoi_2d_binning(x, y, sig, noi, targetSN=target_sn,
                                          cvt=True, plot=vor_plot, pixelsize=1,
                                          quiet=False, n_cpu=7)
end = time.time()
print("Voronoi bining completed in %.2f mins" % ((end-start) / 60))
bin_num, x_node, y_node, x_bar, y_bar, bin_sn, n_pix, scale = res

outfile = cube.replace(".fits", "_%d_bins_voronoi.dat" % target_sn)
print(outfile)
with open(outfile, "w+") as file:
    file.write("#bin_num\tx\ty\n")
    [file.write("%i\t%i\t%i\n" % (bin, x_, y_)) for bin, x_, y_ in zip(bin_num, x, y)]
vor_output = np.column_stack([x, y, x_bar[bin_num], y_bar[bin_num], bin_num])
end = time.time()
print("Voronoi bining completed in %.2f mins" % ((end-start) / 60))
bin_num, x_node, y_node, x_bar, y_bar, bin_sn, n_pix, scale = res

outfile = cube.replace(".fits", "_%d_bins_voronoi.dat" % target_sn)
out_bins = "#bin_num\tx\ty\n"
strings = [("%i\t%i\t%i" % (bin, x_, y_)) for bin, x_, y_ in zip(bin_num, x, y)]
strings  = "\n".join(strings)
out_bins += strings
with open(outfile, "w+") as file:
    file.write(out_bins)
# Fitting procedure
outshape = data_cube[0 ,: , :].shape
outmaps = create_out_maps(els.lines, outshape)
outmaps["redchi"] = np.ones(outshape) * np.nan
outmaps["residuals"] = np.ones(data_cube.shape) * np.nan

outpars = ["amplitude", "fwhm", "center"]
outfile = data_cube[0 ,: , :]
sigma = args.sigma
line_model = create_line_model(els)
cont_prefix = 'cont'
cont_model = PolynomialModel(degree=1, prefix="%s_" %cont_prefix)
line_model += cont_model
#residual_cube[:, y, x] = result.residual
for bn in np.sort(np.unique(bin_num)).astype("int"):
    print(f'\rProcessing bin %d/%d' % (bn , len(bin_sn)), end="")
    # take all indexes with bin = e.g. 0
    idx = vor_output[:,4] == bn
    x, y, x_mean, y_mean, _bn = vor_output[idx].T
    xint, yint= x.astype("int"), y.astype("int")
    xrav = np.ravel(np.asarray(xint))
    yrav = np.ravel(np.asarray(yint))
    bin_spectrum = get_weighted_spectrum(data_cube, xint, yint)
    mean, median, stddev = sigma_clipped_stats(bin_spectrum[:, 1], sigma=3)
    line_model.set_param_hint("%s_c0" % cont_prefix, value=median, vary=True, min=0)
    line_model.set_param_hint("%s_c1" % cont_prefix, value=0, vary=True)

    for index, wav_range in enumerate(wav_ranges):
        fluxes = ma.masked_where((bin_spectrum[:, 0] < wav_range[0]) & (bin_spectrum[:,0] > wav_range[1]), bin_spectrum[:,1])
        for linename in groups[index]:
            # assign sigma and amplitude base on data
            refline = els.lines[linename].ref
            if refline is None:
                centroid = els.lines[linename].wave * (1 + redshift)
                central_wavelength = np.argmin(np.abs(bin_spectrum[:,0] - centroid))
                min_ind = central_wavelength - 4
                max_ind = central_wavelength + 4
                peak_index = np.argmax(bin_spectrum[min_ind:max_ind,1]) + min_ind
                peak_wavelength = bin_spectrum[peak_index, 0]
                init_amplitude = (bin_spectrum[peak_index, 1] - median) * np.sqrt(2 * np.pi * (sigma ** 2))
                line_model.set_param_hint("%s_amplitude" % linename, value=init_amplitude, vary=True, min=0)
            else:
                line_model.set_param_hint("%s_factor" % linename, value=els.lines[linename].th, min=els.lines[linename].low, max=els.lines[linename].up)
                line_model.set_param_hint(name="%s_amplitude" % linename, expr='%s_amplitude * %s_factor' % (refline, linename))
    # mask the range we want to use
    fluxes = ma.masked_where((bin_spectrum[:,0] > wav_ranges[0][0]) & (bin_spectrum[:,0] < wav_ranges[0][1]), bin_spectrum[:,1])
    for wav_range in wav_ranges[1:]:
        fluxes = ma.masked_where((bin_spectrum[:,0] > wav_range[0]) & (bin_spectrum[:,0] < wav_range[1]), fluxes)
    #lmfit does not work well with masked values
    wavelengths = bin_spectrum[:, 0][np.where(fluxes.mask)]
    std_dev = bin_spectrum[:, 2][np.where(fluxes.mask)]
    weights = 1.0 / std_dev
    fluxes = ma.getdata(fluxes)[np.where(fluxes.mask)]

    result = line_model.fit(fluxes, x=wavelengths, weights=weights, nan_policy='omit')
    if args.plot:
        fig, axes = plt.subplots(1, len(wav_ranges), sharey=True)
        fig.suptitle("SN = %.1f, (X, Y) = (%.2f, %.2f)" % (bin_sn[bn], x_bar[bn], y_bar[bn]))

        for ax, wave_range in zip(axes, wav_ranges):
            ax.errorbar(wavelengths, fluxes)
            ax.plot(wavelengths, get_inital_fit(line_model, wavelengths), color='orange', ls='--')
            ax.plot(wavelengths, result.best_fit, color="green")
            ax.set_xlim(wave_range[0], wave_range[1])  # outliers only
        for ax in axes[1:]:
            ax.spines['left'].set_visible(False)
        for ax in axes[:-1]:
            ax.spines['right'].set_visible(False)
        #ax.tick_params(labeltop=False)  # don't put tick labels at the top
        fig.savefig("%s/bin_%d.png" % (outpath, bn))
        plt.close()
    #print("FIT FOR BIN %d" % i)
    #print(result.fit_report())
    best_fit_params = result.params
    compute_errors = False
    if result.errorbars and compute_errors:
        conf = result.conf_interval(sigmas=[1])
    # store params in maps
    # the signal is the same, whether we do the integral ourselves or we take it directly from the height value
    for line in els.lines:
        signal = best_fit_params.get("%s_height" % (line))
        central_wavelength_idx = np.argmin(np.abs(wavelengths - best_fit_params.get("%s_center" % (line))))
        # central wavelength for HB is 20
        noise = np.sqrt(np.nansum((fluxes - result.best_fit)[central_wavelength_idx - 4: central_wavelength_idx + 4] ** 2))
        outmaps["%s_snr" %(line)][yrav, xrav] = signal / noise
        for par in outpars:
            param = "%s_%s" % (line, par)
            outmaps[param][yrav, xrav] = best_fit_params.get(param).value
            if result.errorbars  and compute_errors:
                if param in conf:
                    outmaps["%s_e%s" % (line, par)][yrav, xrav] = np.abs(conf[param][0][1] - conf[param][1][1])
                else:
                    outmaps["%s_e%s" % (line, par)][yrav, xrav] = np.nan

            else:
                outmaps["%s_e%s" % (line, par)][yrav, xrav] = np.nan
    outmaps["redchi"][yrav, xrav] = result.redchi
    #outmaps["residuals"][:, yrav, xrav] = result.residual

for line in els.lines:
    for par in outpars:
        # parameter value
        outfile.data = outmaps["%s_%s" % (line, par)]
        outfile.write("%s/%s_%s.fits" % (outpath, line, par))
        # parameter uncertainty
        outfile.data = outmaps["%s_e%s" % (line, par)]
        outfile.write("%s/%s_e%s.fits" % (outpath, line, par))
    outfile.data = outmaps["%s_snr" % (line)]
    outfile.write("%s/%s_snr.fits" % (outpath, line))
outfile.data = outmaps["redchi"]
outfile.write("%s/redchi.fits" % outpath)
print("Output stored in %s" % outpath)
