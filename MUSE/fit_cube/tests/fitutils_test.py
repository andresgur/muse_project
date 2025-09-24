#!/usr/bin/env python
# -*- coding: utf-8 -*
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lmfit.models import PolynomialModel, GaussianModel
import numpy as np
import unittest
import matplotlib.pyplot as plt
from fitutils import fit_spectrum, SpectrumMaskedError
from mpdaf.obj import WaveCoord, Spectrum
from line import Lines

CATALOG_LINES = Lines()

class TestFitLineUtils(unittest.TestCase):
    
    def test_fit_spectrum_flat_cont(self):
        step = 1.25
        crval  = 4599.6875
        wave1 = WaveCoord(cdelt=step, crval=crval, cunit="angstrom", shape=3801)
        cont_prefix = 'cont_'
        model = PolynomialModel(degree=0, prefix=cont_prefix) + GaussianModel(prefix='HBETA_')
        contvalue = 10
        model.set_param_hint('cont_c0', value=contvalue)
        model.set_param_hint('HBETA_amplitude', value=100)
        redshift = 0.001676
        model.set_param_hint('HBETA_center', value=4861.33 * (1 + redshift))
        model.set_param_hint('HBETA_sigma', value=3)
        params = model.make_params()
        fluxes = model.eval(x=wave1.coord(), params=params)
        noise = np.random.normal(0, 0, size=fluxes.size)
        fluxes += noise
        data = Spectrum(data=fluxes, wave=wave1, var=noise**2)
        data.wave = wave1


        fit_lines = {"HBETA": CATALOG_LINES.lines["HBETA"]}

        fit_results, conf = fit_spectrum(data, redshift=redshift, fit_lines=fit_lines, sigma=3, degree=0, wavelengths=wave1.coord())

        best_values = fit_results.params
        param = best_values[f"{cont_prefix}c{0}"]
        self.assertAlmostEqual(param.value, contvalue, delta=1)
        param = best_values[f"HBETA_amplitude"]
        self.assertAlmostEqual(param.value, 100, delta=1)
        param = best_values[f"HBETA_center"]
        self.assertAlmostEqual(param.value, 4861.33 * (1 + redshift), delta=1)
        param = best_values[f"HBETA_sigma"]
        self.assertAlmostEqual(param.value, 3, delta=1)
        param = best_values[f"HBETA_vel"]
        self.assertAlmostEqual(param.value, 0, delta=1)



    def test_fit_spectrum_linear_cont(self):
            step = 1.25
            crval  = 4599.6875
            wave1 = WaveCoord(cdelt=step, crval=crval, cunit="angstrom", shape=3801)
            cont_prefix = 'cont_'
            model = PolynomialModel(degree=1, prefix=cont_prefix) + GaussianModel(prefix='HBETA_')
            contvalue = 33144
            model.set_param_hint('cont_c0', value=contvalue)
            cont2value = -4
            model.set_param_hint('cont_c1', value=cont2value)
            hb_amplitude = 73837
            model.set_param_hint('HBETA_amplitude', value=hb_amplitude)
            redshift = 0.001676
            model.set_param_hint('HBETA_center', value=4861.33 * (1 + redshift))
            model.set_param_hint('HBETA_sigma', value=3)
            params = model.make_params()
            print(params)
            fluxes = model.eval(x=wave1.coord(), params=params)
            noise = np.random.normal(0, 0, size=fluxes.size)
            fluxes += noise
            data = Spectrum(data=fluxes, wave=wave1, var=noise**2)
            data.wave = wave1

            data.plot()
            #plt.show()

            fit_lines = {"HBETA": CATALOG_LINES.lines["HBETA"]}

            fit_results, conf = fit_spectrum(data, redshift=redshift, fit_lines=fit_lines, sigma=3, degree=1, wavelengths=wave1.coord())

            best_values = fit_results.params
            param = best_values[f"{cont_prefix}c{0}"]
            self.assertAlmostEqual(param.value, contvalue, delta=1)
            param = best_values[f"{cont_prefix}c{1}"]
            self.assertAlmostEqual(param.value, cont2value, delta=1)
            param = best_values[f"HBETA_amplitude"]
            self.assertAlmostEqual(param.value, hb_amplitude, delta=1)
            param = best_values[f"HBETA_center"]
            self.assertAlmostEqual(param.value, 4861.33 * (1 + redshift), delta=1)
            param = best_values[f"HBETA_sigma"]
            self.assertAlmostEqual(param.value, 3, delta=1)
            param = best_values[f"HBETA_vel"]
            self.assertAlmostEqual(param.value, 0, delta=1)



if __name__ == '__main__':
    unittest.main()