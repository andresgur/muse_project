#!/usr/bin/env python
# -*- coding: utf-8 -*
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lineutils import peculiar_velocity, compute_shift, ckms, get_instrument_FWHM, correct_FWHM
import numpy as np
import unittest
import matplotlib.pyplot as plt

class TestLineUtils(unittest.TestCase):
    
    def test_zero_velocity(self):
        """Test that zero peculiar velocity gives observed wavelength equal to reference."""
        reference = 6562.8  # H-alpha
        z_sys = 0.001
        
        # When velocity = 0, observed wavelength should be reference * (1 + z_sys)
        wavelength = reference * (1 + z_sys)
        
        vp, evp = peculiar_velocity(reference, wavelength, z_sys=z_sys)

        self.assertAlmostEqual(vp, 0, delta=1e-10, msg=f"Expected velocity ~0, got {vp}")
        self.assertEqual(evp, 0, "Error should be 0")


    def test_error_propagation_zero_velocity(self):
        """Test that zero peculiar velocity gives observed wavelength equal to reference."""
        reference = 6562.8  # H-alpha
        
        z_sys = 0.001
        
        # When velocity = 0, observed wavelength should be reference * (1 + z_sys)
        wavelength = reference * (1 + z_sys)
        ewavelength = 0.11
        
        vp, evp = peculiar_velocity(reference, wavelength, z_sys=z_sys, ewavelength=ewavelength)
        wavelengths = np.random.normal(wavelength, ewavelength, size=1000000)
        velocities = [peculiar_velocity(reference, wl, z_sys=z_sys)[0] for wl in wavelengths]

        perce = np.percentile(velocities, [50- 68.269 / 2., 50])
        evp_mc = perce[1] - perce[0]

        self.assertAlmostEqual(vp, 0, delta=1e-10, msg=f"Expected velocity {0}, got {vp}")
        self.assertAlmostEqual(evp, evp_mc, delta=1e-2, msg="Error is not working!")

    
    def test_error_propagation_velocity(self):
        """Test that zero peculiar velocity gives observed wavelength equal to reference."""
        reference = 6562.8  # H-alpha
        
        z_sys = 0.001
        
        # When velocity = 0, observed wavelength should be reference * (1 + z_sys)
        wavelength = reference * (1 + z_sys) + 0.5
        z = wavelength / reference - 1

        zp = (1 + z) / (1 + z_sys) - 1
        # non relativistic case, but should be ok for small vps
        vpexpected = zp * ckms
        
        ewavelength = 0.11
        vp, evp = peculiar_velocity(reference, wavelength, z_sys=z_sys, ewavelength=ewavelength)
        wavelengths = np.random.normal(wavelength, ewavelength, size=1000000)
        velocities = [peculiar_velocity(reference, wl, z_sys=z_sys)[0] for wl in wavelengths]

        perce = np.percentile(velocities, [50- 68.269 / 2., 50])
        evp_mc = perce[1] - perce[0]

        self.assertAlmostEqual(vp, vpexpected, delta=1e-3, msg=f"Expected velocity {vpexpected}, got {vp}")
        self.assertAlmostEqual(evp, evp_mc, delta=1e-2, msg="Error propagation is not working in peculiar_velocity!")

    def test_compute_shift_zero_velocity(self):
        """Test that zero velocity gives observed wavelength equal to reference * (1 + z_sys)."""
        reference = 6562.8  # H-alpha
        z_sys = 0.001
        velocity = 0.0
        
        # When velocity = 0, observed wavelength should be reference * (1 + z_sys)
        expected_wavelength = reference * (1 + z_sys)
        
        wavelength, ewavelength = compute_shift(reference, z_sys, velocity)
        
        self.assertAlmostEqual(wavelength, expected_wavelength, delta=1e-10, 
                              msg=f"Expected wavelength {expected_wavelength}, got {wavelength}")
        self.assertEqual(ewavelength, 0, "Error should be 0 when velocity error is 0")

    def test_compute_shift_error_propagation_zero_velocity(self):
        """Test error propagation for compute_shift with zero velocity."""
        reference = 6562.8  # H-alpha
        z_sys = 0.001
        velocity = 0.0
        evelocity = 5.0  # km/s
        
        wavelength, ewavelength = compute_shift(reference, z_sys, velocity, evelocity)
        
        # Monte Carlo validation
        velocities = np.random.normal(velocity, evelocity, size=1000000)
        wavelengths = [compute_shift(reference, z_sys, v)[0] for v in velocities]
        
        perce = np.percentile(wavelengths, [50 - 68.269 / 2., 50])
        ewavelength_mc = perce[1] - perce[0]
        
        expected_wavelength = reference * (1 + z_sys)
        self.assertAlmostEqual(wavelength, expected_wavelength, delta=1e-10, 
                              msg=f"Expected wavelength {expected_wavelength}, got {wavelength}")
        self.assertAlmostEqual(ewavelength, ewavelength_mc, delta=1e-2, 
                              msg="Error propagation is not working in compute_shift!")

    def test_compute_shift_error_propagation_velocity(self):
        """Test error propagation for compute_shift with non-zero velocity."""
        reference = 6562.8  # H-alpha
        z_sys = 0.001
        velocity = 100.0  # km/s
        evelocity = 5.0  # km/s
        
        wavelength, ewavelength = compute_shift(reference, z_sys, velocity, evelocity)
        
        # Monte Carlo validation
        velocities = np.random.normal(velocity, evelocity, size=1000000)
        wavelengths = [compute_shift(reference, z_sys, v)[0] for v in velocities]
        
        perce = np.percentile(wavelengths, [50 - 68.269 / 2., 50])
        ewavelength_mc = perce[1] - perce[0]
        
        # For small velocities, approximate expected wavelength
        v_c = velocity / ckms
        zp_plus_1 = ((1 + v_c) / (1 - v_c)) ** 0.5
        expected_wavelength = reference * zp_plus_1 * (1 + z_sys)
        
        self.assertAlmostEqual(wavelength, expected_wavelength, delta=1e-6, 
                              msg=f"Expected wavelength {expected_wavelength}, got {wavelength}")
        self.assertAlmostEqual(ewavelength, ewavelength_mc, delta=1e-2, 
                              msg="Error propagation is not working in compute_shift!")
        
    def get_instrument_FWHM_error(self):
        # Example test for get_instrument_FWHM
        wavelength = 7000.0
        ewavelength = 0.5

        wavelengths = np.random.normal(wavelength, ewavelength, size=1000000)

        FWHMs = [get_instrument_FWHM(wl)[0] for wl in wavelengths]
        perce = np.percentile(FWHMs, [50 - 68.269 / 2., 50])
        eFWHM_mc = perce[1] - perce[0]

        FWHM, eFWHM = get_instrument_FWHM(wavelength, ewavelength)

        self.assertAlmostEqual(eFWHM, eFWHM_mc, delta=1e-2, 
                               msg="Error propagation is not working in get_instrument_FWH!")


    def test_correct_FWHM_err_FWHM_obs(self):
        FWHM_obs = 350.0
        eFWHM_obs = 10.0
        FWHM_inst = 80.0
        eFWHM_inst = 0

        FWHM_vel_corrected, eFWHM_vel_corrected = correct_FWHM(FWHM_obs, FWHM_inst, eFWHM_obs, eFWHM_inst)

        # Monte Carlo validation
        FWHM_obs_samples = np.random.normal(FWHM_obs, eFWHM_obs, size=10000000)

        FWHM_corrected = correct_FWHM(FWHM_obs_samples, FWHM_inst, 0, 0)[0]

        perce = np.percentile(FWHM_corrected, [50 - 68.269 / 2., 50])
        eFWHM_mc = perce[1] - perce[0]

        self.assertAlmostEqual(eFWHM_vel_corrected, eFWHM_mc, delta=1e-2, 
                               msg="Error propagation in correct_FWHM is not working when eFWHM_obs > 0!")
        

    def test_correct_FWHM_err_FWHM_inst(self):
        FWHM_obs = 350.0
        eFWHM_obs = 0.0
        FWHM_inst = 140.0
        eFWHM_inst = 5.0

        FWHM_vel_corrected, eFWHM_vel_corrected = correct_FWHM(FWHM_obs, FWHM_inst, eFWHM_obs, eFWHM_inst)

        # Monte Carlo validation
        FWHM_inst_samples = np.random.normal(FWHM_inst, eFWHM_inst, size=15000000)

        FWHM_corrected = correct_FWHM(FWHM_obs, FWHM_inst_samples, 0, 0)[0]

        #plt.hist(FWHM_corrected)
        #plt.show()
        perce = np.percentile(FWHM_corrected, [50 - 68.269 / 2., 50, 50 + 68.269 / 2.])
        eFWHM_mc = perce[1] - perce[0]
        eFWHM_mc2 = (perce[2] - perce[1]) 

        self.assertAlmostEqual(eFWHM_vel_corrected, (eFWHM_mc + eFWHM_mc2) / 2., delta=1e-2, 
                            msg="Error propagation in correct_FWHM not working when eFWHM_inst >0!")



if __name__ == '__main__':
    unittest.main()