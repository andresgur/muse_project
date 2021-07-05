# @Author: Andrés Gúrpide <agurpide>
# @Date:   30-04-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 30-04-2021
# Config file for the BPT diagram code


lineratio_paths = {
    "o3_hb": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted/lineratios/OIII5007_HBETAratio.fits",
    "n2_ha": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted/lineratios/NII6583_HALPHAratio.fits",
    "oI_ha": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted/lineratios/OI6300_HALPHAratio.fits",
    "s2_ha": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted/lineratios/SII_HALPHAratio.fits"
}
# type can be law (Law et al. 2021) or kewley (Kewley et al. 2006)
type = "law"

# path for a 2D fits file to give intensity to the individual points in the BPT, it can be set to None and it will be ignored
z_values = "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted/camel_1_n2ha/cleaned_images/cleancamel_1_n2ha_ssmooth_disp_indep_HALPHA.fits"
