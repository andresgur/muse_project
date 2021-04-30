# @Author: Andrés Gúrpide <agurpide>
# @Date:   30-04-2021
# @Email:  agurpidelash@irap.omp.eu
# @Last modified by:   agurpide
# @Last modified time: 30-04-2021
# Config file for the BPT diagram code


lineratio_paths = {
    "o3_hb": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted_hst_uncorrected/lineratios/OIII5007_HBETAratio.fits",
    "n2_ha": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted_hst_uncorrected/lineratios/NII6583_HALPHA.fits",
    "oI_ha": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted_hst_uncorrected/lineratios/OI6300_HALPHA.fits",
    "s2_ha": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted_hst_uncorrected/lineratios/SII_HALPHA.fits"
}
# type can be law or kewley
bpt_type = {
    "type": "law"
}
# usually to give intensity to the individual points in the BPT, it can be left blank
dispersion_fits = {
    "file": "/home/agurpide/optical_data/NGC1313/nancleancubes/coordadjusted_hst_uncorrected/cleancamel_1_n2ha_ssmooth_wavedisp_indep_HALPHAcorrected"
}
