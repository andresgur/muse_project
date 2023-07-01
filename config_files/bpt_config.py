# @Author: Andrés Gúrpide <andresgur>
# @Date:   30-04-2021
# @Email:  andresgurlash@irap.omp.eu
# @Last modified by:   andresgur
# @Last modified time: 30-04-2021
# Config file for the BPT diagram code

[Paths]
ulx_dir:/home/andresgur/optical_data/NGC247/
base_dir:${ulx_dir}/lineratios
o3_hb: ${base_dir}/OIII5007_HBETAratio.fits
n2_ha: ${base_dir}/NII6583_HALPHAratio.fits
oI_ha: ${base_dir}/OI6300_HALPHAratio.fits
s2_ha: ${base_dir}/SII_HALPHAratio.fits


[bpt]
# diagram can be law (Law et al. 2021) or kewley (Kewley et al. 2006)
diagram = law
[z_values]
# path for a 2D fits file to give intensity to the individual points in the BPT, it can be set to None and it will be ignored
z_values = ${Paths:ulx_dir}/camel_1_n2ha/cleaned_images/cleancamel_1_n2ha_ssmooth_wavedisp_indep_HALPHAcorrected_wcsupdated.fits
