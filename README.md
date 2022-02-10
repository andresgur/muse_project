# muse_project

Various scripts for multi-wavelength analysis of ULX and their environments

~~~~~~~~~~~~~~~~~~
~    Chandra     ~
~~~~~~~~~~~~~~~~~~
marxsim_script.py --> Script to be run from the command line to simulate N PSFs based on a given chandra observation. Ciao and marx need to be properly setup before running the script. Use python marxsim_script.py -h to obtain the available command line arguments

check_extension.py --> Script to be run after marxsim_script.py has been run to analyse the output. Use python marxsim_script.py -h to obtain the available command line arguments

~~~~~~~~~~~~~~~~~
~      HST      ~ 
~~~~~~~~~~~~~~~~~
aperture_photometry.py --> Script to retrieve magnitudes, fluxes from a given ds9 region (it automatically reads the necessary keywords from the header to perform the corrections). A finite aperture correction is needed for small PSFs. Background subtraction is also possible (must be added a second line in the ds9 region file).

astro_corr.py --> Script to refine the astrometry of several HST images (more than one is possible, and the solution will be global to all of them) using the gaia catalogue. It uses tweakreg https://drizzlepac.readthedocs.io/en/latest/tweakreg.html.

~~~~~~~~~~~~~~~~~
~      MUSE     ~
~~~~~~~~~~~~~~~~~
----reddening----

deredden.py --> Script to deredden flux line maps, it needs as input the balmer decrement map (Halpha/Hbeta) and it will automatically localize the flux maps

deredden_momcheva.py --> Script to deredden flux line maps based on appendix of http://arxiv.org/abs/1207.5479, it needs as input the hbeta and halpha directories and it will automatically look for the line flux maps and deredden them (python deredden_momcheva.py -a halpha_path -b hbeta_path -Rv 4.05 -i 2.86). It creates a map of E(B-V) and one with the uncertainty on E(B-V)

----utils----

cutfromregion.py --> Script to cut images from an input region ds9 file.

image_stats.py --> Retrieve statistical information (min, max, mean) of a region or entire image by default

----MAPPINGS-----

read_mappings.py --> Script to obtain predicted line ratios from the mapping libraries (python read_mappings.py <files> (for instance V_*_b0_001_s_lines.txt V_*b1_s_lines.txt [MVQP]_*b[0e]_s_lines.txt [MVQ]_*b10_s_lines.txt T_*b0_001_s_lines.txt T_*b1_s_* T_*b10_s_lines.txt). The script will read the paths to the BPT diagrams from bpt_config.py file

read_mappingsV.py --> Script to obtain predicted line ratios from the mapping V libraries. The script will read the paths to the BPT diagrams from bpt_config.py file

-------BPT-------

bpt_newer.py -> bpt diagram from Law et al. 2021, now with geometry

~~~~~~~~~~~~~~~~~
~      Misc     ~
~~~~~~~~~~~~~~~~~
image_stats.py --> Script to obtain some stats from a particular region of an image (python image_stats.py image.fits -region ds9reigonfile.reg

