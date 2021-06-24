# muse_project
Scripts to deal with some subtleties regarding MUSE data

bpt_newer.py -> bpt diagram from Law et al. 2021, now with geometry

marxsim_script.py --> Script to be run from the command line to simulate N PSFs based on a given chandra observation. Ciao and marx need to be properly setup before running the script. Use python marxsim_script.py -h to obtain the available command line arguments

check_extension.py --> Script to be run after marxsim_script.py has been run to analyse the output. Use python marxsim_script.py -h to obtain the available command line arguments

deredden.py --> Script to deredden flux line maps, it needs as input the balmer decrement map (Halpha/Hbeta) and it will automatically localize the flux maps

image_stats.py --> Script to obtain some stats from a particular region of an image (python image_stats.py image.fits -region ds9reigonfile.reg

read_mappings.py --> Script to obtain predicted line ratios from the mapping libraries (python read_mappings.py <files> (for instance V_*_b0_001_s_lines.txt V_*b1_s_lines.txt [MVQP]_*b[0e]_s_lines.txt [MVQ]_*b10_s_lines.txt T_*b0_001_s_lines.txt T_*b1_s_* T_*b10_s_lines.txt). The script will read the paths to the BPT diagrams from bpt_config.py file
