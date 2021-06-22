# muse_project
Scripts to deal with some subtleties regarding MUSE data

bpt_newer.py -> bpt diagram from Law et al. 2021, now with geometry

marxsim_script.py --> Script to be run from the command line to simulate N PSFs based on a given chandra observation. Ciao and marx need to be properly setup before running the script. Use python marxsim_script.py -h to obtain the available command line arguments

check_extension.py --> Script to be run after marxsim_script.py has been run to analyse the output. Use python marxsim_script.py -h to obtain the available command line arguments

deredden.py --> Script to deredden flux line maps, it needs as input the balmer decrement map (Halpha/Hbeta) and it will automatically localize the flux maps

image_stats.py --> Script to obtain some stats from a particular region of an image (python image_stats.py image.fits -region ds9reigonfile.reg
