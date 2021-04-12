##script to extract the radial profile of a certain source and check if it is emitted

ra=1058.9441
dec=857.75425
event_image=MAST_2019-08-14T0354_656N.fits
nannuli=10
innersource=2
outersource=15
background=15
innerback=15
outerback=25
source=ixo27

printf 'Creating source annuli'

step=$(echo "scale=2; ($outersource - $innersource)/$nannuli" | bc)
start=$innersource
echo $step
echo $start
for (( c=$innersource; c<=$outersource; c++ )); do
	end=$(echo "scale=2; ($start + $step)" | bc)
	echo "annulus($ra,$dec,$innersource,$end)" >> annuli_source.reg
	start=$end
done

echo "annulus($ra,$dec,$innerback,$outerback)" >> annuli_back.reg

echo "Check if you are happy with the regions. Otherwise modify them and resave them"

read -p "Done? Enter y to continue y"
if [[ $REPLY =~ ^[Yy]$ ]]; then
	infile="$event_image[bin sky=@annuli_source.reg]"
	outfile=$source"r_profile.fits"
	bkg="$event_image[bin sky=@annuli_back.reg]"
	punlearn dmextract
 	dmextract "$infile" $outfile bkg="$bkg" opt=generic clobber=yes
	python_script="from pycrates import read_file; import matplotlib.pyplot as plt; tab = read_file('r_profile.fits'); xx = tab.get_column('rmid').values; yy = tab.get_column('sur_bri').values; ye = tab.get_column('sur_bri_err').values; plt.errorbar(xx,yy,yerr=ye, marker='o');plt.xscale('log');plt.yscale('log'); plt.xlabel('Radius (pixel)');plt.ylabel('Surface brightness (photons/cm**2/pixel**2/s)'); plt.show()"
	python -c "${python_script}"

# TODO: create simulate PSF with MARX
fi






