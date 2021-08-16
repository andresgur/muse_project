# Prior to execution a file named detected_sources_broadband.csv has to be created if you want to see the sources on ds9 -> to create it, open with fv detected_sources_broadband.fits (the output sources from wavdetect) --> select all --> export as csv


# Source match radius in arcseconds
radius=10 # approx distance from a reference HST compared to chandra sources found experimentally, default is 12


#Source positions from the two source files can only match each other if they are no farther apart than the value specified as the source match radius. Only a single reference source can be within the source match radius of a input source, and only a single input source can be within the source match radius of a reference source. If multiple sources are found in either case, the source pair is not used to calculate the transform. 

residtype=0 # 0 --> upper limit source error ratio above residfac eliminated 1 -- > upper limit on the average
ds9_path=$HOME 

# event file or image of the FOV to create the catalog of external sources with ds9 (it can also be the event file itself from the repro dir)
fov_file=sourceimage_broadband.fits

res_radius=3

if [ $# -ne 1 ] ; then
	param1="<ref_catalog>"
	echo "Provide reference catalog for coordinate update. Usage $param1. This script has to be launched from detectedsources dir. E.g. 'gaia' '2mass'. See https://cxc.cfa.harvard.edu/ciao/ahelp/wcs_match.html#plist.refsrcfile on how to provide a custom reference catalog"
	exit
fi


printf "Using reference catalog: \n \t $1  \n"

if [ $1 == "gaia" ] ; then
	$ds9_path/ds9 $fov_file -scale log -cmap heat -catalog $1 -catalog maxrows 2000 -catalog filter "abs(\$pmRA) < 20 && abs(\$pmDE) <20 && \$e_RA_ICRS < 20 && \$e_DE_ICRS < 20 && \$Dup==0" -catalog export csv $1_chandra.csv -catalog import csv detected_sources_broadband.csv -catalog symbol color green -catalog match error $radius arcsec &
	# convert units https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html
	dmtcalc "$1_chandra.csv[opt skip=1]" $1_catalog.fits expr="ra_err=col4/1000/3600" clobber=yes
	dmtcalc $1_catalog.fits gaia_catalog.fits expr="dec_err=col6/1000/3600" clobber=yes
	dmtcalc $1_catalog.fits gaia_catalog.fits expr="dec=col5" clobber=yes
	dmtcalc $1_catalog.fits gaia_catalog.fits expr="ra=col3" clobber=yes
	ref_cat="gaia_catalog.fits"
outname=$1
elif [ $1 == "2mass" ];then

	$ds9_path/ds9 $fov_file -scale log -cmap heat -catalog $1 -catalog maxrows 2000 -catalog export csv $1_chandra.csv -catalog import csv detected_sources_broadband.csv -catalog symbol color green -catalog match error $radius arcsec &
	# there are no uncertainties as far as I can tell
	ref_cat="$1_chandra.csv[opt skip=1][cols ra=col1,dec=col2]"
	outname=$1
else
	ref_cat=$1
	outname=$(basename -- "$1")
	$ds9_path/ds9 $fov_file -scale log -cmap heat -catalog import csv $1 -catalog maxrows 2000 -catalog export csv $outname\_chandra.csv -catalog import csv detected_sources_broadband.csv -catalog symbol color green -catalog match error $radius arcsec &

fi

match_out=$outname\_chandra_match.fits

chandra_ev=$(ls ../repro/*repro_evt2.fits | head -1) 
 
out_dir=../corr_corrected_$outname

mkdir $out_dir

chandra_asol=$(ls ../repro/*asol*.fits | head -1)

ev_out=astrocorr_evt.fits

asol_out=astrocorr_asol1.fits

mask_out=astrocorr_mask1.fits

dmcopy $chandra_asol $out_dir/$asol_out op=all clobber=yes

dmcopy $chandra_ev $out_dir/$ev_out op=all clobber=yes

chandra_sources=$(ls detected_sources_broadband.fits | head -1) 

echo "Chandra sources: $chandra_sources"

corrected_out=$out_dir/$ev_out

# match catalogs 
printf "Matching catalogs $chandra_sources and reference $ref_cat...\n"
punlearn wcs_match
# use any not corrected file for wcs but be consistent for wcs_match and wcs_update
wcs_match "$chandra_sources[exclude RA_ERR=0, DEC_ERR=0]" "$ref_cat" "$out_dir/$match_out" wcsfile=$chandra_ev radius=$radius clobber=yes residlim=$res_radius  verbose=2 logfile="$out_dir/match.log" method='rst' residtype=$residtype residfac=2
# method rst for everything, or trans for only transaltion
cat $out_dir/'match.log'

read -p "Are you happy with the match? y/n"
if [[ $REPLY =~ ^[Yy]$ ]]; then
 
echo "Updating asol file...$out_dir/$asol_out"

wcs_update infile=$chandra_asol transformfile=$out_dir/$match_out outfile=$out_dir/$asol_out clobber=yes verbose=1 wcsfile=$chandra_ev

echo "Updating wcs header in $corrected_out"

wcs_update infile=$corrected_out transformfile=$out_dir/$match_out verbose=2 clobber=yes outfile="" wcsfile=$chandra_ev

echo "Updating wcs header for detect sources output"

declare -a arr=("soft" "medium" "hard" "broadband")

for band in "${arr[@]}" ; do 
	# update region files using ds9 (seems not to be working anymore :()
	ds9 $corrected_out -regions sources_$band.reg -regions system wcs -regions save $out_dir/sources_$band\_fk5.reg -exit
		
	for file in $(ls *$band*\.fits)
		 do echo "Updating WCS header $file"
		 dmcopy $file $out_dir/$file op=all clobber=yes
		 wcs_update infile=$out_dir/$file transformfile=$out_dir/$match_out verbose=2 clobber=yes outfile="" wcsfile=$chandra_ev
		 done
done


# update event file; asol keyword with the corrected asol file; mask file
dmhedit $corrected_out file= op=add key=ASOLFILE value=$asol_out

fi

