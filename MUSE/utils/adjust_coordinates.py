# @Author: Andrés Gúrpide <agurpide>
# @Date:   10-06-2019
# @Email:  agurpidelash@soton.ac.uk
# @Last modified by:   agurpide
# @Last modified time: 10-12-2021
# Script to adjust the coordinates of a cube from an input image by cross-correlating it( Preferabley HST)


# imports
import argparse
from mpdaf.obj import Cube, Image
import os
import logging
from astropy.io import fits
import gc
from astroquery.simbad import Simbad
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from math import cos, pi
from astropy.table import Column
import numpy as np

def search_sources(coords, radius, obs_epoch=Time('J2000'), equinox=2000, catalog='gaia'):
    if catalog == 'simbad':
        my_simbad = Simbad()
        # add the object type to the output results
        my_simbad.add_votable_fields('otype')
        result_table = my_simbad.query_region(coords,
                                              radius=radius, epoch=obs_epoch.jyear_str, equinox=equinox)
    elif catalog == 'gaia':
        j = Gaia.cone_search(coordinate=coords, radius=radius)
        result_table = j.get_results()
        pm_ra_cosdec = (result_table['pmra'] * np.cos(result_table['dec'].to("rad"))) * u.mas / u.yr
        # write coordinates of the sources propagated to the new epoch
        c = SkyCoord(ra=result_table['ra'], dec=result_table['dec'],
                     pm_ra_cosdec=pm_ra_cosdec, pm_dec=result_table['pmdec'],
                     radial_velocity=result_table["radial_velocity"],
                     obstime=Time(result_table['ref_epoch'].value, format='decimalyear'))

        # get position of the stars in the observation position
        c.apply_space_motion(new_obstime=obs_epoch)
        result_table.add_column(Column(data=c.ra.value, format='E', name='ra_%s' % obs_epoch.jyear_str))
        result_table.add_column(Column(data=c.dec.value, format='E', name='dec_%s' % obs_epoch.jyear_str))
    return result_table

# read arguments
ap = argparse.ArgumentParser(description='Adjust coordinates from HST image')
ap.add_argument("input_cube", nargs=1, help="Cube to be corrected", type=str)
ap.add_argument("--hst", nargs=1, help="HST image", type=str)
ap.add_argument("-n", "--nbright", nargs="?", help="Number of bright stars to use in the astrometric error calculation. Default 5", type=int, default=5)
ap.add_argument("-r", "--radius", nargs="?", help="Radius (in minutes) to look for stars for the astromic error calculation. Default 0.5", type=float, default=0.5)
ap.add_argument("-c", "--catalog", nargs="?", help="Catalog for the calculation of the astrometric error", choices=['gaia', 'simbad'], default='gaia')
ap.add_argument("-o", "--outdir", nargs='?', type=str, help="Output directory", default="coordadjusted")

# parse args
args = ap.parse_args()
nsigma = 3

# create output dir
if not os.path.isdir(args.outdir):
    os.mkdir(args.outdir)

# logger
scriptname = os.path.basename(__file__)

logger = logging.getLogger(scriptname)
logger.setLevel(logging.DEBUG)
out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')

# handler
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.DEBUG)
stream_handler.setFormatter(out_format)
logger.addHandler(stream_handler)


input_cube = args.input_cube[0]
hst_image = args.hst[0]
outcube = '%s/corr%s' % (args.outdir, input_cube)

if os.path.isfile(input_cube) and os.path.isfile(hst_image):
    logger.info("Loading cube %s" % input_cube)
    logger.info("Loading succesful")
    hdul = fits.open(hst_image)
    if hdul[0].header["TELESCOP"] == "HST":
        # http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=HST/ACS_WFC.F555W&&mode=browse&gname=HST&gname2=ACS_WFC#filter
        instrument = hdul[0].header["INSTRUME"]
        if "FILTER" in hdul[0].header:
            hst_filter = hdul[0].header["FILTER"]
        elif "FILTNAM1" in hdul[0].header:
            hst_filter = hdul[0].header["FILTNAM1"]
        else:
            hst_filter = hdul[0].header["FILTER1"]
        filter_key = "%s_%s" % (instrument, hst_filter)
        logger.info("Constructing image in the %s filter" % filter_key)
        # free memory
        hdul.close()      

        cube = Cube(input_cube)

        try:
            cube_image = cube.get_band_image(filter_key)
        except ValueError:
            print(f"Filter {filter_key} not found, attempting to use ACS filter version instead")
            filter_key = filter_key.replace("WFPC2", "ACS")
            try:
                cube_image = cube.get_band_image(filter_key)
            except ValueError:
                print(f"Filter {filter_key} not found, attempting to use WFC3 filter version instead")
                filter_key = filter_key.replace("ACS", "WFC3")
            try:
                cube_image = cube.get_band_image(filter_key)
            except ValueError:
                print(f"Filter {filter_key} not found, reverting to white light image")
                cube_image = cube.sum(axis=0)
                filter_key = "white_light"
        print(f"Successfully created {filter_key} image")
            
    else:
        cube = Cube(input_cube)
        filter_key = "white_light"
        print("Image format not HST, using cube white light image for the coordinate correction")
        cube_image = cube.sum(axis=0)
    del cube
    wlname = "%s/%s" % (args.outdir, input_cube.replace(".fits", "_imgcorr_%s.fits" % filter_key))
    # free memory

    hst_im = Image(hst_image)
    logger.info("Estimating coordinate offset...")
    dy, dx = cube_image.estimate_coordinate_offset(hst_im, nsigma=nsigma)
    cube_image.adjust_coordinates(hst_im, nsigma=nsigma, inplace=True)
    cube_image.write(wlname)
    # update cube WCS header
    cube = Cube(input_cube)
    cube.wcs.set_crpix1(cube_image.wcs.get_crpix1())
    cube.wcs.set_crpix2(cube_image.wcs.get_crpix2())
    cube.wcs.set_cd(cube_image.wcs.get_cd())
    pixel_scale = cube.primary_header['HIERARCH ESO OCS IPS PIXSCALE']
    cube.primary_header["ASTRO_CORR_PX"] = ("%.2f, %.2f" % (dy, dx), 'Coordinate shift dy, dx in pixels')
    print(dy * pixel_scale)
    cube.primary_header["ASTRO_CORR_AS"] = ("%.2f, %.2f" % ((-dy * pixel_scale), (-dx * pixel_scale)), 'Coordinate shift dy, dx in arcseconds')
    cube.primary_header['COMMENT'] = "Coordinates corrected using %s image" % hst_image
    cube.write(outcube)
    del cube
    logger.info("Output successfully stored to %s" % args.outdir)
    wfile = fits.open(wlname)
    header = wfile[0].header
    radesyskeyword = "RADECSYS" if "RADECSYS" in header else "RADESYS"
    radecsys = (header[radesyskeyword].lower())
    ra_center = header["RA"]
    dec_center = header["DEC"]
    obs_epoch = Time(header['DATE-OBS'])
    equinox = header['EQUINOX']
    wfile.close()
    radius = args.radius  /60./ 180. * pi * u.rad
    coords = SkyCoord(ra=ra_center, dec=dec_center, unit=("deg", "deg"), frame=radecsys)
    result_table = search_sources(coords, radius, obs_epoch, equinox, args.catalog)
    # sort table by mag and get the five brightests stars
    result_table.sort("phot_g_mean_mag")
    result_table = result_table[:args.nbright]

    # calculate error in shift
    FWHM = header["SKY_RES"]
    errra = []
    errdec = []
    star_data = []
    
    # Create output file for star analysis
    output_file = os.path.join(args.outdir, f"astro_metric_error_{filter_key}.txt")
    
    for row in result_table:
        starcenter = (row["dec"], row["ra"])
        
        imgfit = cube_image.subimage(starcenter, size=2 * FWHM, unit_size="arcsec", unit_center="deg")
        print(f"\nFitting star {row['source_id']} at RA: {row['ra']:.5f} Dec: {row['dec']:.5f}\n")
        print("MUSE")
        muse_fit = imgfit.gauss_fit(center=starcenter, flux=np.max(imgfit.data), fwhm=FWHM, circular=True,
                                    cont=0, fit_back=True, peak=True, factor=1, weight=True, verbose=True, full_output=False)
        # for the HST we do not need that big of an image
        imgfit = hst_im.subimage(starcenter, size=1.2, unit_size="arcsec", unit_center="deg", minsize=0.1)
        print("HST")
        hst_fit = imgfit.gauss_fit(center=starcenter, flux=np.max(imgfit.data), fwhm=0.11, circular=True,
                                    cont=0, fit_back=True, peak=True, factor=1, weight=True, verbose=True, full_output=False)
        if hst_fit.flux<=0 or muse_fit.flux<=0:
            logger.warning("Star %s either not found or too crowded in one of the images, skipping..." % row["source_id"])
            continue
        errors = np.abs(muse_fit.center - hst_fit.center) 
        mean_dec = (muse_fit.center[0] + hst_fit.center[0]) / 2.0
        # account for the widening of the error in RA due to the declination (this is valid only for small offsets!)
        errra.append(abs(errors[1] * cos(mean_dec * np.pi / 180.)))
        errdec.append(errors[0])
        
        # Collect star data for output
        star_info = {
            'star_id': row["source_id"],
            'ra': row["ra"],
            'dec': row["dec"],
            'gmag': row["phot_g_mean_mag"],
            'muse_ra': muse_fit.center[1],
            'muse_dec': muse_fit.center[0],
            'hst_ra': hst_fit.center[1],
            'hst_dec': hst_fit.center[0],
            'err_ra': errors[1],
            'err_dec': errors[0]
        }
        star_data.append(star_info)
    
    # Write star analysis to file
    with open(output_file, 'w') as f:
        # Write header
        f.write("star\tra\tdec\tGmag\tmuse_ra\tmuse_dec\thst_ra\thst_dec\terr_ra\terr_dec\n")
        
        # Write data for each star
        for star in star_data:
            f.write(f"{star['star_id']}\t{star['ra']:.6f}\t{star['dec']:.6f}\t{star['gmag']:.3f}\t"
                   f"{star['muse_ra']:.6f}\t{star['muse_dec']:.6f}\t{star['hst_ra']:.6f}\t{star['hst_dec']:.6f}\t"
                   f"{star['err_ra']:.6f}\t{star['err_dec']:.6f}\n")
            
        f.write(f"#Errors in arcsec\n#dRa: {np.mean(errra) * 3600.:.4f}\n#dDec {np.mean(errdec) * 3600.:.4f}")
    
    logger.info(f"Star analysis saved to {output_file}")
    
    # Print summary to console
    print("\nStar Analysis Summary:")
    print("star\tra\tdec\tGmag\tmuse_ra\tmuse_dec\thst_ra\thst_dec\terr_ra\terr_dec")
    for star in star_data:
        print(f"{star['star_id']}\t{star['ra']:.6f}\t{star['dec']:.6f}\t{star['gmag']:.3f}\t"
              f"{star['muse_ra']:.6f}\t{star['muse_dec']:.6f}\t{star['hst_ra']:.6f}\t{star['hst_dec']:.6f}\t"
              f"{star['err_ra']:.6f}\t{star['err_dec']:.6f}")

    


else:
    logger.error("Input cube %s or HST image %s not found" % (input_cube, hst_image))
 