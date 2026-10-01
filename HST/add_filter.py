#!/usr/bin/env python3
# @Author: Andrés Gúrpide <agurpide>
# @Date:   05-08-2025
# @Email:  agurpidelash@soton.ac.uk
# Script to add HST filter throughput to MPDAF filter list

import argparse
import os
import logging
from astropy.io import fits
from astropy.time import Time
import stsynphot as stsyn
import numpy as np

def get_band(hst_hdul):
    """
    Extract filter band information from HST FITS header
    
    Parameters:
    -----------
    hst_hdul : astropy.io.fits.HDUList
        Opened HST FITS file
        
    Returns:
    --------
    stsyn.SpectralElement
        Filter throughput band
    """
    date = hst_hdul[0].header["DATE-OBS"]
    obsdate = Time(date).mjd
    
    if "PHOTMODE" in hst_hdul[1].header:
        obs_mode = hst_hdul[1].header["PHOTMODE"]
    else:
        obs_mode = hst_hdul[0].header["PHOTMODE"]
    
    if "WFC3" in obs_mode or "ACS SBC" in obs_mode:
        keywords = obs_mode.split(" ")
        obsmode = "%s,%s,%s,mjd#%.2f" % (keywords[0], keywords[1], keywords[2], obsdate)
    elif "WFPC2" in obs_mode:
        keywords = obs_mode.split(",")
        # ['WFPC2', '1', 'A2D7', 'F656N', '', 'CAL']
        obsmode = "%s,%s,%s,%s,%s" % (keywords[0], keywords[1], keywords[2], keywords[3], keywords[5])
    elif "WFC1" in obs_mode or "HRC" in obs_mode:
        keywords = obs_mode.split(" ")
        obsmode = "%s,%s,%s,%s" % (keywords[0], keywords[1], keywords[2], keywords[3])
    

    return stsyn.band(obsmode)

def get_filter_name(hst_hdul):
    """
    Extract filter name from HST FITS header
    
    Parameters:
    -----------
    hst_hdul : astropy.io.fits.HDUList
        Opened HST FITS file
        
    Returns:
    --------
    str
        Filter name in format INSTRUMENT_FILTER
    """
    try:
        instrument = hst_hdul[0].header["INSTRUME"]
    except KeyError:
        raise KeyError(f"INSTRUME keyword not found in FITS header {hst_hdul.filename()}")
    if "FILTER" in hst_hdul[0].header:
        hst_filter = hst_hdul[0].header["FILTER"]
    elif "FILTNAM1" in hst_hdul[0].header:
        hst_filter = hst_hdul[0].header["FILTNAM1"]
    else:
        hst_filter = hst_hdul[0].header["FILTER1"]
    
    return f"{instrument}_{hst_filter}"

def add_filter_to_fits(filter_name, band, filter_fits_path):
    """
    Add filter throughput to MPDAF filter list FITS file
    
    Parameters:
    -----------
    filter_name : str
        Name of the filter (e.g., WFPC2_F555W)
    band : stsyn.SpectralElement
        Filter throughput band
    filter_fits_path : str
        Path to the filter_list.fits file
    """
    # Sample the filter throughput
    # Use wavelength range that covers typical HST filters (2000-11000 Angstroms)
    wavelengths = band.waveset
    throughput = band(wavelengths)
    wavelengths = wavelengths[throughput > 0]  # Filter out zero throughput values
    throughput = throughput[throughput > 0]  # Filter out zero throughput values

    # Create new HDU with filter data
    col1 = fits.Column(name='lambda', format='E', array=wavelengths.value, unit='ANGSTROM')
    col2 = fits.Column(name='throughput', format='E', array=throughput.value)
    
    new_hdu = fits.BinTableHDU.from_columns([col1, col2], name=filter_name)
    new_hdu.header['EXTNAME'] = filter_name
    new_hdu.header['COMMENT'] = f'Filter throughput for {filter_name}'
    
    # Read existing FITS file or create new one
    if os.path.exists(filter_fits_path):
        with fits.open(filter_fits_path, mode='update') as hdul:
            # Check if filter already exists
            existing_names = [hdu.name for hdu in hdul if hasattr(hdu, 'name')]
            if filter_name in existing_names:
                logger.warning(f"Filter {filter_name} already exists in {filter_fits_path}, skipping...")
                return False
            
            # Append new HDU
            hdul.append(new_hdu)
            hdul.flush()
            logger.info(f"Successfully added filter {filter_name} to {filter_fits_path}")
            return True
    else:
        raise FileNotFoundError(f"Filter file {filter_fits_path} does not exist.")

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description='Add HST filter throughput to MPDAF filter list')
    parser.add_argument('hst_images', nargs='+', help='Path(s) to HST FITS image(s)')
    parser.add_argument('--filter-fits', 
                       default='/home/andresgur/anaconda3/envs/muse/lib/python3.13/site-packages/mpdaf/obj/filters/filter_list.fits',
                       help='Path to MPDAF filter_list.fits file (default: /home/andresgur/anaconda3/envs/muse/lib/python3.13/site-packages/mpdaf/obj/filters/filter_list.fits)')
    
    args = parser.parse_args()
    
    # Setup logging
    global logger
    scriptname = os.path.basename(__file__)
    logger = logging.getLogger(scriptname)
    logger.setLevel(logging.DEBUG)
    out_format = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(out_format)
    logger.addHandler(stream_handler)
    
    # Check if HST images exist
    missing_files = []
    for hst_image in args.hst_images:
        if not os.path.isfile(hst_image):
            missing_files.append(hst_image)
    
    if missing_files:
        logger.error(f"HST image(s) not found: {', '.join(missing_files)}")
        return 1
    
    processed_filters = []
    failed_filters = []
    
    # Process each HST image
    for hst_image in args.hst_images:
        try:
            # Open HST image
            logger.info(f"Processing HST image {hst_image}")
            with fits.open(hst_image) as hdul:
                # Check if it's an HST image
                if hdul[0].header.get("TELESCOP") != "HST":
                    logger.error(f"Image {hst_image} is not from HST telescope, skipping...")
                    failed_filters.append(hst_image)
                    continue
                
                # Get filter information
                filter_name = get_filter_name(hdul)
                logger.info(f"Detected filter: {filter_name}")
                
                # Get filter throughput band
                logger.info("Extracting filter throughput...")
                band = get_band(hdul)
                
                # Add filter to FITS file
                success = add_filter_to_fits(filter_name, band, args.filter_fits)
                
                if success:
                    logger.info(f"Filter {filter_name} successfully added to MPDAF filter list")
                    processed_filters.append(filter_name)
                else:
                    failed_filters.append(hst_image)
                    
        except Exception as e:
            logger.error(f"Error processing HST image {hst_image}: {str(e)}")
            failed_filters.append(hst_image)
    
    # Print summary
    logger.info(f"\nProcessing Summary:")
    logger.info(f"Successfully processed {len(processed_filters)} filters: {', '.join(processed_filters)}")
    if failed_filters:
        logger.warning(f"Failed to process {len(failed_filters)} images: {', '.join(failed_filters)}")
        return 1
    else:
        return 0

if __name__ == "__main__":
    exit(main())
