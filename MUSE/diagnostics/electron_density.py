import pyneb as pn
from mpdaf.obj import Image
import numpy as np
import argparse
import sys, time

def calculate_electron_density(sii_ratio_map, T, ext=1):
    """
    Calculate electron density map using PyNeb diagnostics.
    
    Parameters:
    -----------
    sii_ratio_map : str or mpdaf.obj.Image
        Path to the SII ratio map file or MPDAF Image object
    T : float or str
        Electron temperature in K or path to FITS file containing SIII 6312/9069 ratio map

    Returns:
    --------
    ne_map : mpdaf.obj.Image
        Electron density map
    tem_map : mpdaf.obj.Image (optional, if SIII map provided)
        Electron temperature map
    """
    
    # Load SII ratio map if it's a file path
    ratio_map = Image(sii_ratio_map, ext)
    
    # Get the ratio data
    ratio_data = ratio_map.data
    
    if isinstance(T, float):
        # Single temperature case
        # Initialize output array
        ne_data = np.full(ratio_data.shape, np.nan, dtype=np.float64)
        # Handle NaN and invalid values
        valid_mask = np.isfinite(ratio_data) & (ratio_data > 0) & ~np.isnan(ratio_data)
        # Flatten arrays for easier processing
        ratio_flat = ratio_data[valid_mask]

        # Calculate electron density for each pixel
        print(f"Calculating electron density for {np.count_nonzero(valid_mask)} pixels... This may take a while...")
        start = time.time()
        # PyNeb expects ratio as 6716/6731
        s2 = pn.Atom("S", 2)
        ne_data[valid_mask] = s2.getTemDen(int_ratio=ratio_flat, tem=T, wave1=6716, wave2=6731)
        end = time.time()
        print(f"Calculation completed in {end - start:.2f} seconds.")
        
        # Create output image with same properties as input
        ne_map = ratio_map.copy()
        ne_map.data = ne_data
        
        # Update header information
        ne_map.primary_header['BUNIT'] = ('cm-3', 'Electron density')
        ne_map.primary_header['HISTORY'] = f'Electron density calculated from SII ratio at T={T}K'
        
        return ne_map

    else:
        # SIII diagnostic case - simultaneous Te and Ne calculation
        print("Loading SIII ratio map for simultaneous Te/Ne calculation...")
        siii_map = Image(T, ext)
        siii_data = siii_map.data
        
        # Initialize diagnostics
        diags = pn.Diagnostics()
        
        # Initialize output arrays
        ne_data = np.full(ratio_data.shape, np.nan, dtype=np.float64)
        tem_data = np.full(ratio_data.shape, np.nan, dtype=np.float64)
        
        # Find valid pixels for both diagnostics
        validratio = np.isfinite(ratio_data) & (ratio_data > 0)
        valids3 = np.isfinite(siii_data) & (siii_data > 0)
        valid_mask = validratio & valids3
        
        # Get valid data
        s2_flat = ratio_data[valid_mask]
        siii_flat = siii_data[valid_mask]
        
        print(f"Processing {len(s2_flat)} valid pixels for simultaneous Te/Ne calculation...")
        start = time.time()

        # Solve for Te and Ne
        # be careful here cause the diganostic in pyneb is defined the other way around, so we need to invert the SII ratio
        tem_flat_results, ne_flat_results = diags.getCrossTemDen(diag_den='[SII] 6731/6716', diag_tem='[SIII] 6312/9069', 
                                                                 value_tem=siii_flat, value_den=1 / s2_flat)
        
        # Map results back to 2D arrays
        tem_data[valid_mask] = tem_flat_results
        ne_data[valid_mask] = ne_flat_results

        # assume 10,000 K for invalid SIII ratios
        T_assume = 10000
        nonsiiidata = ~valids3 & validratio

        print(f"Processing remaining {np.count_nonzero(nonsiiidata)} of unavailable [SIII] information assuming T = {T_assume}")
        tem_data[nonsiiidata] = T_assume
        s2 = pn.Atom("S", 2)
        ne_data[nonsiiidata] = s2.getTemDen(int_ratio=ratio_data[nonsiiidata], tem=T_assume, wave1=6716, wave2=6731)

        end = time.time()
        print(f"Calculation completed in {end - start:.2f} seconds.")
        
        # Create output images
        ne_map = ratio_map.copy()
        ne_map.data = ne_data
        ne_map.primary_header['BUNIT'] = ('cm-3', 'Electron density')
        ne_map.primary_header['HISTORY'] = 'Electron density from simultaneous SII+SIII diagnostics'
        
        tem_map = ratio_map.copy()
        tem_map.data = tem_data
        tem_map.primary_header['BUNIT'] = ('K', 'Electron temperature')
        tem_map.primary_header['HISTORY'] = 'Electron temperature from simultaneous SII+SIII diagnostics'
        
        return ne_map, tem_map

def main():
    """Main function to handle command line arguments and run the calculation."""
    
    parser = argparse.ArgumentParser(description='Calculate electron density map from SII ratio map using PyNeb diagnostics')    
    parser.add_argument('SII6716_SII6731_ratio', type=str, help='Path to the SII 6713/6731 ratio map FITS file')    
    parser.add_argument('-T', "--temperature", help='Electron temperature in K or path to SIII 6312/9069 ratio FITS file', default=10000)   
    parser.add_argument('--ext', default=1, help='Extension of the fits file to read', type=int)
    parser.add_argument('-o', '--output', default='map', help='Output rootname for the electron density map (default: map)')
    args = parser.parse_args()
    
    try:
        # Determine if temperature is a file or value
        try:
            T_value = float(args.temperature)
            print(f"Using fixed temperature: {T_value} K")
            is_file = False
        except ValueError:
            print(f"Using SIII ratio map: {args.temperature}")
            T_value = args.temperature
            is_file = True
        
        # Calculate electron density
        print(f"Loading SII ratio map: {args.SII6716_SII6731_ratio}")
        
        result = calculate_electron_density(args.SII6716_SII6731_ratio, T_value, args.ext)
        
        if is_file:
            # Simultaneous calculation - two outputs
            ne_map, tem_map = result
            
            # Save both maps
            ne_output = f"ne_{args.output}.fits"
            tem_output = f"tem_{args.output}.fits"
            
            print(f"Saving electron density map to: {ne_output}")
            ne_map.write(ne_output)
            
            print(f"Saving electron temperature map to: {tem_output}")
            tem_map.write(tem_output)
            
            # Show statistics
            for data, name, unit in [(ne_map.data, "Density", "cm⁻³"), (tem_map.data, "Temperature", "K")]:
                valid_data = data[np.isfinite(data)]
                if len(valid_data) > 0:
                    print(f"\nElectron {name} Statistics:")
                    print(f"  Min {name.lower()}: {np.min(valid_data):.2e} {unit}")
                    print(f"  Max {name.lower()}: {np.max(valid_data):.2e} {unit}")
                    print(f"  Mean {name.lower()}: {np.mean(valid_data):.2e} {unit}")
                    print(f"  Median {name.lower()}: {np.median(valid_data):.2e} {unit}")
                    print(f"  Valid pixels: {len(valid_data)}")
        else:
            # Single temperature calculation - one output
            ne_map = result
            output = f"ne_{args.output}_T{T_value}.fits"
            print(f"Saving electron density map to: {output}")
            ne_map.write(output)
            
            # Show statistics
            valid_data = ne_map.data[np.isfinite(ne_map.data)]
            if len(valid_data) > 0:
                print("\nElectron Density Statistics:")
                print(f"  Min density: {np.min(valid_data):.2e} cm⁻³")
                print(f"  Max density: {np.max(valid_data):.2e} cm⁻³")
                print(f"  Mean density: {np.mean(valid_data):.2e} cm⁻³")
                print(f"  Median density: {np.median(valid_data):.2e} cm⁻³")
                print(f"  Valid pixels: {len(valid_data)}")
            else:
                print("No valid density values calculated.")
        
        print("Done!")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()