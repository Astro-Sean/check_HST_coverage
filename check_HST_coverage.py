#!/usr/bin/env python3
"""
Check for public HST images at a given RA/DEC or TNS name
Usage: python check_hst_coverage.py --ra 123.456 --dec 45.678
       python check_hst_coverage.py --tns AT2023xyz
       python check_hst_coverage.py --ra 123.456 --dec 45.678 --download --plot
"""

import argparse
import sys
import os
import matplotlib.pyplot as plt
from astroquery.mast import Observations
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits

def resolve_tns_name(tns_name):
    """Resolve TNS name to coordinates using SIMBAD."""
    try:
        from astroquery.simbad import Simbad
        # Reset to default fields
        Simbad.reset_votable_fields()
        result = Simbad.query_object(tns_name)
        if result is not None and len(result) > 0:
            # Coordinates are in lowercase ra/dec columns (ICRS, degrees by default)
            ra = float(result['ra'][0])
            dec = float(result['dec'][0])
            print(f"Resolved TNS name '{tns_name}' to RA={ra:.6f}, DEC={dec:.6f}")
            return ra, dec
        else:
            print(f"Could not resolve TNS name '{tns_name}'")
            return None, None
    except Exception as e:
        print(f"Error resolving TNS name: {e}")
        return None, None

def point_in_polygon(ra, dec, s_region):
    """Check if a point (ra, dec) is within the polygon defined by s_region."""
    try:
        # Parse s_region (format: "POLYGON ICRS ra1 dec1 ra2 dec2 ...")
        if not s_region or s_region == '--':
            return False
        
        parts = s_region.split()
        if parts[0] != 'POLYGON' or len(parts) < 4:
            return False
        
        # Extract coordinates (skip "POLYGON" and coordinate system like "ICRS")
        start_idx = 2 if len(parts) > 2 and parts[1].upper() in ['ICRS', 'J2000'] else 1
        coords = []
        for i in range(start_idx, len(parts), 2):
            if i + 1 < len(parts):
                try:
                    ra_coord = float(parts[i])
                    dec_coord = float(parts[i + 1])
                    coords.append((ra_coord, dec_coord))
                except (ValueError, IndexError):
                    continue
        
        if len(coords) < 3:
            return False
        
        # Ray casting algorithm for point-in-polygon
        n = len(coords)
        inside = False
        
        x, y = ra, dec
        p1x, p1y = coords[0]
        
        for i in range(n + 1):
            p2x, p2y = coords[i % n]
            
            if y > min(p1y, p2y):
                if y <= max(p1y, p2y):
                    if x <= max(p1x, p2x):
                        if p1y != p2y:
                            xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                        if p1x == p2x or x <= xinters:
                            inside = not inside
            p1x, p1y = p2x, p2y
        
        return inside
    except Exception as e:
        print(f"Error parsing s_region: {e}")
        return False

def check_hst_coverage(ra, dec, radius=0.1):
    """Check for HST observations at given coordinates."""
    print(f"\nSearching for HST observations at RA={ra:.6f}, DEC={dec:.6f}")
    print(f"Search radius: {radius} deg")
    print("="*60)
    
    # Query HST observations at this position
    coord = SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    obs_table = Observations.query_region(coord, radius=f"{radius} deg")
    
    # Filter for HST only
    hst_obs = obs_table[obs_table['obs_collection'] == 'HST']
    
    print(f"\nTotal observations in region: {len(obs_table)}")
    print(f"HST observations: {len(hst_obs)}")
    
    if len(hst_obs) == 0:
        print("\nNo HST observations found at this location.")
        return None
    
    # Convert to pandas for easier handling
    df = hst_obs.to_pandas()
    
    # Filter for observations where target is within the image footprint
    print("\nFiltering for observations that overlap with target coordinates...")
    df['overlaps'] = df.apply(lambda row: point_in_polygon(ra, dec, row.get('s_region', '--')), axis=1)
    df_overlap = df[df['overlaps'] == True]
    
    print(f"Observations overlapping target: {len(df_overlap)}")
    
    if len(df_overlap) == 0:
        print("\nNo HST observations overlap with target coordinates.")
        print("Showing all observations in search radius instead:")
        df_overlap = df
    else:
        df = df_overlap
    
    # Display relevant information
    print(f"\n{'='*80}")
    print(f"{'Instrument':<15} {'Filter':<20} {'Obs ID':<20} {'RA':<12} {'DEC':<12}")
    print(f"{'='*80}")
    
    for _, row in df.iterrows():
        instrument = row.get('instrument_name', 'N/A')
        filters = row.get('filters', 'N/A')
        obs_id = row.get('obs_id', 'N/A')
        s_ra = row.get('s_ra', ra)
        s_dec = row.get('s_dec', dec)
        
        print(f"{instrument:<15} {filters:<20} {obs_id:<20} {s_ra:<12.6f} {s_dec:<12.6f}")
    
    print(f"{'='*80}")
    
    # Summary
    print(f"\nSummary:")
    print(f"  Unique instruments: {df['instrument_name'].nunique()}")
    print(f"  Unique filters: {df['filters'].nunique()}")
    print(f"  Total HST observations: {len(df)}")
    
    return df

def download_hst_images(obs_table, output_dir="hst_images", max_images=1, file_type="drz"):
    """Download HST products using astroquery.
    
    Args:
        obs_table: DataFrame of HST observations
        output_dir: Directory to save downloaded files
        max_images: Maximum number of files to download
        file_type: File type to download (default: 'drz' for drizzle-combined images)
    """
    print(f"\nDownloading HST products (max {max_images} files, type: {file_type.upper()})...")
    print("="*60)
    
    # obs_table is already a DataFrame from check_hst_coverage
    df = obs_table
    
    # Filter for observations with data URLs
    df_with_data = df[df['dataURL'].notna() & (df['dataURL'] != '--')]
    
    if len(df_with_data) == 0:
        print("No data products available.")
        return []
    
    print(f"Found {len(df_with_data)} observations with data products")
    
    # Convert back to astropy Table for astroquery
    from astropy.table import Table
    obs_table_subset = Table.from_pandas(df_with_data.head(max_images * 3))  # Get more obs to find DRZ files
    
    try:
        # Get product list
        print("Getting product list from MAST...")
        products = Observations.get_product_list(obs_table_subset)
        
        print(f"Found {len(products)} products")
        
        # Filter for science products
        science_products = products[products['productType'] == 'SCIENCE']
        
        if len(science_products) == 0:
            print("No science products found. Downloading all products...")
            science_products = products
        
        # Filter by file type (e.g., drz, flt, raw)
        file_type_lower = file_type.lower()
        if 'productFilename' in science_products.colnames:
            mask = [fname.lower().endswith(f'_{file_type_lower}.fits') for fname in science_products['productFilename']]
            filtered_products = science_products[mask]
            
            if len(filtered_products) == 0:
                print(f"No {file_type.upper()} files found. Falling back to first available science product...")
                filtered_products = science_products[:1]
            else:
                print(f"Found {len(filtered_products)} {file_type.upper()} files")
                science_products = filtered_products
        else:
            print("Warning: productFilename column not found, downloading all science products")
        
        # Limit to max_images products
        science_products = science_products[:max_images]
        
        print(f"Downloading {len(science_products)} products...")
        
        # Download products
        manifest = Observations.download_products(
            science_products,
            download_dir=output_dir,
            curl_flag=False
        )
        
        downloaded_files = []
        # Handle manifest - it could be a pandas DataFrame or astropy Table
        try:
            if hasattr(manifest, 'iterrows'):
                for index, row in manifest.iterrows():
                    if row['Status'] == 'COMPLETE':
                        local_path = row['Local Path']
                        downloaded_files.append(local_path)
                        print(f"Downloaded: {os.path.basename(local_path)}")
            elif hasattr(manifest, '__iter__'):
                for row in manifest:
                    if row['Status'] == 'COMPLETE':
                        local_path = row['Local Path']
                        downloaded_files.append(local_path)
                        print(f"Downloaded: {os.path.basename(local_path)}")
            else:
                print(f"Manifest type: {type(manifest)}")
                print(f"Manifest: {manifest}")
        except Exception as e:
            print(f"Error processing manifest: {e}")
            # Fallback: check for files in the download directory matching file type
            import glob
            pattern = f"**/*_{file_type}.fits"
            fits_files = glob.glob(os.path.join(output_dir, pattern), recursive=True)
            # Sort by modification time to get most recent downloads
            fits_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
            # Limit to max_images
            downloaded_files = fits_files[:max_images]
            for f in downloaded_files:
                print(f"Found file: {os.path.basename(f)}")
        
        # Clean up file paths - move files to output directory root
        cleaned_files = []
        for filepath in downloaded_files:
            if os.path.exists(filepath):
                # Extract filename
                filename = os.path.basename(filepath)
                # New path in output directory root
                new_path = os.path.join(output_dir, filename)
                
                # Move if not already in root
                if filepath != new_path:
                    import shutil
                    shutil.move(filepath, new_path)
                    # Clean up empty directories
                    old_dir = os.path.dirname(filepath)
                    while old_dir != output_dir and os.path.exists(old_dir):
                        try:
                            os.rmdir(old_dir)
                            old_dir = os.path.dirname(old_dir)
                        except OSError:
                            break
                
                cleaned_files.append(new_path)
        
        print(f"\nSuccessfully downloaded {len(cleaned_files)} files to {output_dir}")
        return cleaned_files
        
    except Exception as e:
        print(f"Error downloading products: {e}")
        print("\nFalling back to providing MAST portal links:")
        print(f"{'='*100}")
        print(f"{'Instrument':<15} {'Filter':<20} {'Obs ID':<25} {'MAST Portal URL'}")
        print(f"{'='*100}")
        
        for _, row in df_with_data.head(max_images).iterrows():
            instrument = row.get('instrument_name', 'N/A')
            filters = row.get('filters', 'N/A')
            obs_id = row.get('obs_id', 'N/A')
            
            # Construct MAST portal URL
            mast_url = f"https://mast.stsci.edu/portal/Mashup/Mast.ServicePortalv3?inputsearchtype=ObsID&obsid={obs_id}"
            
            print(f"{instrument:<15} {filters:<20} {obs_id:<25} {mast_url}")
        
        print(f"{'='*100}")
        return []


def plot_hst_images(image_files, output_file="hst_mosaic.png", target_ra=None, target_dec=None):
    """Create a two-panel plot from downloaded HST FITS images."""
    print(f"\nCreating plot from {len(image_files)} images...")
    print("="*60)
    
    if len(image_files) == 0:
        print("No images to plot.")
        return
    
    # Plot first image only
    img_file = image_files[0]
    
    try:
        # Read FITS file
        with fits.open(img_file) as hdul:
            # Get the science data extension (usually extension 1)
            data = hdul[1].data
            # Get metadata from primary header (extension 0)
            header = hdul[0].header
            
            # Get WCS information from science extension
            from astropy.wcs import WCS
            wcs = WCS(hdul[1].header)
            
            # Apply zscale scaling using astropy
            from astropy.visualization import ZScaleInterval
            zscale_interval = ZScaleInterval()
            z1, z2 = zscale_interval.get_limits(data)
            
            # Create figure with WCSAxes for proper coordinate display
            from astropy.visualization.wcsaxes import WCSAxes
            fig = plt.figure(figsize=(16, 7))
            
            # Left panel: Full image with WCS
            ax1 = fig.add_subplot(1, 2, 1, projection=wcs)
            im1 = ax1.imshow(data, origin='lower', cmap='viridis', vmin=z1, vmax=z2)
            ax1.set_title('Full Image', fontsize=14)
            ax1.grid(color='white', linestyle=':', linewidth=0.5, alpha=0.5)
            
            # Format RA/DEC labels
            ax1.coords['ra'].set_axislabel('RA', fontsize=12)
            ax1.coords['dec'].set_axislabel('DEC', fontsize=12)
            ax1.coords['ra'].format_unit = 'deg'
            ax1.coords['dec'].format_unit = 'deg'
            
            # Add crosshair at target coordinates if provided
            if target_ra is not None and target_dec is not None:
                # Convert RA/DEC to pixel coordinates
                x_target, y_target = wcs.world_to_pixel_values(target_ra, target_dec)
                
                # Draw crosshair
                ax1.axvline(x_target, color='red', linestyle='--', linewidth=1, alpha=0.7)
                ax1.axhline(y_target, color='red', linestyle='--', linewidth=1, alpha=0.7)
                ax1.plot(x_target, y_target, 'r+', markersize=15, markeredgewidth=2)
            
            plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
            
            # Right panel: 5 arcsec cutout
            ax2 = fig.add_subplot(1, 2, 2)
            if target_ra is not None and target_dec is not None:
                # Convert 5 arcsec to degrees
                cutout_radius_deg = 5.0 / 3600.0  # 5 arcsec in degrees
                
                # Get pixel coordinates of target
                x_target, y_target = wcs.world_to_pixel_values(target_ra, target_dec)
                
                # Convert cutout radius to pixels using pixel scale
                if 'CDELT1' in header and 'CDELT2' in header:
                    pixel_scale = abs(header['CDELT1'])  # degrees per pixel
                    cutout_radius_pix = cutout_radius_deg / pixel_scale
                else:
                    # Default: assume ~0.1 arcsec/pixel for HST
                    cutout_radius_pix = 5.0 / 0.1  # 50 pixels
                
                # Extract cutout
                x_min = int(max(0, x_target - cutout_radius_pix))
                x_max = int(min(data.shape[1], x_target + cutout_radius_pix))
                y_min = int(max(0, y_target - cutout_radius_pix))
                y_max = int(min(data.shape[0], y_target + cutout_radius_pix))
                
                cutout = data[y_min:y_max, x_min:x_max]
                
                # Display cutout with zscale
                zscale_interval_cut = ZScaleInterval()
                z1_cut, z2_cut = zscale_interval_cut.get_limits(cutout)
                im2 = ax2.imshow(cutout, origin='lower', cmap='viridis', vmin=z1_cut, vmax=z2_cut)
                ax2.set_title(f'5 arcsec Cutout', fontsize=14)
                
                # Calculate RA/DEC offsets for axes in arcseconds
                import numpy as np
                ny, nx = cutout.shape
                
                # Create coordinate grids for the cutout
                y_indices, x_indices = np.indices((ny, nx))
                # Convert to image coordinates
                x_img = x_indices + x_min
                y_img = y_indices + y_min
                
                # Convert to world coordinates
                from astropy.coordinates import SkyCoord
                from astropy import units as u
                world_coords = wcs.pixel_to_world(x_img.flatten(), y_img.flatten())
                ra_array = world_coords.ra.deg.reshape(ny, nx)
                dec_array = world_coords.dec.deg.reshape(ny, nx)
                
                # Get RA/DEC at the center of the cutout (this should be the target position)
                center_ra = ra_array[ny//2, nx//2]
                center_dec = dec_array[ny//2, nx//2]
                
                # Calculate offsets from cutout center in arcseconds
                dra_arcsec = (ra_array - center_ra) * 3600.0 * np.cos(np.radians(center_dec))
                ddec_arcsec = (dec_array - center_dec) * 3600.0
                
                # Set tick positions and labels with nice round values
                # Get the range of arcsecond values
                x_min_arcsec = dra_arcsec[ny//2, 0]
                x_max_arcsec = dra_arcsec[ny//2, -1]
                y_min_arcsec = ddec_arcsec[0, nx//2]
                y_max_arcsec = ddec_arcsec[-1, nx//2]
                
                # Ensure min < max for proper tick generation
                x_min_arcsec, x_max_arcsec = min(x_min_arcsec, x_max_arcsec), max(x_min_arcsec, x_max_arcsec)
                y_min_arcsec, y_max_arcsec = min(y_min_arcsec, y_max_arcsec), max(y_min_arcsec, y_max_arcsec)
                
                # Create tick values at 1 arcsec intervals centered on 0
                def nice_ticks(min_val, max_val):
                    """Generate tick values at 1 arcsec intervals centered on 0."""
                    # Generate ticks from floor(min) to ceil(max) at 1 arcsec steps
                    tick_min = np.floor(min_val)
                    tick_max = np.ceil(max_val)
                    ticks = np.arange(tick_min, tick_max + 1, 1.0)
                    return ticks
                
                x_tick_values = nice_ticks(x_min_arcsec, x_max_arcsec)
                y_tick_values = nice_ticks(y_min_arcsec, y_max_arcsec)
                
                # Find pixel positions corresponding to these tick values
                def find_pixel_positions(tick_values, arcsec_array, axis):
                    """Find pixel indices for given tick values."""
                    positions = []
                    for tick in tick_values:
                        # Find the pixel closest to this tick value
                        if axis == 'x':
                            # For x-axis, use middle row
                            diff = np.abs(arcsec_array[ny//2, :] - tick)
                        else:
                            # For y-axis, use middle column
                            diff = np.abs(arcsec_array[:, nx//2] - tick)
                        pos = np.argmin(diff)
                        positions.append(pos)
                    return positions
                
                x_tick_indices = find_pixel_positions(x_tick_values, dra_arcsec, 'x')
                y_tick_indices = find_pixel_positions(y_tick_values, ddec_arcsec, 'y')
                
                # Sort tick indices and corresponding values together to ensure proper ordering
                x_sorted = sorted(zip(x_tick_indices, x_tick_values))
                y_sorted = sorted(zip(y_tick_indices, y_tick_values))
                x_tick_indices, x_tick_values = zip(*x_sorted) if x_sorted else ([], [])
                y_tick_indices, y_tick_values = zip(*y_sorted) if y_sorted else ([], [])
                
                ax2.set_xticks(x_tick_indices)
                ax2.set_yticks(y_tick_indices)
                ax2.set_xticklabels([f'{x:.1f}' for x in x_tick_values])
                ax2.set_yticklabels([f'{y:.1f}' for y in y_tick_values])
                ax2.set_xlabel('ΔRA (arcsec)', fontsize=12)
                ax2.set_ylabel('ΔDEC (arcsec)', fontsize=12)
                
                # Add hollow circle at center
                # Check for FWHM in header, otherwise use default
                fwhm = None
                for key in header.keys():
                    if key and 'FWHM' in key.upper():
                        fwhm = header[key]
                        break
                
                center_x = cutout.shape[1] / 2
                center_y = cutout.shape[0] / 2
                
                if fwhm is not None and fwhm > 0:
                    # Scale circle to 1.7 * FWHM
                    # Convert FWHM (presumably in arcsec) to pixels
                    if 'CDELT1' in header and 'CDELT2' in header:
                        pixel_scale = abs(header['CDELT1']) * 3600.0  # arcsec per pixel
                        circle_radius_pix = (1.7 * fwhm) / pixel_scale
                        # Convert to markersize (roughly proportional to diameter in points)
                        markersize = circle_radius_pix * 2
                    else:
                        markersize = 10  # default
                else:
                    # Default circle size
                    markersize = 10
                
                ax2.plot(center_x, center_y, 'ro', markersize=markersize, markeredgewidth=2, fillstyle='none')
                
                plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
            else:
                ax2.text(0.5, 0.5, 'Target coordinates\nnot provided', ha='center', va='center', 
                         transform=ax2.transAxes, fontsize=14)
                ax2.axis('off')
            
            # Add observation info to the plot
            obs_date = header.get('DATE-OBS', 'Unknown')
            filter_name = header.get('FILTER', 'Unknown')
            instrument = header.get('INSTRUME', 'Unknown')
            
            # Format date
            try:
                from datetime import datetime
                if obs_date != 'Unknown':
                    dt = datetime.strptime(obs_date, '%Y-%m-%d')
                    formatted_date = dt.strftime('%d %B %Y')
                else:
                    formatted_date = 'Unknown'
            except:
                formatted_date = obs_date
            
            # Add info text
            info_text = f'HST {instrument} | {filter_name} | {formatted_date}'
            fig.suptitle(info_text, fontsize=16, y=0.98)
            
            # Adjust layout to prevent label overlap
            plt.tight_layout(rect=[0, 0, 1, 0.96])
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            print(f"Plot saved to: {output_file}")
            
            # Verify file was created
            if os.path.exists(output_file):
                file_size = os.path.getsize(output_file) / 1024  # KB
                print(f"  File size: {file_size:.1f} KB")
            plt.close()
            
    except Exception as e:
        print(f"Error reading {img_file}: {e}")
        import traceback
        traceback.print_exc()

def main():
    parser = argparse.ArgumentParser(description='Check for HST coverage at a given position')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--ra', type=float, help='Right Ascension in degrees')
    group.add_argument('--tns', type=str, help='TNS name to resolve')
    parser.add_argument('--dec', type=float, help='Declination in degrees (required with --ra)')
    parser.add_argument('--radius', type=float, default=0.1, help='Search radius in degrees (default: 0.1)')
    parser.add_argument('--download', action='store_true', help='Download HST images')
    parser.add_argument('--plot', action='store_true', help='Create plot from downloaded images')
    parser.add_argument('--max-images', type=int, default=1, help='Maximum number of images to download (default: 1)')
    parser.add_argument('--file-type', type=str, default='drz', help='File type to download (default: drz, options: drz, flt, raw, etc.)')
    parser.add_argument('--output-dir', type=str, default='hst_images', help='Output directory for images (default: hst_images)')
    
    args = parser.parse_args()
    
    # Get coordinates and set output directory
    if args.tns:
        ra, dec = resolve_tns_name(args.tns)
        if ra is None or dec is None:
            sys.exit(1)
        # Use TNS name as output directory (sanitize)
        output_dir = args.tns.replace(' ', '_').replace('/', '_')
        print(f"Output directory: {output_dir}")
    else:
        if args.dec is None:
            print("Error: --dec is required when using --ra")
            sys.exit(1)
        ra, dec = args.ra, args.dec
        output_dir = args.output_dir
    
    # Create output directory if TNS name was provided
    if args.tns:
        os.makedirs(output_dir, exist_ok=True)
    
    # Check HST coverage
    result = check_hst_coverage(ra, dec, radius=args.radius)
    
    if result is not None:
        print(f"\nHST coverage found! {len(result)} observations available.")
        
        # Download images if requested
        if args.download:
            downloaded_files = download_hst_images(result, output_dir=output_dir, max_images=args.max_images, file_type=args.file_type)
            
            # Plot images if requested
            if args.plot and downloaded_files:
                for img_file in downloaded_files:
                    # Read header to get filter and date
                    try:
                        with fits.open(img_file) as hdul:
                            header = hdul[0].header
                            filter_name = header.get('FILTER', 'Unknown')
                            obs_date = header.get('DATE-OBS', 'Unknown')
                            
                            # Format date for folder name
                            if obs_date != 'Unknown':
                                try:
                                    from datetime import datetime
                                    dt = datetime.strptime(obs_date, '%Y-%m-%d')
                                    date_str = dt.strftime('%Y%m%d')
                                except:
                                    date_str = obs_date.replace('-', '')
                            else:
                                date_str = 'Unknown'
                            
                            # Create subfolder name
                            subfolder_name = f"{filter_name}_{date_str}"
                            subfolder_path = os.path.join(output_dir, subfolder_name)
                            os.makedirs(subfolder_path, exist_ok=True)
                            
                            # Get base filename without extension
                            base_fits_name = os.path.basename(img_file)
                            base_name = os.path.splitext(base_fits_name)[0]
                            
                            # Move FITS file to subfolder if not already there
                            new_fits_path = os.path.join(subfolder_path, base_fits_name)
                            if img_file != new_fits_path and not os.path.exists(new_fits_path):
                                import shutil
                                shutil.move(img_file, new_fits_path)
                                print(f"Moved {base_fits_name} to {subfolder_name}/")
                                img_file = new_fits_path
                            
                            # Create PNG with same name as FITS file
                            plot_filename = f"{base_name}.png"
                            plot_path = os.path.join(subfolder_path, plot_filename)
                            
                            plot_hst_images([img_file], output_file=plot_path, 
                                           target_ra=ra, target_dec=dec)
                    except Exception as e:
                        print(f"Error processing {img_file}: {e}")
                        continue
            elif args.plot:
                print("\nNo images downloaded, skipping plot.")
    else:
        print(f"\nNo HST coverage at this location.")

if __name__ == "__main__":
    main()
