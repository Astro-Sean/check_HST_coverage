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
    
    # Sort by observation date (oldest first)
    if 't_obs_start' in df.columns:
        df = df.sort_values('t_obs_start', ascending=True)
        print(f"\nSorting observations by date (oldest first)...")
    elif 't_min' in df.columns:
        df = df.sort_values('t_min', ascending=True)
        print(f"\nSorting observations by date (oldest first)...")
    
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

def download_hst_images(obs_table, output_dir="hst_images", max_images=1, file_type="flt"):
    """Download HST products using astroquery.
    
    Args:
        obs_table: DataFrame of HST observations
        output_dir: Directory to save downloaded files
        max_images: Maximum number of files to download
        file_type: File type to download (default: 'flt' for flat-fielded images)
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
        
        # Skip files that already exist
        import glob
        already_have = []
        to_download_mask = []
        for i, row in enumerate(science_products):
            fname = row['productFilename']
            existing = glob.glob(os.path.join(output_dir, '**', fname), recursive=True)
            if existing:
                print(f"Already exists, skipping: {fname}")
                already_have.append(existing[0])  # take first match only
                to_download_mask.append(False)
            else:
                to_download_mask.append(True)
        
        science_products = science_products[to_download_mask]
        
        downloaded_files = list(already_have)
        
        if len(science_products) == 0:
            print("All files already downloaded.")
            return downloaded_files
        
        print(f"Downloading {len(science_products)} products...")
        
        # Download products
        manifest = Observations.download_products(
            science_products,
            download_dir=output_dir,
            curl_flag=False
        )
        
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


def plot_hst_images(image_files, output_file="hst_mosaic.png", target_ra=None, target_dec=None, clean_cosmic_rays=False):
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
            
            # Clean cosmic rays if requested
            if clean_cosmic_rays:
                try:
                    from astroscrappy import detect_cosmics
                    import numpy as np
                    
                    # Get parameters from header for informed selection
                    gain = header.get('GAIN', 1.0)  # electrons/ADU
                    readnoise = header.get('READNOISE', 5.0)  # electrons
                    exptime = header.get('EXPTIME', 1.0)  # seconds
                    
                    print(f"  Cleaning cosmic rays (gain={gain}, readnoise={readnoise}, exptime={exptime})...")
                    
                    # Run LA Cosmic with header-informed parameters using astroscrappy
                    cr_cleaned_data, cr_mask = detect_cosmics(
                        data,
                        gain=gain,
                        readnoise=readnoise,
                        sigclip=5.0,
                        sigfrac=0.3,
                        objlim=5.0,
                        satlevel=65535.0
                    )
                    
                    data = cr_cleaned_data
                    print(f"  Cosmic ray cleaning complete")
                except ImportError:
                    print("  Warning: astroscrappy not installed. Skipping cosmic ray cleaning.")
                    print("  Install with: pip install astroscrappy")
                except Exception as e:
                    print(f"  Warning: Cosmic ray cleaning failed: {e}")
                    print("  Continuing with original data")
            
            # Apply zscale scaling using astropy
            from astropy.visualization import ZScaleInterval
            zscale_interval = ZScaleInterval()
            z1, z2 = zscale_interval.get_limits(data)
            
            # Create figure with single WCSAxes panel
            from astropy.visualization.wcsaxes import WCSAxes
            fig = plt.figure(figsize=(10, 9))
            
            # Full image with WCS
            ax1 = fig.add_subplot(1, 1, 1, projection=wcs)
            im1 = ax1.imshow(data, origin='lower', cmap='viridis', vmin=z1, vmax=z2)
            ax1.grid(color='white', linestyle=':', linewidth=0.5, alpha=0.5)
            
            # Format RA/DEC labels
            ax1.coords['ra'].set_axislabel('RA (deg)', fontsize=12)
            ax1.coords['dec'].set_axislabel('Dec (deg)', fontsize=12)
            ax1.coords['ra'].set_ticklabel(size=10)
            ax1.coords['dec'].set_ticklabel(size=10)
            
            # Check if target is within bounds for cutout
            target_in_bounds = False
            if target_ra is not None and target_dec is not None:
                # Convert RA/DEC to pixel coordinates
                x_target, y_target = wcs.world_to_pixel_values(target_ra, target_dec)
                
                # Check if target is within the detector bounds
                ny, nx = data.shape
                if 0 <= x_target < nx and 0 <= y_target < ny:
                    target_in_bounds = True
                else:
                    print(f"Warning: Target coordinates ({target_ra:.6f}, {target_dec:.6f}) fall outside the image bounds.")
                    print(f"  Pixel coordinates: ({x_target:.1f}, {y_target:.1f}) outside [0:{nx}, 0:{ny}]")
            
            # Add colorbar with units from header
            bunit = header.get('BUNIT', '')
            if bunit:
                label = f'Intensity ({bunit})'
            else:
                label = 'Intensity'
            plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, label=label)
            
            # Inset axes in top-right corner: 5 arcsec wide (2.5 arcsec radius) cutout
            from mpl_toolkits.axes_grid1.inset_locator import inset_axes
            ax2 = inset_axes(ax1, width='35%', height='35%', loc='upper right')
            if target_ra is not None and target_dec is not None and target_in_bounds:
                # Convert 2.5 arcsec radius to degrees (gives 5 arcsec wide cutout)
                cutout_radius_deg = 2.5 / 3600.0
                
                # Get pixel coordinates of target
                x_target, y_target = wcs.world_to_pixel_values(target_ra, target_dec)
                
                # Convert cutout radius to pixels using pixel scale derived from WCS
                import astropy.units as u
                pixel_scale_deg = wcs.proj_plane_pixel_scales()[0].to(u.deg).value
                cutout_radius_pix = cutout_radius_deg / pixel_scale_deg
                
                # Extract cutout
                x_min = int(max(0, x_target - cutout_radius_pix))
                x_max = int(min(data.shape[1], x_target + cutout_radius_pix))
                y_min = int(max(0, y_target - cutout_radius_pix))
                y_max = int(min(data.shape[0], y_target + cutout_radius_pix))
                
                cutout = data[y_min:y_max, x_min:x_max]
                
                # Add rectangle on full image to mark cutout region
                from matplotlib.patches import Rectangle, ConnectionPatch
                rect = Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                               linewidth=2, edgecolor='red', facecolor='none')
                ax1.add_patch(rect)
                
                # Display cutout with zscale
                zscale_interval_cut = ZScaleInterval()
                z1_cut, z2_cut = zscale_interval_cut.get_limits(cutout)
                im2 = ax2.imshow(cutout, origin='lower', cmap='viridis', vmin=z1_cut, vmax=z2_cut)
                
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
                
                # Map arcsec values to pixel positions using linear interpolation
                # dra_arcsec at middle row: pixel 0 -> left edge, pixel nx-1 -> right edge
                dra_row = dra_arcsec[ny//2, :]        # shape (nx,): arcsec at each x pixel
                ddec_col = ddec_arcsec[:, nx//2]      # shape (ny,): arcsec at each y pixel
                
                # Pixel arrays (0-indexed)
                x_pixels = np.arange(nx, dtype=float)
                y_pixels = np.arange(ny, dtype=float)
                
                # Arcsec range (in pixel order, may be reversed for RA)
                x_arcsec_range = (dra_row[0], dra_row[-1])
                y_arcsec_range = (ddec_col[0], ddec_col[-1])
                
                # Create tick values at 1 arcsec intervals within the actual range
                def nice_ticks(arcsec_start, arcsec_end, step=1.0):
                    """Generate tick values at given step within the range."""
                    lo, hi = min(arcsec_start, arcsec_end), max(arcsec_start, arcsec_end)
                    return np.arange(np.ceil(lo / step) * step, np.floor(hi / step) * step + step * 0.5, step)
                
                x_tick_arcsec = nice_ticks(*x_arcsec_range)
                y_tick_arcsec = nice_ticks(*y_arcsec_range)
                
                # Convert arcsec tick values to pixel positions via linear interpolation
                # np.interp requires xp to be increasing, so sort by arcsec value
                x_tick_pix = np.interp(x_tick_arcsec,
                                       sorted([dra_row[0], dra_row[-1]]),
                                       [0, nx - 1] if dra_row[0] < dra_row[-1] else [nx - 1, 0])
                y_tick_pix = np.interp(y_tick_arcsec,
                                       sorted([ddec_col[0], ddec_col[-1]]),
                                       [0, ny - 1] if ddec_col[0] < ddec_col[-1] else [ny - 1, 0])
                
                # Set ticks and labels
                ax2.set_xticks(x_tick_pix)
                ax2.set_yticks(y_tick_pix)
                ax2.set_xticklabels([f'{x:.1f}' for x in x_tick_arcsec])
                ax2.set_yticklabels([f'{y:.1f}' for y in y_tick_arcsec])
                ax2.set_xlabel('Δα (")', fontsize=8)
                ax2.set_ylabel('Δδ (")', fontsize=8)
                ax2.tick_params(labelsize=7)
                
                # Move x-axis to top and y-axis to right to avoid connector overlap
                ax2.xaxis.tick_top()
                ax2.xaxis.set_label_position('top')
                ax2.yaxis.tick_right()
                ax2.yaxis.set_label_position('right')
                
                # Set axis limits to match the cutout
                ax2.set_xlim(-0.5, nx - 0.5)
                ax2.set_ylim(-0.5, ny - 0.5)
                
                # Draw connector lines from cutout rectangle corners to inset corners
                # Connect to bottom-left and top-right corners to avoid axis labels
                from matplotlib.patches import ConnectionPatch
                # bottom-left of rect -> bottom-left of inset
                con1 = ConnectionPatch(
                    xyA=(x_min, y_min), coordsA=ax1.transData,
                    xyB=(0, 0),         coordsB=ax2.transData,
                    color='red', lw=1.5, linestyle='--', zorder=100)
                # top-right of rect -> top-right of inset
                con2 = ConnectionPatch(
                    xyA=(x_max, y_max), coordsA=ax1.transData,
                    xyB=(nx - 1, ny - 1), coordsB=ax2.transData,
                    color='red', lw=1.5, linestyle='--', zorder=100)
                fig.add_artist(con1)
                fig.add_artist(con2)
                
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
                    pixel_scale_arcsec = pixel_scale_deg * 3600.0  # arcsec per pixel
                    circle_radius_pix = (1.7 * fwhm) / pixel_scale_arcsec
                    markersize = circle_radius_pix * 2
                else:
                    # Default circle size
                    markersize = 10
                
                ax2.plot(center_x, center_y, 'ro', markersize=markersize, markeredgewidth=2, fillstyle='none')
            
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
            
            # Add observation info as text inside the plot (bottom-left)
            info_text = f'HST {instrument} | {filter_name} | {formatted_date}'
            ax1.text(0.01, 0.01, info_text, transform=ax1.transAxes, fontsize=10,
                     color='white', va='bottom', ha='left',
                     bbox=dict(facecolor='black', alpha=0.4, edgecolor='none', pad=3))
            
            # Add compass rose in bottom-right showing North and East (red)
            # Use axes fraction coordinates for reliable positioning
            compass_x, compass_y = 0.85, 0.10  # bottom-right in axes fraction (avoiding inset and edges)
            arrow_len = 0.06  # length in axes fraction (smaller to fit within image)
            
            # Draw North arrow (pointing up in axes coordinates)
            ax1.annotate('N', xy=(compass_x, compass_y + arrow_len), xytext=(compass_x, compass_y),
                        xycoords='axes fraction', textcoords='axes fraction',
                        arrowprops=dict(arrowstyle='->', color='red', lw=0.5),
                        ha='center', va='bottom', color='red', fontsize=16, zorder=300)
            
            # Draw East arrow (pointing right in axes coordinates)
            ax1.annotate('E', xy=(compass_x + arrow_len * 1.5, compass_y), xytext=(compass_x, compass_y),
                        xycoords='axes fraction', textcoords='axes fraction',
                        arrowprops=dict(arrowstyle='->', color='red', lw=0.5),
                        ha='left', va='center', color='red', fontsize=16, zorder=300)
            
            # Adjust layout
            plt.tight_layout()
            
            # Only save plot if target is within bounds
            if target_in_bounds or target_ra is None or target_dec is None:
                plt.savefig(output_file, dpi=150, bbox_inches='tight')
                print(f"Plot saved to: {output_file}")
                
                # Verify file was created
                if os.path.exists(output_file):
                    file_size = os.path.getsize(output_file) / 1024  # KB
                    print(f"  File size: {file_size:.1f} KB")
            else:
                print(f"Skipping plot save: Target outside image bounds")
                plt.close()
                return
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
    parser.add_argument('--file-type', type=str, default='flt', help='File type to download (default: flt, options: flt, drz, crj, etc.)')
    parser.add_argument('--output-dir', type=str, default='hst_images', help='Output directory for images (default: hst_images)')
    parser.add_argument('--clean-cosmic-rays', action='store_true', help='Clean cosmic rays using ccdproc (requires ccdproc installation)')
    
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
                                           target_ra=ra, target_dec=dec, 
                                           clean_cosmic_rays=args.clean_cosmic_rays)
                    except Exception as e:
                        print(f"Error processing {img_file}: {e}")
                        continue
            elif args.plot:
                print("\nNo images downloaded, skipping plot.")
    else:
        print(f"\nNo HST coverage at this location.")

if __name__ == "__main__":
    main()
