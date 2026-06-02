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
    print(f"  Total HST observations: {len(hst_obs)}")
    
    return df

def download_hst_images(obs_table, output_dir="hst_images", max_images=5):
    """Download HST products using astroquery."""
    print(f"\nDownloading HST products (max {max_images} observations)...")
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
    obs_table_subset = Table.from_pandas(df_with_data.head(max_images))
    
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
        
        print(f"Downloading {len(science_products)} products...")
        
        # Download products
        manifest = Observations.download_products(
            science_products,
            download_dir=output_dir,
            curl_flag=False
        )
        
        downloaded_files = []
        for _, row in manifest.iterrows():
            if row['Status'] == 'COMPLETE':
                downloaded_files.append(row['Local Path'])
                print(f"Downloaded: {os.path.basename(row['Local Path'])}")
        
        print(f"\nSuccessfully downloaded {len(downloaded_files)} files to {output_dir}")
        return downloaded_files
        
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

def plot_hst_images(image_files, output_file="hst_mosaic.png"):
    """Create a plot from downloaded HST JPEG images."""
    print(f"\nCreating plot from {len(image_files)} images...")
    print("="*60)
    
    if len(image_files) == 0:
        print("No images to plot.")
        return
    
    from PIL import Image
    
    # Determine grid size
    n_images = min(len(image_files), 9)  # Max 3x3 grid
    n_cols = min(3, n_images)
    n_rows = (n_images + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(12, 4*n_rows))
    if n_images == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for i, img_file in enumerate(image_files[:n_images]):
        try:
            img = Image.open(img_file)
            axes[i].imshow(img)
            axes[i].set_title(os.path.basename(img_file))
            axes[i].axis('off')
        except Exception as e:
            print(f"Error reading {img_file}: {e}")
            axes[i].text(0.5, 0.5, 'Error', ha='center', va='center')
            axes[i].axis('off')
    
    # Hide unused subplots
    for i in range(n_images, len(axes)):
        axes[i].axis('off')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Check for HST coverage at a given position')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--ra', type=float, help='Right Ascension in degrees')
    group.add_argument('--tns', type=str, help='TNS name to resolve')
    parser.add_argument('--dec', type=float, help='Declination in degrees (required with --ra)')
    parser.add_argument('--radius', type=float, default=0.1, help='Search radius in degrees (default: 0.1)')
    parser.add_argument('--download', action='store_true', help='Download HST images')
    parser.add_argument('--plot', action='store_true', help='Create plot from downloaded images')
    parser.add_argument('--max-images', type=int, default=5, help='Maximum number of images to download (default: 5)')
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
            downloaded_files = download_hst_images(result, output_dir=output_dir, max_images=args.max_images)
            
            # Plot images if requested
            if args.plot and downloaded_files:
                plot_hst_images(downloaded_files, output_file=os.path.join(output_dir, "hst_mosaic.png"))
            elif args.plot:
                print("\nNo images downloaded, skipping plot.")
    else:
        print(f"\nNo HST coverage at this location.")

if __name__ == "__main__":
    main()
