#!/usr/bin/env python3
"""
Build HST observation watchlist for LSST cross-matching
Queries MAST for HST observations and extracts RA, DEC, filter, and instrument
"""

import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from astroquery.mast import Observations

def fetch_hst_observations(limit=10000):
    """Fetch HST observations from MAST."""
    print(f"\nFetching HST observations (limit={limit})...")
    
    # Query HST observations for image data

    obs_table = Observations.query_criteria( # Exact match on data product typ
                                          
                                            dataproduct_type = 'image',
                                            dataRights = 'PUBLIC',
                                            calib_level = ['1','2'],
                                            t_min = ">60000")  # Range match on minimum wavelength
    print(obs_table.to_pandas()[:1].T)
    print(obs_table.columns)
    
    # Limit results for speed
    if len(obs_table) > limit:
        obs_table = obs_table[:limit]
        print(f"Limited to {limit} observations for speed")
    
    print(f"Total HST observations retrieved: {len(obs_table):,}")
    return obs_table

def process_hst_observations(obs_table):
    """Process HST observations to extract RA, DEC, filter, and instrument."""
    print("\nProcessing HST observations...")
    
    # Convert to pandas DataFrame
    df = obs_table.to_pandas()
    
    # Check for required columns
    required_cols = ['s_ra', 's_dec', 'instrument_name', 'filters']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        print(f"Available columns: {df.columns.tolist()}")
        return None
    
    # Extract basic info
    output = df[required_cols].copy()
    output.columns = ["RA", "DEC", "Instrument", "Filter"]
    
    # Round coordinates
    output["RA"] = output["RA"].round(6)
    output["DEC"] = output["DEC"].round(6)
    
    print(f"Processed {len(output):,} observations")
    if len(output) > 0:
        print(f"  RA range : {output['RA'].min():.2f} – {output['RA'].max():.2f} deg")
        print(f"  Dec range: {output['DEC'].min():.2f} – {output['DEC'].max():.2f} deg")
        print(f"  Instruments: {output['Instrument'].nunique()} unique")
        print(f"  Filters: {output['Filter'].nunique()} unique")
    else:
        print("  No observations to process")

    return output

def main():
    print("="*60)
    print("HST Observation Watchlist Builder")
    print("="*60)

    # Fetch observations
    obs_table = fetch_hst_observations(limit=10000)

    # Process observations
    output = process_hst_observations(obs_table)

    if output is None:
        print("Error processing observations.")
        return

    # Save to CSV
    out_path = "hst_observations_watchlist.csv"
    output.to_csv(out_path, index=False)
    print(f"\nWritten {len(output):,} entries to '{out_path}'")
    print(f"\nFormat: RA, DEC, Instrument, Filter")

if __name__ == "__main__":
    main()
