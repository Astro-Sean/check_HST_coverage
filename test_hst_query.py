#!/usr/bin/env python3
"""
Basic HST query test using astroquery examples
"""

from astroquery.mast import Observations

print("Testing HST query with query_object...")
print("="*60)

# Test 1: Query by object name (like the documentation example)
print("\nTest 1: Query M31 region")
obs_table = Observations.query_object("M31", radius="0.1 deg")
print(f"Found {len(obs_table)} observations")
print(f"Columns: {obs_table.colnames}")
print(f"\nFirst 3 rows:")
print(obs_table[:3])

# Test 2: Query HST specifically using query_criteria
print("\n" + "="*60)
print("Test 2: Query HST observations using query_criteria")
obs_hst = Observations.query_criteria(obs_collection="HST", dataproduct_type="image")
print(f"Found {len(obs_hst)} HST observations")
if len(obs_hst) > 0:
    print(f"\nFirst HST observation:")
    print(obs_hst[0])
    print(f"\nAvailable columns: {obs_hst.colnames}")
    
    # Extract RA, DEC, instrument, filter
    df = obs_hst.to_pandas()
    print(f"\nSample of RA, DEC, Instrument, Filter:")
    sample = df[['s_ra', 's_dec', 'instrument_name', 'filters']].head(5)
    print(sample)

print("\n" + "="*60)
print("Test complete!")
