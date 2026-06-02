# Nearby Galaxy Watchlist for Lasair

Watchlist CSVs of nearby galaxies for transient searches on [Lasair](https://lasair-lsst.lsst.ac.uk/) (ZTF/LSST alert broker).

## Output Files

Watchlists are organized in distance-based folders and filtered for VLT visibility (dec ≤ 35°):

```
10Mpc/
  └── nearby_galaxies_watchlist.csv (3,956 entries)

50Mpc/
  └── nearby_galaxies_watchlist.csv (35,021 entries)

100Mpc/
  ├── nearby_galaxies_watchlist.csv (112,448 entries)
  └── nearby_galaxies_watchlist_2.csv (100,000 entries)

200Mpc/
  ├── nearby_galaxies_watchlist.csv (437,711 entries)
  └── nearby_galaxies_watchlist_2.csv (100,000 entries)
```

- **nearby_galaxies_watchlist.csv**: VLT-visible galaxies (dec ≤ 35°) for this distance
- **nearby_galaxies_watchlist_2.csv**: First 100k entries (created only when full list exceeds 100k)

## Catalogs Used

| Catalog | Source | Galaxies | Distance Cut | Radius Method |
|---------|--------|----------|-------------|---------------|
| **HECATE v2** | VizieR `J/MNRAS/548/G522` | 179,318 | ≤ 200 Mpc | R1 semi-major axis (arcmin → arcsec) |
| **GLADE+** | VizieR `VII/291/gladep` | 367,093 | dL ≤ 200 Mpc | Distance-based: (5 kpc / dL) × 206.265 arcsec |
| **CLU** | VizieR `J/ApJ/880/7/table2` | 436 | ≤ 200 Mpc | Distance-based: (5 kpc / Dist) × 206.265 arcsec |

Deduplication: 1 arcsec spatial grid, priority HECATE > CLU > GLADE+.  
Radius bounds: 5 – 600 arcsec.  
Distance-based radius assumes ~10 kpc galaxy diameter (5 kpc semi-major axis).

## CSV Format

```
RA,DEC,ID,Radius
10.684684,41.268978,NGC0224,600.0
```

- `RA`, `DEC` — J2000 decimal degrees
- `ID` — galaxy identifier
- `Radius` — cone radius in **arcseconds** (5–600 arcsec, based on galaxy size)

## Uploading to Lasair

### Watchlist Names and Descriptions

Use these suggested names and descriptions when creating watchlists in Lasair:

| File | Name | Description |
|------|------|-------------|
| `10Mpc/nearby_galaxies_watchlist.csv` | Nearby Galaxies 10 Mpc | VLT-visible galaxies within 10 Mpc from HECATE, GLADE+, and CLU catalogs (3,956 sources) |
| `50Mpc/nearby_galaxies_watchlist.csv` | Nearby Galaxies 50 Mpc | VLT-visible galaxies within 50 Mpc from HECATE, GLADE+, and CLU catalogs (35,021 sources) |
| `100Mpc/nearby_galaxies_watchlist.csv` | Nearby Galaxies 100 Mpc | VLT-visible galaxies within 100 Mpc from HECATE, GLADE+, and CLU catalogs (112,448 sources) |
| `100Mpc/nearby_galaxies_watchlist_2.csv` | Nearby Galaxies 100 Mpc (compact) | VLT-visible galaxies within 100 Mpc - first 100k sources (100,000 sources) |
| `200Mpc/nearby_galaxies_watchlist.csv` | Nearby Galaxies 200 Mpc | VLT-visible galaxies within 200 Mpc from HECATE, GLADE+, and CLU catalogs (437,711 sources) |
| `200Mpc/nearby_galaxies_watchlist_2.csv` | Nearby Galaxies 200 Mpc (compact) | VLT-visible galaxies within 200 Mpc - first 100k sources (100,000 sources) |

### Web upload (≤ 100,000 entries)
1. Log in at https://lasair-lsst.lsst.ac.uk/
2. Go to **Watchlists → Create New Watchlist**
3. Upload one of the compact files:
   - `100Mpc/nearby_galaxies_watchlist_2.csv` (100k entries)
   - `200Mpc/nearby_galaxies_watchlist_2.csv` (100k entries)
4. Set default radius = 60 (arcsec)

### Smaller distance cuts (direct web upload)
The 10 Mpc and 50 Mpc full lists are small enough for direct upload:
- `10Mpc/nearby_galaxies_watchlist.csv` (3,956 entries)
- `50Mpc/nearby_galaxies_watchlist.csv` (35,021 entries)

### Large watchlists (> 100k entries — email required)
For the full 100 Mpc and 200 Mpc lists, use the official process:
1. Create an **empty** watchlist via the web interface
2. Note the watchlist ID from the URL
3. Email [lasair-help@mlist.is.ed.ac.uk](mailto:lasair-help@mlist.is.ed.ac.uk) with:
   - Your watchlist ID
   - The CSV file (or a download link)

## Regenerating

```bash
python3 build_nearby_galaxy_watchlist.py
```

Requires: `astroquery`, `astropy`, `pandas`, `numpy`, `requests`

## Visualization

Each distance folder includes a `skymap.png` showing:
- VLT-visible galaxies only (dec ≤ 35°) in Aitoff projection
- Galaxies plotted with color indicating distance (0 to distance limit)
- Colorbar showing distance in Mpc
- Dense regions correspond to the galactic plane and galaxy clusters
- Warmer colors (yellow/green) indicate more distant galaxies, cooler colors (purple/blue) indicate closer galaxies

## Notes

- **LSST Lasair** URL: https://lasair-lsst.lsst.ac.uk/
- **ZTF Lasair** URL: https://lasair-ztf.lsst.ac.uk/
- Median cone radius: ~19 arcsec (HECATE entries), 60 arcsec (GLADE+/CLU entries)
- Large nearby galaxies (M31, LMC, M101, etc.) are capped at 600 arcsec (10 arcmin)
- GLADE+ completeness at 200 Mpc is ~50% in B-band; HECATE is more complete for well-studied galaxies
