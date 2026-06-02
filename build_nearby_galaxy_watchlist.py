"""
Build a Lasair watchlist CSV of nearby galaxies from three catalogs:
  - HECATE v2  (VizieR J/MNRAS/548/G522) : ~179k galaxies, d <= 200 Mpc, D25 radii
  - GLADE+     (VizieR VII/291/gladep)    : ~450k galaxies, dL <= 200 Mpc, fallback radius
  - CLU        (VizieR J/ApJ/880/7/table2): ~473 northern-sky galaxies, d <= 200 Mpc

Output: nearby_galaxies_watchlist.csv  (name, ra, dec, radius_arcsec)
Radius = R1 (semi-major axis in arcmin) where available from HECATE, else 60 arcsec.
"""

import warnings
warnings.filterwarnings("ignore")

import sys
import time
import os
import numpy as np
import pandas as pd
import requests
import io
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

FALLBACK_RADIUS_ARCSEC = 60.0   # 1 arcmin fallback
D_MAX_MPC              = 200.0
VLT_DEC_MAX            = 35.0    # VLT visibility limit (degrees)
TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"

# Distance-based radius: assume typical galaxy size ~10 kpc
# angular_size (arcsec) = (size_kpc / dist_mpc) * 206.265
# For 10 kpc galaxy: radius = 5 kpc = (5 / dist_mpc) * 206.265 arcsec
def dist_to_radius_arcsec(dist_mpc, fallback=FALLBACK_RADIUS_ARCSEC):
    """Estimate radius from distance assuming ~10 kpc galaxy size."""
    try:
        d = float(dist_mpc)
        if d > 0 and np.isfinite(d):
            r = (5.0 / d) * 206.265  # 5 kpc semi-major axis in arcsec
            return round(min(max(r, 5.0), 600.0), 1)  # cap 5-600 arcsec
    except Exception:
        pass
    return fallback

# ── helpers ───────────────────────────────────────────────────────────────────

def tap_query(adql, timeout=300):
    """Execute a synchronous VizieR TAP query and return a DataFrame."""
    params = {
        "REQUEST": "doQuery",
        "LANG":    "ADQL",
        "FORMAT":  "csv",
        "QUERY":   adql,
    }
    resp = requests.get(TAP_URL, params=params, timeout=timeout)
    resp.raise_for_status()
    # Check for VOTABLE error embedded in 200-response
    if "QUERY_STATUS" in resp.text and "ERROR" in resp.text:
        raise ValueError(f"TAP error: {resp.text[:300]}")
    t = pd.read_csv(io.StringIO(resp.text))
    return t


def r1_to_radius_arcsec(r1_arcmin, fallback=FALLBACK_RADIUS_ARCSEC):
    """
    R1 is the semi-major axis in arcmin (HECATE convention).
    Convert to arcsec. Uses fallback if missing/zero.
    """
    try:
        r = float(r1_arcmin) * 60.0  # arcmin -> arcsec
        if np.isfinite(r) and r >= 3.0:
            return round(min(r, 600.0), 1)  # cap at 10 arcmin
    except Exception:
        pass
    return fallback


def clean_df(df):
    """Drop rows missing RA/Dec, enforce valid coordinate ranges."""
    df = df.dropna(subset=["ra", "dec"]).copy()
    df["ra"]  = pd.to_numeric(df["ra"],  errors="coerce")
    df["dec"] = pd.to_numeric(df["dec"], errors="coerce")
    df = df.dropna(subset=["ra", "dec"])
    df = df[(df["ra"] >= 0) & (df["ra"] < 360)]
    df = df[(df["dec"] >= -90) & (df["dec"] <= 90)]
    df.reset_index(drop=True, inplace=True)
    return df


# ── 1. HECATE v2 ──────────────────────────────────────────────────────────────

def fetch_hecate(d_max=D_MAX_MPC, m_min=-17.0):
    print(f"\n[HECATE v2] Querying VizieR TAP for galaxies within {d_max} Mpc …")
    adql = (
        'SELECT OBJNAME, RAJ2000, DEJ2000, Dist, R1 '
        'FROM "J/MNRAS/548/G522/hecatev2" '
        f'WHERE Dist <= {d_max} AND Dist > 0'
    )
    try:
        t = tap_query(adql, timeout=300)
        print(f"  [HECATE] Raw rows: {len(t)}")
        print(f"  Columns: {list(t.columns)}")
    except Exception as e:
        print(f"  [HECATE] Query failed: {e}")
        return pd.DataFrame()

    t["name"]     = t["OBJNAME"].astype(str).str.strip().str.strip('"')
    t["ra"]       = pd.to_numeric(t["RAJ2000"], errors="coerce")
    t["dec"]      = pd.to_numeric(t["DEJ2000"], errors="coerce")
    t["dist_mpc"] = pd.to_numeric(t["Dist"],    errors="coerce")
    t["radius_arcsec"] = t["R1"].apply(r1_to_radius_arcsec)
    t["catalog"]  = "HECATE"

    t = clean_df(t[["name", "ra", "dec", "dist_mpc", "radius_arcsec", "catalog"]])
    print(f"  [HECATE] After cleaning: {len(t)} galaxies")
    r_stats = t["radius_arcsec"]
    print(f"  Radius range: {r_stats.min():.1f}–{r_stats.max():.1f} arcsec  "
          f"(median {r_stats.median():.1f})")
    return t


# ── 2. GLADE+ ─────────────────────────────────────────────────────────────────

def fetch_glade(d_max=D_MAX_MPC, m_min=-17.0):
    print(f"\n[GLADE+] Querying VizieR TAP for galaxies within {d_max} Mpc …")
    adql = (
        'SELECT GWGC, "HyperLEDA", "2MASS", RAJ2000, DEJ2000, dL '
        'FROM "VII/291/gladep" '
        f'WHERE dL <= {d_max} AND Type = \'G\''
    )
    try:
        t = tap_query(adql, timeout=600)
        print(f"  [GLADE+] Raw rows: {len(t)}")
        print(f"  Columns: {list(t.columns)}")
    except Exception as e:
        print(f"  [GLADE+] Query failed: {e}")
        return pd.DataFrame()

    def best_name(row):
        for col in ["GWGC", "HyperLEDA", "2MASS"]:
            if col in row.index:
                val = str(row[col]).strip().strip('"').strip("'")
                if val and val not in ("-", "--", "nan", ""):
                    return val
        return f"GLADE_{row.name}"

    t["name"]     = t.apply(best_name, axis=1)
    t["ra"]       = pd.to_numeric(t["RAJ2000"], errors="coerce")
    t["dec"]      = pd.to_numeric(t["DEJ2000"], errors="coerce")
    t["dist_mpc"] = pd.to_numeric(t["dL"],      errors="coerce")
    t["radius_arcsec"] = t["dist_mpc"].apply(dist_to_radius_arcsec)
    t["catalog"]  = "GLADE+"

    t = clean_df(t[["name", "ra", "dec", "dist_mpc", "radius_arcsec", "catalog"]])
    print(f"  [GLADE+] After cleaning: {len(t)} galaxies")
    return t


# ── 3. CLU ────────────────────────────────────────────────────────────────────

def fetch_clu(d_max=D_MAX_MPC, m_min=-17.0):
    print(f"\n[CLU] Querying VizieR TAP for galaxies within {d_max} Mpc …")
    adql = (
        'SELECT CLU, RAJ2000, DEJ2000, Dist '
        'FROM "J/ApJ/880/7/table2" '
        f'WHERE Dist <= {d_max}'
    )
    try:
        t = tap_query(adql, timeout=120)
        print(f"  [CLU] Raw rows: {len(t)}")
    except Exception as e:
        print(f"  [CLU] Query failed: {e}")
        return pd.DataFrame()

    t["name"]     = t["CLU"].astype(str).str.strip().str.strip('"')
    t["ra"]       = pd.to_numeric(t["RAJ2000"], errors="coerce")
    t["dec"]      = pd.to_numeric(t["DEJ2000"], errors="coerce")
    t["dist_mpc"] = pd.to_numeric(t["Dist"],    errors="coerce")
    t["radius_arcsec"] = t["dist_mpc"].apply(dist_to_radius_arcsec)
    t["catalog"]  = "CLU"

    t = clean_df(t[["name", "ra", "dec", "dist_mpc", "radius_arcsec", "catalog"]])
    print(f"  [CLU] After cleaning: {len(t)} galaxies")
    return t


# ── 4. Visualization ───────────────────────────────────────────────────────────

def plot_skymap(df, distance_mpc, output_path):
    """Create skymap with points colored by distance (VLT-visible only)."""
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection='aitoff')

    # Convert RA to radians (-180 to 180 for aitoff)
    ra_rad = np.radians(df['ra'] - 180)
    dec_rad = np.radians(df['dec'])

    # Color by distance
    norm = Normalize(vmin=0, vmax=distance_mpc)
    cmap = plt.cm.viridis

    scatter = ax.scatter(ra_rad, dec_rad, c=df['dist_mpc'],
                        s=1, alpha=0.6, cmap=cmap, norm=norm)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, orientation='horizontal',
                       pad=0.05, aspect=40)
    cbar.set_label('Distance (Mpc)')

    ax.set_title(f'VLT-Visible Nearby Galaxies (≤ {distance_mpc} Mpc, dec ≤ 35°) - {len(df):,} galaxies',
                 fontsize=14)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Skymap saved to '{output_path}'")


# ── 5. Merge + deduplicate ────────────────────────────────────────────────────

def deduplicate(df, sep_arcsec=1.0):
    """
    Remove spatial duplicates within sep_arcsec using a coarse coordinate grid.
    Priority order: HECATE (best radii) > CLU > GLADE+
    """
    priority = {"HECATE": 0, "CLU": 1, "GLADE+": 2}
    df = df.copy()
    df["_prio"] = df["catalog"].map(priority).fillna(9).astype(int)
    df = df.sort_values("_prio")

    grid = sep_arcsec / 3600.0
    df["_ra_r"]  = (df["ra"]  / grid).round().astype(int)
    df["_dec_r"] = (df["dec"] / grid).round().astype(int)
    df = df.drop_duplicates(subset=["_ra_r", "_dec_r"], keep="first")
    df = df.drop(columns=["_prio", "_ra_r", "_dec_r"])
    return df.reset_index(drop=True)


# ── 5. Main ───────────────────────────────────────────────────────────────────

def build_watchlist(distance_mpc):
    """Build watchlist for a specific distance cut."""
    print(f"\n{'='*60}")
    print(f"Building watchlist for D <= {distance_mpc} Mpc")
    print(f"{'='*60}")
    t0 = time.time()
    frames = []

    hecate = fetch_hecate(d_max=distance_mpc)
    if not hecate.empty:
        frames.append(hecate)

    glade = fetch_glade(d_max=distance_mpc)
    if not glade.empty:
        frames.append(glade)

    clu = fetch_clu(d_max=distance_mpc)
    if not clu.empty:
        frames.append(clu)

    if not frames:
        print(f"\n[ERROR] All catalogs returned empty for {distance_mpc} Mpc.")
        return None

    combined = pd.concat(frames, ignore_index=True)
    print(f"\n[Merge] Total before dedup: {len(combined):,}")

    merged = deduplicate(combined, sep_arcsec=1.0)
    print(f"[Merge] After dedup (1 arcsec grid): {len(merged):,}")

    # Enforce radius bounds
    merged["radius_arcsec"] = merged["radius_arcsec"].clip(5.0, 600.0)

    # Filter for VLT visibility (dec <= 35°)
    merged_vlt = merged[merged["dec"] <= VLT_DEC_MAX].copy()
    print(f"[VLT filter] {len(merged_vlt):,} galaxies visible from VLT (dec <= {VLT_DEC_MAX}°)")

    # Per-catalog summary
    print("\n── Catalog breakdown ──────────────────────────")
    for cat, grp in merged.groupby("catalog"):
        print(f"  {cat:10s}: {len(grp):7,d} galaxies")
    print(f"  {'TOTAL':10s}: {len(merged):7,d} galaxies")

    # Build Lasair watchlist CSV with radius - VLT filtered
    output = merged_vlt[["ra", "dec", "name", "radius_arcsec"]].copy()
    output.columns = ["RA", "DEC", "ID", "Radius"]
    output["RA"]     = output["RA"].round(6)
    output["DEC"]    = output["DEC"].round(6)
    output["Radius"] = output["Radius"].round(1)

    elapsed = time.time() - t0
    print(f"\n[Output] Total entries (VLT-visible): {len(output):,}")
    print(f"  RA range : {output['RA'].min():.2f} – {output['RA'].max():.2f} deg")
    print(f"  Dec range: {output['DEC'].min():.2f} – {output['DEC'].max():.2f} deg")
    print(f"  Radius   : {output['Radius'].min():.1f} – {output['Radius'].max():.1f} arcsec"
          f"  (median {output['Radius'].median():.1f})")
    print(f"  Elapsed  : {elapsed:.0f}s")

    return output, merged_vlt


def main():
    distance_cuts = [10, 50, 100, 200]  # Mpc
    all_outputs = {}
    all_merged = {}

    for d_mpc in distance_cuts:
        output, merged = build_watchlist(d_mpc)
        if output is not None:
            all_outputs[d_mpc] = output
            all_merged[d_mpc] = merged

    # Write all files to distance-based folders
    print(f"\n{'='*60}")
    print("Writing output files")
    print(f"{'='*60}")

    for d_mpc, output in all_outputs.items():
        # Create folder for this distance
        folder_name = f"{d_mpc}Mpc"
        os.makedirs(folder_name, exist_ok=True)

        # Full watchlist for this distance
        out_path = os.path.join(folder_name, "nearby_galaxies_watchlist.csv")
        output.to_csv(out_path, index=False)
        print(f"  Written {len(output):,} entries to '{out_path}'")

        # Compact 100k version if needed (append _2 since it exceeds limit)
        if len(output) > 100000:
            output_100k = output.head(100000)
            out_path_100k = os.path.join(folder_name, "nearby_galaxies_watchlist_2.csv")
            output_100k.to_csv(out_path_100k, index=False)
            print(f"    + {len(output_100k):,} entries to '{out_path_100k}'")

        # Create skymap with distance coloring
        merged_df = all_merged[d_mpc]
        plot_path = os.path.join(folder_name, "skymap.png")
        plot_skymap(merged_df, d_mpc, plot_path)

    print(f"\nUpload the compact files to your Lasair watchlist at:")
    print("  https://lasair-ztf.lsst.ac.uk/watchlists/")


if __name__ == "__main__":
    main()
