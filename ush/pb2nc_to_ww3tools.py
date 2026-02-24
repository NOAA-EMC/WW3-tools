"""
Rewrite daily MET pb2nc NetCDF files to daily ww3tools readable NetCDF files.

Input pattern:
  jason3_b031_xx124_yyyymmdd.nc

Output pattern:
  Altimeter_JASON3_yyyymmdd.nc

Outputs (1D over time):
  time, lat, lon, hs, hs_cal, wsp, wsp_cal

Mapping:
  time     : seconds since epoch
  lat, lon : observation location
  hs      <- KBSW
  wsp     <- WSPA
  hs_cal  <- NaN
  wsp_cal <- NaN
"""

import os
import re
import numpy as np
import pandas as pd
from datetime import datetime, timezone
from netCDF4 import Dataset, chartostring

# =========================
# USER-DEFINED SETTINGS
# =========================
WORKDIR = "/scratch3/NCEPDEV/marine/Ming.Chen/ursa/wind_eval"            # working directory contains WW3-tools

SAT_NAME = "JASON3"
HS_MNEM  = "KBSW"
WSP_MNEM = "WSPA"

# ========================
PB2NC_DIR = os.path.join(WORKDIR, "processsatdata", "pb2nc_out")
OUT_DIR   = os.path.join(WORKDIR, "processsatdata", "pb2nc_altimeter")

os.makedirs(OUT_DIR, exist_ok=True)
# Find input files
infiles = sorted([
    os.path.join(PB2NC_DIR, f) for f in os.listdir(PB2NC_DIR)
    if re.match(r"jason3_b031_xx124_\d{8}\.nc$", f)
])

if not infiles:
    print("No pb2nc files found.")
    exit()

# =======================
for infile in infiles:
    match = re.search(r"(\d{8})\.nc$", os.path.basename(infile))
    if not match:
        continue
    ymd = match.group(1)
    outfile = os.path.join(OUT_DIR, f"Altimeter_{SAT_NAME}_{ymd}.nc")

    print(f"\nProcessing {os.path.basename(infile)} → {os.path.basename(outfile)}")

    with Dataset(infile, "r") as nc:
        # Read string tables
        obs_var_list = [str(s).strip().strip("\x00") for s in chartostring(nc.variables["obs_var"][:])]
        vld_table    = [str(s).strip().strip("\x00") for s in chartostring(nc.variables["hdr_vld_table"][:])]

        m2vid = {name.strip(): i for i, name in enumerate(obs_var_list) if name.strip()}

        wsp_vid = m2vid.get(WSP_MNEM)
        hs_vid  = m2vid.get(HS_MNEM)

        if wsp_vid is None:
            print(f"  Skip: missing {WSP_MNEM}")
            continue

        # Load core arrays
        obs_val = nc.variables["obs_val"][:].astype(np.float32)
        obs_vid = nc.variables["obs_vid"][:].astype(np.int64)
        obs_hid = nc.variables["obs_hid"][:].astype(np.int64)
        hdr_vld = nc.variables["hdr_vld"][:].astype(np.int64)
        hdr_lat = nc.variables["hdr_lat"][:].astype(np.float32)
        hdr_lon = nc.variables["hdr_lon"][:].astype(np.float32)

        # Fix 1-based → 0-based if needed
        if obs_vid.min() == 1: obs_vid -= 1
        if obs_hid.min() == 1: obs_hid -= 1
        if hdr_vld.min() == 1: hdr_vld -= 1

        # ── WSPA as MASTER ──────────────────────────────────────────────────
        mask_master = (obs_vid == wsp_vid)
        if not np.any(mask_master):
            print(f"  Skip: no {WSP_MNEM} observations")
            continue

        hid_master = obs_hid[mask_master]
        wsp        = obs_val[mask_master]
        lat        = hdr_lat[hid_master]
        lon        = hdr_lon[hid_master]
        # lon = ((lon + 180) % 360) - 180   # uncomment if you prefer -180/180

        hv_master = hdr_vld[hid_master]

        # Parse times
        time_strs = [vld_table[int(j)].split()[0] if 0 <= int(j) < len(vld_table) else "" for j in hv_master]
        times = np.full(len(hv_master), np.nan, dtype=np.float64)
        for i, ts in enumerate(time_strs):
            if ts:
                try:
                    dt = datetime.strptime(ts, "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)
                    times[i] = dt.timestamp()
                except:
                    pass

        # ── Collocate KBSW (SWH) ────────────────────────────────────────────
        hs = np.full_like(wsp, np.nan, dtype=np.float32)
        if hs_vid is not None:
            mask_hs = (obs_vid == hs_vid)
            if np.any(mask_hs):
                hid_hs = obs_hid[mask_hs]
                val_hs = obs_val[mask_hs]
                hid_to_hs = dict(zip(hid_hs, val_hs))
                for i, h in enumerate(hid_master):
                    hs[i] = hid_to_hs.get(h, np.nan)

        # ── Build DataFrame for debugging ───────────────────────────────────
        df = pd.DataFrame({
            "time_unix": times,
            "time_str": pd.Series(pd.to_datetime(times, unit="s", utc=True, errors="coerce")).dt.strftime("%Y-%m-%d %H:%M:%S"),
            "lat": lat,
            "lon": lon,
            "wsp": wsp,
            "hs": hs,
        })

        # Sort by time (NaNs at end)
        df = df.sort_values("time_unix", na_position="last").reset_index(drop=True)

        # Debug print
        n_total = len(df)
        n_valid_hs = df["hs"].notna().sum()
        print(f"  Total points (WSPA master): {n_total}")
        print(f"  Valid times: {df['time_unix'].notna().sum()}")
        print(f"  Valid wsp  : {df['wsp'].notna().sum()}  (should be nearly all)")
        print(f"  Valid hs   : {n_valid_hs} / {n_total}  ({n_valid_hs / n_total:.1%} matched)")
        print("\n  First 8 rows:")
        print(df.head(8)[["time_str", "lat", "lon", "wsp", "hs"]].to_string(index=False))

        # Optional: save CSV for inspection
        # df.to_csv(f"debug_wspa_master_{ymd}.csv", index=False)

    # ── Write NetCDF ────────────────────────────────────────────────────────
    n = len(df)
    hs_cal  = np.full(n, np.nan, dtype=np.float32)
    wsp_cal = np.full(n, np.nan, dtype=np.float32)

    with Dataset(outfile, "w", format="NETCDF4") as out:
        out.createDimension("time", n)

        vtime = out.createVariable("time", "f8", ("time",))
        vlat  = out.createVariable("lat",  "f4", ("time",))
        vlon  = out.createVariable("lon",  "f4", ("time",))
        vwsp  = out.createVariable("wsp",  "f4", ("time",))
        vhs   = out.createVariable("hs",   "f4", ("time",))
        vwspc = out.createVariable("wsp_cal", "f4", ("time",))
        vhsc  = out.createVariable("hs_cal",  "f4", ("time",))

        vtime.units = "seconds since 1970-01-01 00:00:00 UTC"
        vlat.units  = "degrees_north"
        vlon.units  = "degrees_east"
        vwsp.units  = "m s-1"
        vhs.units   = "m"
        vwspc.units = "m s-1"
        vhsc.units  = "m"

        vwsp.long_name  = "wind_speed"
        vhs.long_name   = "significant_wave_height"
        vwspc.long_name = "calibrated_wind_speed (not available)"
        vhsc.long_name  = "calibrated_significant_wave_height (not available)"

        vtime[:] = df["time_unix"].values
        vlat[:]  = df["lat"].values
        vlon[:]  = df["lon"].values
        vwsp[:]  = df["wsp"].values
        vhs[:]   = df["hs"].values
        vwspc[:] = wsp_cal
        vhsc[:]  = hs_cal

        out.title  = "Track-style altimeter data from MET pb2nc (WSPA master)"
        out.source = "MET pb2nc (preBUFR)"
        out.note   = "hs_cal and wsp_cal filled with NaN"
