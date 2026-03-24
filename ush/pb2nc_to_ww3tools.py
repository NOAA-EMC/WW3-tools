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
WORKDIR = "/scratch3/NCEPDEV/marine/Ming.Chen/ursa/preBUFR_filter/"            # working directory contains WW3-tools

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
        obs_var_list = [str(s).strip().strip("\x00") for s in chartostring(nc.variables["obs_var"][:])]        # data name
        vld_table    = [str(s).strip().strip("\x00") for s in chartostring(nc.variables["hdr_vld_table"][:])]  # timestamps

        m2vid = {name.strip(): i for i, name in enumerate(obs_var_list) if name.strip()}                       # add indexes in data name

        wsp_vid = m2vid.get(WSP_MNEM)    # number represents WSP (wind speed)
        hs_vid  = m2vid.get(HS_MNEM)     # number represents HS  (significant wave height)

        if wsp_vid is None:
            print(f"  Skip: missing {WSP_MNEM}")
            continue

        # Load core arrays
        obs_val = nc.variables["obs_val"][:].astype(np.float32)   # mixed data
        obs_vid = nc.variables["obs_vid"][:].astype(np.int64)     # variable ID
        obs_hid = nc.variables["obs_hid"][:].astype(np.int64)     # connection key to connect val/vid to hdr_lat/hdr_lon/hdr_vld
        hdr_vld = nc.variables["hdr_vld"][:].astype(np.int64)     # time index for lat and lon
        hdr_lat = nc.variables["hdr_lat"][:].astype(np.float32)   # latitude
        hdr_lon = nc.variables["hdr_lon"][:].astype(np.float32)   # longitude

        # Fix 1-based → 0-based if needed
        if obs_vid.min() == 1: obs_vid -= 1
        if obs_hid.min() == 1: obs_hid -= 1
        if hdr_vld.min() == 1: hdr_vld -= 1

        # create WSP table
        mask_wsp = (obs_vid == wsp_vid)                           # find indexes of WSP for mapping WSP data in obs_val
        if not np.any(mask_wsp):
            print(f"  Skip: no {WSP_MNEM} observations")
            continue

        hid_wsp = obs_hid[mask_wsp]                               # connection keys of WSP for mapping lat, lon, and time
        # create dataframe for WSP
        df_wsp = pd.DataFrame({
            "hid": hid_wsp,
            "lat": hdr_lat[hid_wsp],                              # lat using connection key matching WSP
            "lon": hdr_lon[hid_wsp],                              # lon using connection key matching WSP
            "vld": hdr_vld[hid_wsp],                              # Store the time index here
            "wsp": obs_val[mask_wsp]                              # WSP obs values
        })

        # Create the HS table separately to avoid losing data when using WSP as the master variable
        mask_hs = (obs_vid == hs_vid)
        if not np.any(mask_hs):
            print(f"  Skip: no {HS_MNEM} observations")
            continue

        hid_hs  = obs_hid[mask_hs]                                # connection keys of Hs for mapping lat, lon, and time
        # create dataframe for HS
        df_hs = pd.DataFrame({
            "hid": hid_hs,
            "lat": hdr_lat[hid_hs],                              # lat using connection key matching Hs
            "lon": hdr_lon[hid_hs],                              # lon using connection key matching Hs
            "vld": hdr_vld[hid_hs],                              # Store the time index here
            "hs":  obs_val[mask_hs]                              # Hs obs values
        })

        # merge WSP and HS
        df = pd.merge(df_wsp, df_hs, on=["hid", "lat", "lon", "vld"], how="outer")

        # fill in the timestamps based on vld
        time_strs = [vld_table[int(j)].split()[0] if 0 <= int(j) < len(vld_table) else "" for j in df["vld"]]
        times = np.full(len(df), np.nan, dtype=np.float64)       # convert yyyymmdd_hhmmss to unix timestamps
        for i, ts in enumerate(time_strs):
            if ts:
                try:
                    dt = datetime.strptime(ts, "%Y%m%d_%H%M%S").replace(tzinfo=timezone.utc)
                    times[i] = dt.timestamp()
                except:
                    pass

        df["time_unix"] = times
        df["time_str"] = pd.to_datetime(df["time_unix"], unit="s", utc=True, errors="coerce").dt.strftime("%Y-%m-%d %H:%M:%S")

        df = df.drop(columns=["vld"])[["time_unix", "time_str", "lat", "lon", "wsp", "hs"]]

        # sort the table by time (not necessary to move NaNs at end)
        df = df.sort_values(by="time_unix", ascending=True)
        df = df.reset_index(drop=True)

        # QC - Physical Bounds Clipping: Wind Speed (m/s) 0.5 - 100; SWH (m) 0 - 50
        df.loc[(df['wsp'] < 0.5) | (df['wsp'] > 100), 'wsp'] = np.nan   # values below 0.5 m/s sensor noise floor
        df.loc[(df['hs'] < 0) | (df['hs'] > 50), 'hs']     = np.nan

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
        vlat  = out.createVariable("latitude",  "f4", ("time",))
        vlon  = out.createVariable("longitude",  "f4", ("time",))
        vwsp  = out.createVariable("wsp",  "f4", ("time",))
        vhs   = out.createVariable("hs",   "f4", ("time",))
        vwspc = out.createVariable("wsp_cal", "f4", ("time",))
        vhsc  = out.createVariable("hs_cal",  "f4", ("time",))

        vtime.units    = "seconds since 1970-01-01 00:00:00 UTC"
        vtime.calendar = "standard"
        vlat.units     = "degrees_north"
        vlon.units     = "degrees_east"
        vwsp.units     = "m s-1"
        vhs.units      = "m"
        vwspc.units    = "m s-1"
        vhsc.units     = "m"

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
