"""
Merge daily track-style NetCDF files into monthly NetCDF.

Input daily files:
  Altimeter_JASON3_YYYYMMDD.nc

Output monthly files:
  Altimeter_JASON3_YYYYMM.nc

Variables merged (1D over time):
  time, lat, lon, hs, hs_cal, wsp, wsp_cal
"""

# =========================
# USER-DEFINED SETTINGS
# =========================
WORKDIR = "/work2/noaa/marine/ming.chen/issue107"  # working directory contains WW3-tools
SAT_NAME = "JASON3"
EXTRA_DAYS = 16                                    # Merge extra days for convenience in forecast-day evaluation

IN_DIR = "/work2/noaa/marine/ming.chen/GFS_Retro_Data/processsatdata/out/JASON3" # add option for input directory, default: WORKDIR/processsatdata/out/SATELLITE

# could define merge in specific year and month, set e.g. "202508"
# set to "None" to merge all months found in IN_DIR
TARGET_YYYYMM = "202508"

# ========================

import os
import re
import numpy as np
from netCDF4 import Dataset

# =========================
# DIRECTORIES
# =========================
if "IN_DIR" not in globals():
    IN_DIR = os.path.join(WORKDIR, "processsatdata", "out", SAT_NAME)

if "OUT_DIR" not in globals():
    OUT_DIR = os.path.join(WORKDIR, "processsatdata", "altimeter_monthly")

# =========================
# FUNCTIONS
# =========================
def _list_daily_files(in_dir):
    """
    Return list of (fullpath, yyyymm, yyyymmdd) for daily files in directory.
    """
    pat = re.compile(rf"Altimeter_{SAT_NAME}_(\d{{8}})\.nc$")
    out = []
    for f in os.listdir(in_dir):
        m = pat.match(f)
        if not m:
            continue
        yyyymmdd = m.group(1)
        yyyymm = yyyymmdd[:6]
        out.append((os.path.join(in_dir, f), yyyymm, yyyymmdd))
    out.sort(key=lambda x: x[2])
    return out


def _read_daily(path):
    """
    Read required variables from one daily file.
    Returns dict of numpy arrays.
    """
    with Dataset(path, "r") as nc:
        data = {
            "time":    nc.variables["time"][:].astype(np.float64),
            "lat":     nc.variables["latitude"][:].astype(np.float32),
            "lon":     nc.variables["longitude"][:].astype(np.float32),
            "hs":      nc.variables["hs"][:].astype(np.float32),
            "hs_cal":  nc.variables["hs_cal"][:].astype(np.float32),
            "wsp":     nc.variables["wsp"][:].astype(np.float32),
            "wsp_cal": nc.variables["wsp_cal"][:].astype(np.float32),
        }
        # Grab units if present (safe defaults otherwise)
        units = {
            "time": getattr(nc.variables["time"], "units", "seconds since 1970-01-01 00:00:00 UTC"),
            "lat":  getattr(nc.variables["latitude"], "units", "degrees_north"),
            "lon":  getattr(nc.variables["longitude"], "units", "degrees_east"),
            "hs":   getattr(nc.variables["hs"], "units", "m"),
            "wsp":  getattr(nc.variables["wsp"], "units", "m s-1"),
            "hs_cal": getattr(nc.variables["hs_cal"], "units", "m"),
            "wsp_cal": getattr(nc.variables["wsp_cal"], "units", "m s-1"),
        }
        gattrs = {}
        for k in ("title", "source", "note"):
            if hasattr(nc, k):
                gattrs[k] = getattr(nc, k)
    return data, units, gattrs


def _write_monthly(outfile, merged, units, gattrs):
    """
    Write one monthly NetCDF file.
    """
    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    with Dataset(outfile, "w", format="NETCDF4") as nc:
        nc.createDimension("time", None)
        nc.createDimension("sname", 1)

        vtime = nc.createVariable("time", "f8", ("time",))
        vlat  = nc.createVariable("latitude", "f4", ("time",))
        vlon  = nc.createVariable("longitude", "f4", ("time",))
        vhs   = nc.createVariable("hs", "f4", ("time",))
        vhsc  = nc.createVariable("hs_cal", "f4", ("time",))
        vwsp  = nc.createVariable("wsp", "f4", ("time",))
        vwspc = nc.createVariable("wsp_cal", "f4", ("time",))

        vtime.standard_name = "time"
        vtime.units = "seconds since 1970-01-01 00:00:00"
        vtime.calendar = "standard"
        vtime.axis = "T"

        vlat.units = "degrees_north"
        vlon.units = "degrees_east"

        vhs.long_name = "sigificant_wave_height"
        vhs.units = "m"
        vhsc.long_name = "calibrated_sigificant_wave_height"
        vhsc.units = "m"

        vwsp.long_name = "wind_speed"
        vwsp.units = "m/s"
        vwspc.long_name = "calibrated_wind_speed"
        vwspc.units = "m/s"

        vtime[:] = merged["time"]
        vlat[:]  = merged["lat"]
        vlon[:]  = merged["lon"]
        vhs[:]   = merged["hs"]
        vhsc[:]  = merged["hs_cal"]
        vwsp[:]  = merged["wsp"]
        vwspc[:] = merged["wsp_cal"]

        # Global attrs
        nc.title  = gattrs.get("title", "Monthly merged altimeter track data")
        nc.source = gattrs.get("source", "Daily track NetCDF merged")
        nc.note   = gattrs.get("note", "")
        nc.history = f"Merged daily files into monthly file: {os.path.basename(outfile)}"

def _next_yyyymm(yyyymm):
    """
    Handling month of Dec. when add 16 days for next month in Jan.
    """
    year = int(yyyymm[:4])
    month = int(yyyymm[4:6])
    if month == 12:
        return f"{year + 1}01"
    return f"{year}{month + 1:02d}"

def _select_files_for_month(daily, target_yyyymm, extra_days=16):
    """
    Select files for one output month:
    - all files in target_yyyymm
    - first extra_days files from the next month (extra_days=16 by default)
    """
    next_yyyymm = _next_yyyymm(target_yyyymm)
    selected = []

    for path, yyyymm, yyyymmdd in daily:
        if yyyymm == target_yyyymm:
            selected.append(path)
        elif yyyymm == next_yyyymm:
            day = int(yyyymmdd[6:8])
            if day <= extra_days:
                selected.append(path)

    return selected

# =========================
# MAIN
# =========================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    daily = _list_daily_files(IN_DIR)
    if not daily:
        print("No daily files found.")
        return

    available_months = sorted({yyyymm for _, yyyymm, _ in daily})

    if TARGET_YYYYMM:
        months_to_merge = [TARGET_YYYYMM] if TARGET_YYYYMM in available_months else []
    else:
        months_to_merge = available_months

    if not months_to_merge:
        print("No matching months found.")
        return

    for yyyymm in months_to_merge:

        paths = _select_files_for_month(daily, yyyymm, EXTRA_DAYS)
        if not paths:
            print(f"No files found for output month {yyyymm}")
            continue

        # Read & concatenate
        chunks = {k: [] for k in ("time","lat","lon","hs","hs_cal","wsp","wsp_cal")}
        units = None
        gattrs = None

        for p in paths:
            data, u, ga = _read_daily(p)
            if units is None:
                units = u
            if gattrs is None:
                gattrs = ga
            for k in chunks:
                chunks[k].append(data[k])

        merged = {k: np.concatenate(chunks[k]) for k in chunks}

        # Sort by time
        order = np.argsort(merged["time"])
        for k in merged:
            merged[k] = merged[k][order]

        outfile = os.path.join(OUT_DIR, f"Altimeter_{SAT_NAME}_{yyyymm}.nc")
        _write_monthly(outfile, merged, units, gattrs)

        print(f"Wrote {os.path.basename(outfile)}  n={len(merged['time'])}  ndays={len(paths)}")
        print("Selected file paths:")
        print("====================")

        for p in paths:
            print(p)

if __name__ == "__main__":
    main()

