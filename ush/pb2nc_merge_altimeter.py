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
WORKDIR = "/scratch3/NCEPDEV/marine/Ming.Chen/ursa/wind_eval"  # working directory contains WW3-tools
SAT_NAME = "JASON3"

# could define merge in specific year and month, set e.g. "202508"
# set to "None" to merge all months found in IN_DIR
TARGET_YYYYMM = None

# ========================

import os
import re
import numpy as np
from netCDF4 import Dataset

# =========================
# DIRECTORIES
# =========================
IN_DIR  = os.path.join(WORKDIR, "processsatdata", "pb2nc_altimeter")
OUT_DIR = os.path.join(WORKDIR, "processsatdata", "pb2nc_altimeter_monthly")

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
            "lat":     nc.variables["lat"][:].astype(np.float32),
            "lon":     nc.variables["lon"][:].astype(np.float32),
            "hs":      nc.variables["hs"][:].astype(np.float32),
            "hs_cal":  nc.variables["hs_cal"][:].astype(np.float32),
            "wsp":     nc.variables["wsp"][:].astype(np.float32),
            "wsp_cal": nc.variables["wsp_cal"][:].astype(np.float32),
        }
        # Grab units if present (safe defaults otherwise)
        units = {
            "time": getattr(nc.variables["time"], "units", "seconds since 1970-01-01 00:00:00 UTC"),
            "lat":  getattr(nc.variables["lat"], "units", "degrees_north"),
            "lon":  getattr(nc.variables["lon"], "units", "degrees_east"),
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

# =========================
# MAIN
# =========================
def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    daily = _list_daily_files(IN_DIR)
    if not daily:
        print("No daily files found.")
        return

    # Group by month
    months = {}
    for path, yyyymm, yyyymmdd in daily:
        if TARGET_YYYYMM and yyyymm != TARGET_YYYYMM:
            continue
        months.setdefault(yyyymm, []).append(path)

    if not months:
        print("No matching months found.")
        return

    for yyyymm, paths in sorted(months.items()):
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


if __name__ == "__main__":
    main()

