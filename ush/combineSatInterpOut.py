"""
combineSatInterpOut.py

PURPOSE:

    Create combined NetCDF files for a specified time period using
interpolated model-to-satellite matchup data.

    The script supports:
        1. one all-in-one combined output file
        2. separate combined output files by forecast day

    These combined files are intended to simplify downstream
post-processing, plotting, and statistical evaluation.

USAGE:
    Example:
        python combineSatInterpOut.py \
            -models retrov17_01 GFSv16 \
            -WORKDIR /scratch3/NCEPDEV/marine/Ming.Chen/ursa/ww3tools \
            -satellites JASON3 CRYOSAT2 SARAL SENTINEL3A \
            -startdate 2024111512 \
            -enddate   2024121718 \
            -interval_hours 6 \
            -max_forecast_day 16 \
            -force_season winter \
            -selected_years 2024
            -INPUTDIR_BASE /work2/noaa/marine/ming.chen/GFS_Retro_Data/data/outinterp

        Disable outputs:
            --no_output_all_in_one
            --no_output_per_day
OUTPUT:
    Output location:
        {WORKDIR}/processsatdata/outcombine/{satellite}/

    Output filename format:
        1. All-in-one combined file
            combined_all_{model}_{season}_{satellite}.nc
        2. Per-forecast-day files
            combined_day{DD}_{model}_{season}_{satellite}.nc

Author and DATE:
 03/12/2026: Ming Chen (ming.chen1@noaa.gov)

"""

import numpy as np
import netCDF4 as nc
import datetime as dt
import os
import xarray as xr
import glob
import argparse

def _parse_yyyymmddhh(s):
    try:
        return dt.datetime.strptime(s, "%Y%m%d%H")
    except ValueError as e:
        raise ValueError(f"Invalid datetime '{s}'. Expected YYYYMMDDHH, e.g. 2024111512") from e

def parse_args():
    parser = argparse.ArgumentParser(
        description="Create combined NetCDF files for easier post processing."
    )

    # Required
    parser.add_argument("-models", nargs="+", required=True,
                        help="Models list, e.g. retrov17_01 GFSv16")
    parser.add_argument("-WORKDIR", required=True,
                        help="Base WORKDIR containing processsatdata/")
    parser.add_argument("-satellites", nargs="+", required=True,
                        help="Satellites list, e.g. JASON3 CRYOSAT2 SARAL SENTINEL3A")
    parser.add_argument("-startdate", required=True,
                        help="Start cycle YYYYMMDDHH")
    parser.add_argument("-enddate", required=True,
                        help="End cycle YYYYMMDDHH")

    # Optional knobs
    parser.add_argument("-interval_hours", type=int, default=6,
                        help="Cycle interval hours (default: 6)")
    parser.add_argument("-max_forecast_day", type=int, default=16,
                        help="Max forecast day (default: 16)")
    parser.add_argument("-force_season", default=None,
                        help="Force season name (winter/summer/hurricane/other) or None (default: auto)")
    parser.add_argument("-INPUTDIR_BASE", default=None,
                        help="Optional input directory base. Default: WORKDIR/processsatdata/outinterp")

    # Outputs: keep original behavior (both True) unless disabled
    parser.add_argument("--no_output_all_in_one", action="store_true",
                        help="Disable all-in-one output file")
    parser.add_argument("--no_output_per_day", action="store_true",
                        help="Disable per-day output files")

    args = parser.parse_args()

    # Assign globals used by the original code
    global models, WORKDIR, satellites
    global startdate, enddate, interval_hours, max_forecast_day
    global force_season
    global output_all_in_one, output_per_day
    global INPUTDIR_BASE, OUTDIR

    models = args.models
    WORKDIR = args.WORKDIR
    satellites = args.satellites

    startdate = _parse_yyyymmddhh(args.startdate)
    enddate = _parse_yyyymmddhh(args.enddate)
    if enddate < startdate:
        raise ValueError("enddate must be >= startdate")

    interval_hours = args.interval_hours
    max_forecast_day = args.max_forecast_day

    # allow passing literal "None" from CLI
    if isinstance(args.force_season, str) and args.force_season.lower() == "none":
        force_season = None
    else:
        force_season = args.force_season

    output_all_in_one = not args.no_output_all_in_one
    output_per_day = not args.no_output_per_day

    if args.INPUTDIR_BASE:
        INPUTDIR_BASE = args.INPUTDIR_BASE
    else:
        INPUTDIR_BASE = os.path.join(WORKDIR, "processsatdata", "outinterp")

    OUTDIR = os.path.join(WORKDIR, "processsatdata", "outcombine")

def determine_season(dt_obj):
    month = dt_obj.month
    day   = dt_obj.day

    # Atlantic hurricane season (NOAA): June 1 – Nov 30
    if (month == 6 and day >= 1) or month in [7,8,9,10] or (month == 11 and day <= 30):
        return "hurricane"

    # Meteorological winter: Dec–Feb
    if month in [12, 1, 2]:
        return "winter"

    # Meteorological summer: Jun–Aug
    if month in [6, 7, 8]:
        return "summer"

    return "other"

def get_inputdir(model):
    return os.path.join(INPUTDIR_BASE, model)

def get_grids_for_model(model):
    if model == "multi1":
        return ['global.0p50', 'alaska.0p16', 'atlocn.0p16', 'epacif.0p16', 'wcoast.0p16',
                'alaska.0p06', 'atlocn.0p06', 'wcoast.0p06']
    elif model == "GFSv16":
        return ['global.0p25']
    elif model == "retrov17_01":
        return ['global.0p25']
    else:
        return ['global.0p25']

def get_max_forecast_days(model, season):
    return max_forecast_day

def main():

    if not os.path.isdir(OUTDIR):
        os.makedirs(OUTDIR)

    now = startdate
    all_cycles = []
    while now <= enddate:
        all_cycles.append(now)
        now += dt.timedelta(hours=interval_hours)

    cycles_to_process = all_cycles

    if not cycles_to_process:
        print("No cycles found between startdate and enddate.")
        return

    print(
        f"Processing {len(cycles_to_process)} cycles "
        f"from {startdate:%Y%m%d%H} to {enddate:%Y%m%d%H} "
        f"(interval {interval_hours}h)"
    )

    for model in models:
        print(f"\n=== Model: {model} ===")
        grids = get_grids_for_model(model)
        inputdir = get_inputdir(model)

        if not os.path.isdir(inputdir):
            print(f"  Directory not found: {inputdir}")
            continue

        cycles_by_season = {}
        for cycle_dt in cycles_to_process:
            season = force_season if force_season else determine_season(cycle_dt)
            cycles_by_season.setdefault(season, []).append(cycle_dt.strftime('%Y%m%d%H'))

        for season_name, date_strs in cycles_by_season.items():
            print(f"  → {season_name}: {len(date_strs)} cycles")

            endday = get_max_forecast_days(model, season_name)

            for sat in satellites:
                data = collect_data(model, sat, date_strs, grids, inputdir, season_name)
                if not data.get('time'):
                    continue

                arr = {k: np.array(v) for k, v in data.items()}
                adjust_wind(model, arr['fhrs'], arr['model_wnd'])
                write_outputs(arr, model, sat, season_name, endday, year=None)

def collect_data(model, sat, date_list, grids, inputdir, season_name):
    collected = {
        'time':[], 'lats':[], 'lons':[], 'fhrs':[], 'fhrsall':[],
        'obs_hs':[], 'obs_hs_cal':[], 'obs_wnd':[], 'obs_wnd_cal':[],
        'model_hs':[], 'model_wnd':[]
    }
    for dstr in date_list:
        cycle = process_one_cycle(model, sat, dstr, grids, inputdir, season_name)
        if cycle:
            for k in collected:
                collected[k].extend(cycle[k])
    return collected


def process_one_cycle(model, sat, date_str, grids, inputdir, season_name):
    base = None
    for g_idx, grid in enumerate(grids):
        # Adjust filename pattern if your files include season name or not
        # Here assuming no season in filename:
        fname = f"{model}_{grid}_{date_str}_{sat}.nc"

        path = os.path.join(inputdir, fname)
        if not os.path.exists(path):
            continue

        with nc.Dataset(path) as ds:
            t   = np.array(ds.variables['time'][:])
            fh  = np.array(ds.variables['fcst_hr'][:])
            la  = np.array(ds.variables['latitude'][:])
            lo  = np.array(ds.variables['longitude'][:])
            ohs = np.array(ds.variables['obs_hs'][:])
            ohc = np.array(ds.variables['obs_hs_cal'][:])
            ow  = np.array(ds.variables['obs_wnd'][:])
            owc = np.array(ds.variables['obs_wnd_cal'][:])
            mhs = np.array(ds.variables['model_hs'][:])
            mw  = np.array(ds.variables['model_wnd'][:])

            ict = ds.getncattr('initial_condition_time')
            fha = (t - ict) / 3600.0

            if g_idx == 0:
                base = {
                    'time':t, 'lats':la, 'lons':lo,
                    'fhrs':fh, 'fhrsall':fha,
                    'obs_hs':ohs, 'obs_hs_cal':ohc,
                    'obs_wnd':ow, 'obs_wnd_cal':owc,
                    'model_hs':mhs, 'model_wnd':mw
                }
            else:
                if base is None: continue
                base['model_hs'] = np.where(~np.isnan(mhs), mhs, base['model_hs'])
                base['model_wnd'] = np.where(~np.isnan(mw), mw, base['model_wnd'])

    return base


def adjust_wind(model, fhrs, model_wnd):
    if model in ["HR1", "HR2"]:
        model_wnd[fhrs < 3] = np.nan
    elif model in ["HR3a", "HR3b"]:
        model_wnd[fhrs < 1] = np.nan


def write_outputs(arr, model, sat, season, endday, year):
    year_str = f"_{year}" if year else ""

    def save(filename_suffix, t, la, lo, ohs, ohc, ow, owc, mhs, mw, fh):
        # All variables as double (float64)
        ds = xr.Dataset({
            'time': xr.DataArray(
                t.astype(np.float64),
                dims=['time'],
                name='time',
                attrs={
                    'standard_name': 'time',
                    'units': 'seconds since 1970-01-01 00:00:00',
                    'calendar': 'standard',
                    'axis': 'T'
                }
            ),
            'latitude': xr.DataArray(
                la.astype(np.float64),
                dims=['time'],
                attrs={'units': 'degree_north'}
            ),
            'longitude': xr.DataArray(
                lo.astype(np.float64),
                dims=['time'],
                attrs={'units': 'degree_east'}
            ),
            'model_hs': xr.DataArray(
                mhs.astype(np.float64),
                dims=['time'],
                attrs={'units': 'm'}
            ),
            'model_wnd': xr.DataArray(
                mw.astype(np.float64),
                dims=['time'],
                attrs={'units': 'm'}
            ),
            'obs_hs': xr.DataArray(
                ohs.astype(np.float64),
                dims=['time'],
                attrs={'units': 'm'}
            ),
            'obs_hs_cal': xr.DataArray(
                ohc.astype(np.float64),
                dims=['time'],
                attrs={'units': 'm'}
            ),
            'obs_wnd': xr.DataArray(
                ow.astype(np.float64),
                dims=['time'],
                attrs={'units': 'm/s'}
            ),
            'obs_wnd_cal': xr.DataArray(
                owc.astype(np.float64),
                dims=['time'],
                attrs={'units': 'm/s'}
            ),
            'fcst_hr': xr.DataArray(
                fh.astype(np.float64),
                dims=['time'],
                attrs={
                    'description': 'Forecast hour relative to initial condition time',
                    'units': 'hours'
                }
            ),
        })

        # Explicitly set _FillValue = NaN for all variables
        encoding = {
            var: {'_FillValue': np.nan, 'dtype': 'float64'}
            for var in ds.variables
        }

        # Global attributes
        ds.attrs['satellite_name'] = sat
        ds.attrs['model_name'] = model

        fullpath = os.path.join(OUTDIR, filename_suffix)
        ds.to_netcdf(
            fullpath,
            format='NETCDF4',
            encoding=encoding,
            engine='netcdf4'
        )
        print(f"    Wrote: {os.path.basename(fullpath)}  ({len(t)} points)")

    if output_all_in_one:
        fname = f"combined_{model}_{season}_{sat}.nc"
        save(fname, arr['time'], arr['lats'], arr['lons'],
             arr['obs_hs'], arr['obs_hs_cal'], arr['obs_wnd'], arr['obs_wnd_cal'],
             arr['model_hs'], arr['model_wnd'], arr['fhrs'])

    if output_per_day:
        day0 = 0
        day = 1
        while day <= endday:
            idx = (arr['fhrsall'] > day0*24) & (arr['fhrsall'] <= day*24)
            if np.any(idx):
                fname = f"combined_day{day:02d}_{model}_{season}_{sat}.nc"
                save(fname,
                     arr['time'][idx], arr['lats'][idx], arr['lons'][idx],
                     arr['obs_hs'][idx], arr['obs_hs_cal'][idx],
                     arr['obs_wnd'][idx], arr['obs_wnd_cal'][idx],
                     arr['model_hs'][idx], arr['model_wnd'][idx], arr['fhrs'][idx])
            day0 = day
            day += 1


if __name__ == '__main__':
    parse_args()
    main()
