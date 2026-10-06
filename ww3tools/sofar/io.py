import subprocess
import os
import glob
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr
import ncepbufr

def convert_sofar_bufr_to_nc(
    bufr_dir, 
    output_dir=None,
    start_date=None,
    end_date=None,
    pattern="sofar_b001_xx121_*"
):
    """
    Converts BUFR marine observation files to CF-1.9 compliant NetCDF files.

    Parameters
    ----------
    bufr_dir : str or Path
        Directory containing BUFR files.
    output_dir : str or Path, optional
        Directory where output .nc files will be saved. Defaults to bufr_dir.
    start_date : str, datetime, or pd.Timestamp, optional
        Start of filtering date range (inclusive), e.g., '2026-01-01'.
    end_date : str, datetime, or pd.Timestamp, optional
        End of filtering date range (inclusive), e.g., '2026-01-05'.
    pattern : str, optional
        Filename glob pattern to match BUFR files. Defaults to 'sofar_b001_xx121_*'.

    Returns
    -------
    list of Path
        List of generated NetCDF file paths.
    """
    bufr_path = Path(bufr_dir)
    out_path = Path(output_dir) if output_dir else bufr_path
    out_path.mkdir(parents=True, exist_ok=True)

    t_start = pd.to_datetime(start_date) if start_date else None
    t_end = pd.to_datetime(end_date) if end_date else None

    file_list = sorted([
        f for f in bufr_path.glob(pattern) 
        if not f.name.endswith('.nc') and f.is_file()
    ])

    if not file_list:
        print(f"No matching BUFR files found in {bufr_dir} with pattern '{pattern}'.")
        return []

    print(f"Found {len(file_list)} BUFR files to process.")

    string_vars = ['LSTN']
    numeric_vars = [
        'CLATH', 'CLONH', 'YEAR', 'MNTH', 'DAYS', 'HOUR', 'MINU',
        'SGWH', 'WDIR', 'WSPD', 'DSDW', 'DDWC', 'AVWP', 'SPWP'
    ]
    variables = string_vars + numeric_vars
    wave_keys = ['SGWH', 'WDIR', 'WSPD', 'DSDW', 'DDWC', 'AVWP', 'SPWP']

    # BUFR Mnemonic to CF Standard Variable Name Mapping
    rename_map = {
        'CLATH': 'latitude',
        'CLONH': 'longitude',
        'LSTN': 'platform_id',
        'SGWH': 'sea_surface_wave_significant_height',
        'WDIR': 'wind_from_direction',
        'WSPD': 'wind_speed',
        'DSDW': 'sea_surface_wave_directional_spread',
        'DDWC': 'sea_surface_wave_from_direction_at_variance_spectral_peak',
        'AVWP': 'sea_surface_wave_mean_period',
        'SPWP': 'sea_surface_wave_period_at_variance_spectral_peak',
    }

    # CF-1.9 Variable Metadata Definitions
    var_metadata = {
        'latitude': {
            'standard_name': 'latitude',
            'long_name': 'latitude',
            'units': 'degrees_north',
            'axis': 'Y'
        },
        'longitude': {
            'standard_name': 'longitude',
            'long_name': 'longitude',
            'units': 'degrees_east',
            'axis': 'X'
        },
        'time': {
            'standard_name': 'time',
            'long_name': 'time of observation',
            'axis': 'T'
        },
        'platform_id': {
            'long_name': 'platform identifier',
            'cf_role': 'trajectory_id'
        },
        'sea_surface_wave_significant_height': {
            'standard_name': 'sea_surface_wave_significant_height',
            'long_name': 'Significant wave height',
            'units': 'm'
        },
        'wind_from_direction': {
            'standard_name': 'wind_from_direction',
            'long_name': 'Wind direction',
            'units': 'degree'
        },
        'wind_speed': {
            'standard_name': 'wind_speed',
            'long_name': 'Wind speed',
            'units': 'm s-1'
        },
        'sea_surface_wave_directional_spread': {
            'standard_name': 'sea_surface_wave_directional_spread',
            'long_name': 'Directional spread of surface waves',
            'units': 'degree'
        },
        'sea_surface_wave_from_direction_at_variance_spectral_peak': {
            'standard_name': 'sea_surface_wave_from_direction_at_variance_spectral_peak',
            'long_name': 'Dominant wave direction at spectral peak',
            'units': 'degree'
        },
        'sea_surface_wave_mean_period': {
            'standard_name': 'sea_surface_wave_mean_period',
            'long_name': 'Average wave period',
            'units': 's'
        },
        'sea_surface_wave_period_at_variance_spectral_peak': {
            'standard_name': 'sea_surface_wave_period_at_variance_spectral_peak',
            'long_name': 'Peak wave period',
            'units': 's'
        }
    }

    converted_files = []

    for filepath in file_list:
        records = []
        try:
            bufr = ncepbufr.open(str(filepath))
            while bufr.advance() == 0:
                while bufr.load_subset() == 0:
                    rec = {}
                    for var in variables:
                        try:
                            if var in string_vars:
                                raw = bufr.read_subset(var)
                                text = raw.view('S8').tobytes().decode('ascii', errors='ignore').replace('\x00', '').strip()
                                rec[var] = text if text else None
                            else:
                                val = bufr.read_subset(var)
                                v = val.flat[0] if val.size > 0 else np.nan
                                
                                if np.ma.is_masked(v) or v >= 1e9:
                                    rec[var] = np.nan
                                else:
                                    rec[var] = float(v)
                        except Exception:
                            rec[var] = None if var in string_vars else np.nan
                    
                    if any(pd.notna(rec[k]) for k in wave_keys if k in rec):
                        records.append(rec)
            bufr.close()
        except Exception as e:
            print(f"Error reading BUFR file {filepath.name}: {e}")
            continue

        if not records:
            print(f"  [{filepath.name}] No valid records extracted. Skipping.")
            continue

        df = pd.DataFrame(records)

        # Assemble Datetime column
        time_cols = ['YEAR', 'MNTH', 'DAYS', 'HOUR']
        if all(col in df.columns for col in time_cols):
            df['time'] = pd.to_datetime({
                'year': df['YEAR'].fillna(1970).astype(int),
                'month': df['MNTH'].fillna(1).astype(int),
                'day': df['DAYS'].fillna(1).astype(int),
                'hour': df['HOUR'].fillna(0).astype(int),
                'minute': df.get('MINU', 0).fillna(0).astype(int)
            })
            df = df.drop(columns=['YEAR', 'MNTH', 'DAYS', 'HOUR', 'MINU'], errors='ignore')

        # Filter by date range if provided
        if 'time' in df.columns:
            if t_start is not None:
                df = df[df['time'] >= t_start]
            if t_end is not None:
                df = df[df['time'] <= t_end]

        if df.empty:
            print(f"  [{filepath.name}] No records fell within date range [{start_date} to {end_date}]. Skipping.")
            continue

        # Reset index to create 1D observation dimension 'obs'
        df = df.reset_index(drop=True)
        df.index.name = 'obs'
        ds = df.to_xarray()

        # Rename BUFR mnemonics to CF standard variable names
        rename_dict = {k: v for k, v in rename_map.items() if k in ds}
        if rename_dict:
            ds = ds.rename(rename_dict)

        # Set time as the coordinate
        ds = ds.set_coords('time')

        # Apply CF-1.9 Global Attributes
        now_str = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        ds.attrs['Conventions'] = 'CF-1.9'
        ds.attrs['featureType'] = 'trajectory'
        ds.attrs['title'] = 'Surface Marine Observation Data Converted from NCEP BUFR'
        ds.attrs['source'] = f'NCEP BUFR observation file: {filepath.name}'
        ds.attrs['history'] = f'{now_str}: Converted to CF-1.9 NetCDF via ww3-tools pipeline.'

        # # Assign CF Variable Attributes & Link Coordinates
        # coord_attr_str = 'time latitude longitude'
        # non_data_vars = {'time', 'latitude', 'longitude', 'platform_id'}

        # for vname, attrs in var_metadata.items():
        #     if vname in ds:
        #         ds[vname].attrs.update(attrs)
        #         if vname not in non_data_vars:
        #             ds[vname].attrs['coordinates'] = coord_attr_str

        # Define explicit CF time encoding
        encoding = {}
        if 'time' in ds:
            encoding['time'] = {
                'units': 'seconds since 1970-01-01 00:00:00',
                'calendar': 'gregorian'
            }

        # Write to NetCDF file
        nc_filename = out_path / f"{filepath.name}.nc"
        ds.to_netcdf(nc_filename, encoding=encoding)
        converted_files.append(nc_filename)
        print(f"  [{filepath.name}] Saved CF-1.9 compliant NetCDF ({len(df)} records) -> {nc_filename.name}")

    print(f"\nDone! Successfully converted {len(converted_files)} files.")
    return converted_files