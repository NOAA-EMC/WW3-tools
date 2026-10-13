import subprocess
import os
import glob
import datetime
from pathlib import Path
from datetime import datetime, timedelta
import sys

import numpy as np
import pandas as pd
import netCDF4 as nc
import xarray as xr
import ncepbufr
import dask
from dask.distributed import Client, LocalCluster

# Setup dask cluster for parallelization
def setup_dask_cluster():
    """
    Initializes a Dask cluster constrained by Slurm's allocated CPUs.
    """
    slurm_cpus = int(os.environ.get('SLURM_CPUS_PER_TASK', 4))
    
    cluster = LocalCluster(
        n_workers=slurm_cpus,     
        threads_per_worker=1,     
        memory_limit='auto'       
    )
    
    client = Client(cluster)
    print(f"Initialized Dask cluster with {slurm_cpus} workers.")
    return client

# Fetch the Sofar data from the HPSS archive
def check_and_fetch_sofar_hpss(
    start_date_str, 
    end_date_str, 
    work_dir,
    data_dir="bufr_data",
    slurm_params=None
):
    """
    Checks for daily Sofar buoy data from start_date to end_date + 16 days.
    If files are missing, submits an HPSS extraction job for those specific 
    files and terminates the current process.
    """
    target_dir = os.path.join(work_dir, data_dir)
    os.makedirs(target_dir, exist_ok=True)
    
    start_date = datetime.strptime(start_date_str, "%Y-%m-%d")
    end_date = datetime.strptime(end_date_str, "%Y-%m-%d") + timedelta(days=16)
    
    num_days = (end_date - start_date).days + 1
    date_list = [start_date + timedelta(days=x) for x in range(num_days)]
    
    missing_dates = []
    archive_base = "/NCEPPROD/hpssprod/runhistory"

    # First loop: Check all required dates to see if the file is present
    for dt in date_list:
        yy = dt.strftime("%Y")
        mm = dt.strftime("%m")
        dd = dt.strftime("%d")
        target_file = f"sofar_b001_xx121_{yy}{mm}{dd}"
        
        if not os.path.exists(os.path.join(target_dir, target_file)):
            missing_dates.append(dt)

    # Evaluate if the missing_dates list is empty
    if not missing_dates:
        print("Success: All required Sofar HPSS files are present in the data directory.")
        return True

    # If files are missing, configure parameters and build extraction commands
    print(f"Missing {len(missing_dates)} files. Generating and submitting HPSS retrieval script...")
    
    slurm_params = slurm_params or {
        "ntasks": 1,
        "queue": "batch",
        "partition": "u1-service",
        "walltime": "01:00:00",
        "account": "marine-cpu",
        "job_name": "get_sofar_hpss"
    }

    extraction_commands = []
    
    # Second loop: Build extraction commands ONLY for the missing dates
    for dt in missing_dates:
        yy = dt.strftime("%Y")
        mm = dt.strftime("%m")
        dd = dt.strftime("%d")
        target_file = f"sofar_b001_xx121_{yy}{mm}{dd}"
        file_atm = f"{archive_base}/rh{yy}/{yy}{mm}/{yy}{mm}{dd}/dcom_{yy}{mm}{dd}.tar"
        
        extraction_commands.append(f"""
echo "Processing {yy}-{mm}-{dd}..."
htar -P -xvf {file_atm} ./b001
if [ -f "b001/xx121" ]; then
    mv b001/xx121 "{target_file}"
    echo "Successfully extracted {target_file}"
else
    echo "Warning: Target xx121 not found in {file_atm}"
    missing_files+=("{target_file}")
fi
rm -rf b001
""")

    job_name = slurm_params.get("job_name", "get_sofar_hpss")
    script_content = f"""#!/bin/sh --login
#SBATCH --ntasks={slurm_params.get("ntasks", 1)}
#SBATCH -q {slurm_params.get("queue", "batch")}
#SBATCH --partition={slurm_params.get("partition", "u1-service")}
#SBATCH -t {slurm_params.get("walltime", "01:00:00")}
#SBATCH -A {slurm_params.get("account", "marine-cpu")}
#SBATCH -J {job_name}
#SBATCH -o {target_dir}/{job_name}_%j.out

cd {target_dir}

missing_files=()

{"".join(extraction_commands)}

if [ ${{#missing_files[@]}} -ne 0 ]; then
    echo "=========================================="
    echo "The following files could NOT be extracted:"
    for file in "${{missing_files[@]}}"; do
        echo " - $file"
    done
    echo "=========================================="
else
    echo "All requested files were successfully extracted."
fi
"""

    script_path = os.path.join(target_dir, "submit_hpss.sh")
    with open(script_path, "w") as f:
        f.write(script_content)

    subprocess.run(["sbatch", script_path], check=True)
    print("HPSS retrieval job submitted to service node. Terminating current process.")
    sys.exit(0)

# Check that model data exists
def check_model_data(model_dir, start_date_str, end_date_str, pattern_prefix="gfs."):
    """
    Verifies that model data directories or files exist for every day between 
    the start and end dates.
    If any data is missing, it prints the missing dates/paths and terminates the job.
    """
    # Check if the base model directory actually exists
    if not os.path.exists(model_dir):
        print(f"Error: The base model directory '{model_dir}' does not exist.")
        print("Terminating the current job.")
        sys.exit(0)

    # Generate list of required dates
    start_date = datetime.strptime(start_date_str, "%Y-%m-%d")
    end_date = datetime.strptime(end_date_str, "%Y-%m-%d")
    
    num_days = (end_date - start_date).days + 1
    date_list = [start_date + timedelta(days=x) for x in range(num_days)]
    
    missing_items = []

    # Check for existing model directories or files for each date
    for dt in date_list:
        date_formatted = dt.strftime("%Y%m%d")
        target_name = f"{pattern_prefix}{date_formatted}"
        target_path = os.path.join(model_dir, target_name)
        
        if not os.path.exists(target_path):
            missing_items.append(target_path)

    # Evaluate if any items are missing
    if not missing_items:
        print(f"Success: All required model data from {start_date_str} to {end_date_str} is present.")
        print("Continuing with the analysis...")
        return True

    # Missing data path: print exactly what is missing and kill the job
    print("Error: The following model data directories or files are missing:")
    for item in missing_items:
        print(f" - {item}")
    
    print("Terminating the current job. Please ensure the model outputs are properly staged.")
    sys.exit(0)


# Convert the sofar data from bufr to netCDF
def convert_sofar_bufr_to_nc(
    work_dir,
    bufr_subdir="bufr_data", 
    nc_subdir="nc_data",
    start_date=None,
    end_date=None,
    pattern="sofar_b001_xx121_*"
):
    """
    Converts BUFR marine observation files to CF-compliant NetCDF files.
    Skips conversion if the target NetCDF file already exists.
    Directories are built relative to the provided work_dir.
    Appends a sixteen-day buffer to the provided end date.
    """
    # Establish relative paths based on the working directory
    bufr_path = Path(work_dir) / bufr_subdir
    out_path = Path(work_dir) / nc_subdir
    out_path.mkdir(parents=True, exist_ok=True)

    t_start = pd.to_datetime(start_date) if start_date else None
    
    # Add a sixteen-day buffer to the end date to capture trailing data
    t_end = pd.to_datetime(end_date) + pd.Timedelta(days=16) if end_date else None

    file_list = sorted([
        f for f in bufr_path.glob(pattern) 
        if not f.name.endswith('.nc') and f.is_file()
    ])

    if not file_list:
        print(f"No matching BUFR files found in {bufr_path} with pattern '{pattern}'.")
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

    converted_files = []

    for filepath in file_list:
        nc_filename = out_path / f"{filepath.name}.nc"
        
        # Check if the NetCDF file already exists to bypass conversion
        if nc_filename.exists():
            print(f"  [{filepath.name}] NetCDF already exists -> {nc_filename.name}. Skipping conversion.")
            converted_files.append(nc_filename)
            continue

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
            print(f"  [{filepath.name}] No records fell within date range [{start_date} to extended end date]. Skipping.")
            continue

        # Reset index to create observation dimension
        df = df.reset_index(drop=True)
        df.index.name = 'obs'
        ds = df.to_xarray()

        # Rename BUFR mnemonics to CF standard variable names
        rename_dict = {k: v for k, v in rename_map.items() if k in ds}
        if rename_dict:
            ds = ds.rename(rename_dict)

        # Set time as the coordinate
        ds = ds.set_coords('time')

        # Apply Global Attributes
        now_str = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        ds.attrs['Conventions'] = 'CF-1.9'
        ds.attrs['featureType'] = 'trajectory'
        ds.attrs['title'] = 'Surface Marine Observation Data Converted from NCEP BUFR'
        ds.attrs['source'] = f'NCEP BUFR observation file: {filepath.name}'
        ds.attrs['history'] = f'{now_str}: Converted to CF-compliant NetCDF via ww3-tools pipeline.'

        # Define explicit time encoding
        encoding = {}
        if 'time' in ds:
            encoding['time'] = {
                'units': 'seconds since 1970-01-01 00:00:00',
                'calendar': 'gregorian'
            }

        # Write to NetCDF file
        ds.to_netcdf(nc_filename, encoding=encoding)
        converted_files.append(nc_filename)
        print(f"  [{filepath.name}] Saved CF-compliant NetCDF ({len(df)} records) -> {nc_filename.name}")

    print(f"\nDone! Successfully checked or converted {len(converted_files)} files.")
    return converted_files

@dask.delayed
def interpolate_single_forecast_cycle(cycle_dir, init_time, ds_sofar_ref):
    """Worker: Processes exactly one forecast cycle and returns the interpolated Dataset."""
    gfs_hour = f"{init_time.hour:02d}"
    f_files = sorted(glob.glob(f"{cycle_dir}/gfs.t{gfs_hour}z.global.0p25.f*.grib2"))
    
    if not f_files:
        return None

    ds_gfs = xr.open_mfdataset(
        f_files,
        engine="cfgrib",
        combine="nested",
        concat_dim="step",
        coords="minimal",
        compat="override",
        chunks={},  
        parallel=False, 
        backend_kwargs={"filter_by_keys": {"typeOfLevel": "surface"}, "indexpath": ""}
    )

    if "time" in ds_gfs.dims:
        ds_gfs = ds_gfs.squeeze("time")
        
    ds_gfs = ds_gfs.assign_coords(
        valid_time=("step", (np.datetime64(init_time) + ds_gfs["step"]).data)
    )
    ds_gfs = ds_gfs.swap_dims({"step": "valid_time"}).sortby("valid_time")

    gfs_vars = ["swh", "ws", "wdir", "mp1", "wvdir"]
    ds_gfs = ds_gfs[gfs_vars]
    
    start_time = ds_gfs.valid_time.min().values
    end_time = ds_gfs.valid_time.max().values
    ds_sofar_sub = ds_sofar_ref.sel(time=slice(start_time, end_time))

    if len(ds_sofar_sub.time) == 0:
        return None  

    ds_gfs = ds_gfs.chunk({"valid_time": -1})
    sofar_dim = list(ds_sofar_sub.sizes.keys())[0]
    ds_sofar_sub = ds_sofar_sub.chunk({sofar_dim: -1})

    ds_model_interp = ds_gfs.interp(
        valid_time=ds_sofar_sub.time,
        latitude=ds_sofar_sub.latitude,
        longitude=ds_sofar_sub.longitude,
        method="linear"
    )

    fcst_hr_vals = (ds_sofar_sub.time.values - np.datetime64(init_time)) / np.timedelta64(1, 'h')
    
    ds_combined = xr.Dataset(
        data_vars={
            "latitude": ("time", ds_sofar_sub["latitude"].data),
            "longitude": ("time", (ds_sofar_sub["longitude"] % 360).data), 
            "fcst_hr": ("time", fcst_hr_vals),
            "model_hs": ("time", ds_model_interp["swh"].data),
            "obs_hs": ("time", ds_sofar_sub["sea_surface_wave_significant_height"].data),
            "model_wnd": ("time", ds_model_interp["ws"].data),
            "obs_wnd": ("time", ds_sofar_sub["wind_speed"].data),
            "model_wdir": ("time", ds_model_interp["wdir"].data),
            "obs_wdir": ("time", ds_sofar_sub["wind_from_direction"].data),
            "model_mwp": ("time", ds_model_interp["mp1"].data),
            "obs_mwp": ("time", ds_sofar_sub["sea_surface_wave_mean_period"].data),
            "model_mwdir": ("time", ds_model_interp["wvdir"].data),
            "obs_mwdir": ("time", ds_sofar_sub["sea_surface_wave_from_direction_at_variance_spectral_peak"].data),
        },
        coords={"time": ds_sofar_sub.time.data},
        attrs={
            "satellite_name": "SOFAR_BUOY",
            "model_name": "GFS_WAVE",
            "initialization_time": str(init_time)
        }
    )
    
    return ds_combined.compute()

@dask.delayed
def process_single_cycle(cycle_dir, init_time, sofar_path, model_name):
    """Worker: Processes exactly one forecast cycle and returns the interpolated Dataset."""
    print(f"[{init_time}] Starting processing...")
    
    # --- FIX: Load the reference data natively inside the worker ---
    # Using open_zarr allows the worker to lazily read only the time chunk it needs
    ds_sofar_ref = xr.open_zarr(sofar_path)
    # ---------------------------------------------------------------
    
    gfs_hour = f"{init_time.hour:02d}"
    f_files = sorted(glob.glob(f"{cycle_dir}/gfs.t{gfs_hour}z.global.0p25.f*.grib2"))
    
    if not f_files:
        print(f"[{init_time}] SKIP: No GRIB files found.")
        return None

    # ... (Keep all your existing GFS loading and cleanup code exactly the same) ...
    ds_gfs = xr.open_mfdataset(
        f_files, engine="cfgrib", combine="nested", concat_dim="step",
        coords="minimal", compat="override", chunks={}, parallel=False, 
        backend_kwargs={"filter_by_keys": {"typeOfLevel": "surface"}, "indexpath": ""}
    )

    if "time" in ds_gfs.dims:
        ds_gfs = ds_gfs.squeeze("time", drop=True)
        
    coords_to_keep = ["latitude", "longitude", "step"]
    coords_to_drop = [c for c in ds_gfs.coords if c not in coords_to_keep]
    if coords_to_drop:
        ds_gfs = ds_gfs.drop_vars(coords_to_drop)
        
    ds_gfs = ds_gfs.assign_coords(valid_time=("step", (np.datetime64(init_time) + ds_gfs["step"]).data))
    ds_gfs = ds_gfs.swap_dims({"step": "valid_time"}).sortby("valid_time")
    ds_gfs = ds_gfs[["swh", "ws", "wdir", "mp1", "wvdir"]]
    
    start_time = ds_gfs.valid_time.min().values
    end_time = ds_gfs.valid_time.max().values
    
    # --- Subset and compute ONLY the small slice needed for this worker ---
    mask = (ds_sofar_ref.time >= start_time) & (ds_sofar_ref.time <= end_time)
    ds_sofar_sub = ds_sofar_ref.isel(obs=mask.compute()).compute()
    # ----------------------------------------------------------------------

    if ds_sofar_sub.sizes["obs"] == 0:
        print(f"[{init_time}] SKIP: No SOFAR observations found for this time window.")
        return None

    print(f"[{init_time}] Interpolating {ds_sofar_sub.sizes['obs']} SOFAR observations...")
    
    ds_gfs = ds_gfs.chunk({"valid_time": -1})
    ds_sofar_sub = ds_sofar_sub.chunk({"obs": -1})

    ds_model_interp = ds_gfs.interp(
        valid_time=ds_sofar_sub.time,
        latitude=ds_sofar_sub.latitude,
        longitude=ds_sofar_sub.longitude,
        method="linear"
    )

    fcst_hr_vals = (ds_sofar_sub.time.values - np.datetime64(init_time)) / np.timedelta64(1, 'h')
    init_time_array = np.full(ds_sofar_sub.sizes["obs"], np.datetime64(init_time))
    
    ds_combined = xr.Dataset(
        data_vars={
            "init_time": ("obs", init_time_array),
            "time": ("obs", ds_sofar_sub.time.data),
            "latitude": ("obs", ds_sofar_sub["latitude"].data),
            "longitude": ("obs", (ds_sofar_sub["longitude"] % 360).data), 
            "fcst_hr": ("obs", fcst_hr_vals),
            
            "model_hs": ("obs", ds_model_interp["swh"].data),
            "obs_hs": ("obs", ds_sofar_sub["sea_surface_wave_significant_height"].data),
            "model_wnd": ("obs", ds_model_interp["ws"].data),
            "obs_wnd": ("obs", ds_sofar_sub["wind_speed"].data),
            "model_wdir": ("obs", ds_model_interp["wdir"].data),
            "obs_wdir": ("obs", ds_sofar_sub["wind_from_direction"].data),
            "model_mwp": ("obs", ds_model_interp["mp1"].data),
            "obs_mwp": ("obs", ds_sofar_sub["sea_surface_wave_mean_period"].data),
            "model_mwdir": ("obs", ds_model_interp["wvdir"].data),
            "obs_mwdir": ("obs", ds_sofar_sub["sea_surface_wave_from_direction_at_variance_spectral_peak"].data),
        },
        attrs={
            "platform_id": "SOFAR Spotter Buoy",
            "model_name": f"{model_name}"
        }
    )
    
    result = ds_combined.compute()
    print(f"[{init_time}] SUCCESS: Completed interpolation.")
    return result


def process_matchup_pipeline(start_date, end_date, sofar_dir, gfs_dir, model_name):
    """Master Pipeline: Loads SOFAR data, queues Dask workers for GFS cycles, and concatenates the result."""
    from dask.distributed import get_client, progress
    
    print(f"--- Starting Pipeline ---")
    print(f"Target Date Range: {start_date} to {end_date}")
    
    print("Loading SOFAR Buoy data...")
    sofar_files = sorted(glob.glob(f"{sofar_dir}/sofar_*.nc"))
    print(f"Found {len(sofar_files)} SOFAR files.")
    
    ds_list = []
    for f in sofar_files:
        ds = xr.open_dataset(f, engine="netcdf4", chunks={})
        
        # Identify the native physical dimension and unify it as 'obs'
        native_dim = [d for d in ds.sizes if d != "time"][0]
        if native_dim != "obs":
            ds = ds.rename_dims({native_dim: "obs"})
            
        # Ensure time, lat, and lon are explicitly set as coordinates
        ds = ds.set_coords(["time", "latitude", "longitude"])
        ds_list.append(ds)
        
    print("Concatenating SOFAR observations (preserving all duplicates)...")
    ds_sofar = xr.concat(ds_list, dim="obs")
    ds_sofar = ds_sofar.sortby("time")
    ds_sofar["longitude"] = (ds_sofar["longitude"] + 360) % 360
    
    # --- FIX: Write reference to fast HPC storage ---
    temp_sofar_path = "./temp_sofar_reference.zarr"
    print(f"Writing SOFAR reference to disk at {temp_sofar_path}...")
    
    # Compute and save as a consolidated Zarr store for lightning-fast parallel reads
    ds_sofar.compute().to_zarr(temp_sofar_path, mode="w", consolidated=True)
    # ------------------------------------------------
    
    print("Scanning GFS forecast cycle directories...")
    cycle_dirs = sorted(glob.glob(f"{gfs_dir}/gfs.*/[0-2][0-9]/products/wave/gridded/global.0p25"))
    print(f"Found {len(cycle_dirs)} total cycles on disk.")
    
    tasks = []
    for cycle_dir in cycle_dirs:
        parts = Path(cycle_dir).parts
        gfs_date = [p for p in parts if p.startswith("gfs.")][0].split(".")[1]
        gfs_hour = parts[parts.index(f"gfs.{gfs_date}") + 1]
        
        init_time = pd.to_datetime(f"{gfs_date} {gfs_hour}:00:00")
        
        if pd.to_datetime(start_date) <= init_time <= pd.to_datetime(end_date):
            # Pass the tiny string path instead of the data!
            tasks.append(process_single_cycle(cycle_dir, init_time, temp_sofar_path, model_name))

    print(f"Filtered to {len(tasks)} forecast cycles within the date range.")
    
    if len(tasks) == 0:
        print("No tasks queued. Exiting.")
        return None

    print(f"Computing in parallel across Dask workers...")
    
    try:
        client = get_client()
        futures = client.compute(tasks)
        progress(futures)
        results = client.gather(futures)
    except ValueError:
        from dask.diagnostics import ProgressBar
        with ProgressBar():
            results = dask.compute(*tasks)

    valid_results = [ds for ds in results if ds is not None]
    
    if not valid_results:
        print("No valid matchups generated (all returned None).")
        return None

    print(f"Successfully processed {len(valid_results)} cycles. Concatenating along 'obs' dimension...")
    ds_final = xr.concat(valid_results, dim="obs")
    
    print(f"Final dataset generated with {len(ds_final.obs)} total matched observations.")
    return ds_final