import glob
from pathlib import Path
import pandas as pd
import numpy as np
import xarray as xr
import dask

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