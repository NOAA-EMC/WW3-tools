import os
from datetime import datetime, timedelta, time
from typing import List
import xarray as xr

def get_corresponding_files(
    base_dir: str,
    file_pattern: str,
    start_date: str,
    end_date: str,
    loop_type: str = "date",
) -> List[str]:
    """Finds existing files in a directory across an inclusive date range.

    Truncates datetime inputs to day precision (YYYY-MM-DD) and includes
    all corresponding days or forecast hours through the end of `end_date`.

    Args:
        base_dir (str): Target directory containing the files.
        file_pattern (str): Filename structure. Uses strftime for 'date' loops,
            and allows {fhr} formatting tags for 'forecast_hour' loops.
        start_date (str): Start date string (e.g., 'YYYY-MM-DD').
        end_date (str): End date string (e.g., 'YYYY-MM-DD').
        loop_type (str): Either 'date' (checks each day) or 'forecast_hour'
            (checks every hour from start_date 00:00 through end_date 23:00).
            Default is 'date'.

    Returns:
        List[str]: A list of file paths that exist on disk within the date range.

    Raises:
        ValueError: If `loop_type` is not 'date' or 'forecast_hour', or if date
            formats cannot be parsed.
    """
    if loop_type not in ["date", "forecast_hour"]:
        raise ValueError("loop_type must be either 'date' or 'forecast_hour'.")

    # Helper to parse string inputs down to the date component (YYYY-MM-DD)
    def parse_to_date(date_str: str):
        for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d", "%Y%m%d%H", "%Y%m%d"):
            try:
                return datetime.strptime(date_str, fmt).date()
            except ValueError:
                pass
        raise ValueError(f"Time format for '{date_str}' not recognized.")

    start_d = parse_to_date(start_date)
    end_d = parse_to_date(end_date)

    if start_d > end_d:
        raise ValueError(f"start_date ({start_d}) cannot be after end_date ({end_d}).")

    start_dt = datetime.combine(start_d, time.min)
    valid_paths = []

    if loop_type == "forecast_hour":
        # End date is inclusive through 23:00:00 of end_d
        end_dt = datetime.combine(end_d, time(23, 59, 59))
        current_dt = start_dt

        while current_dt <= end_dt:
            fhr = int((current_dt - start_dt).total_seconds() // 3600)
            file_name = current_dt.strftime(file_pattern).format(fhr=fhr)
            full_path = os.path.join(base_dir, file_name)

            if os.path.exists(full_path):
                valid_paths.append(full_path)

            current_dt += timedelta(hours=1)
    else:
        # Daily check inclusive of end_d
        current_d = start_d
        while current_d <= end_d:
            current_dt = datetime.combine(current_d, time.min)
            file_name = current_dt.strftime(file_pattern)
            full_path = os.path.join(base_dir, file_name)

            if os.path.exists(full_path):
                valid_paths.append(full_path)

            current_d += timedelta(days=1)

    return valid_paths

def combine_netcdf_files(file_paths: List[str], concat_dim: str = "time") -> xr.Dataset:
  """Opens a list of NetCDF file paths individually with xarray and concatenates
  them into a single Dataset.

  Args:
      file_paths (List[str]): List of NetCDF file paths to process.
      concat_dim (str): Dimension along which to concatenate (default: 'time').

  Returns:
      xr.Dataset: The combined xarray Dataset.
  """
  if not file_paths:
    raise ValueError("The provided file_paths list is empty.")

  # 1. Open each file in the list individually
  datasets = [xr.open_dataset(path) for path in file_paths]

  # 2. Concatenate all opened datasets along the specified dimension
  combined_ds = xr.concat(datasets, dim=concat_dim)

  # 3. Sort along the time dimension to ensure chronological order
  if concat_dim in combined_ds.coords or concat_dim in combined_ds.dims:
    combined_ds = combined_ds.sortby(concat_dim)

  return combined_ds

def clean_sofar_data(ds: xr.Dataset, threshold: float = 500.0) -> xr.Dataset:
    """Cleans Sofar buoy dataset by removing unphysical fill values (e.g., values

    >= 500).

    Args:
        ds (xr.Dataset): Input Sofar buoy dataset.
        threshold (float): Upper boundary threshold for physical values
          (default: 500.0).

    Returns:
        xr.Dataset: Filtered xarray Dataset with unphysical records dropped.
    """
    # Key variables to apply quality control filtering on
    qc_vars = ["hs", "wspd", "wdir", "avwp", "spwp", "ddwc", "dsdw"]

    # Build a combined boolean mask across all present variables
    mask = True
    for var in qc_vars:
      if var in ds.data_vars:
        mask = mask & (ds[var] < threshold)

    # Filter dataset and drop points where the condition fails
    return ds.where(mask, drop=True)

def combine_grib_files(
    file_paths: List[str], 
    engine: str = "cfgrib",
    chunks: Optional[Dict[str, int]] = None
) -> xr.Dataset:
    """Opens a list of GRIB files lazily using Dask, concatenates them along 
    the forecast step, and computes a valid time coordinate.

    Args:
        file_paths (List[str]): List of GRIB file paths to process.
        engine (str): Engine for xarray to use (default: 'cfgrib').
        chunks (dict, optional): Chunk sizes for Dask lazy loading. 
            Example: {"step": -1, "latitude": 100, "longitude": 100}.
            If None, defaults to Dask's auto-chunking ({}).

    Returns:
        xr.Dataset: The combined xarray Dataset with a valid 'time' dimension.
    """
    if not file_paths:
        raise ValueError("The provided file_paths list is empty.")
        
    # If chunks is None, use an empty dict to trigger Dask's default auto-chunking
    if chunks is None:
        chunks = {} 

    # 1. Open each file individually with Dask chunking enabled
    datasets = [
        xr.open_dataset(path, engine=engine, chunks=chunks) 
        for path in file_paths
    ]

    # 2. Concatenate all opened datasets along the forecast 'step' dimension
    with xr.set_options(use_new_combine_kwarg_defaults=True):
        combined_ds = xr.concat(datasets, dim="step")

    # 3. Shift time variable to be the time reference (start of forecast)
    if "time" in combined_ds:
        combined_ds = combined_ds.rename({"time": "time_ref"})

    # 4. Compute a new time variable as the reference time plus the time step
    if "time_ref" in combined_ds and "step" in combined_ds:
        combined_ds["time"] = combined_ds["time_ref"] + combined_ds["step"]
        combined_ds = combined_ds.swap_dims({"step": "time"})

    # 5. Sort chronologically by the new valid time dimension
    if "time" in combined_ds.coords or "time" in combined_ds.dims:
        combined_ds = combined_ds.sortby("time")

    return combined_ds