# SOFAR Spotter Buoy-Model Validation and Verification Pipeline
Interpolates wave model forecasts (e.g., NOAA GFS) with drifting SOFAR Spotter buoy observations.

### Scripts & Inputs 
Run scripts are designed to be run from anywhere — just update the absolute paths in the user inputs.

1. **`get_sofar_data.sh`**: Fetches SOFAR data from HPSS system on Ursa.
2. **`make_sofar_matchup.sh`**: Runs the parallelized 3D interpolation (lat, lon, time) and generates the final matchup files. Output is converted netCDF files for all SOFAR bufr files and both a netCDF and Zarr matchup file. 
3. **`plot_sofar_validation.sh`**: Computes forecast-hour windowed statistics from one or more matchup files and saves comparison plots. Settings are read from `plot_sofar_validation_config.yaml`.

**Required Input Adjustments:**
Edit the variable block at the top of *both* scripts before running:
* `START_DATE` / `END_DATE`: Date range for forecast initialization times. Format is ('YYYY-MM-DD') for both.
* `WORK_DIR`: Output directory for final NetCDF/Zarr files. Subdirectories will be made for netCDF and zarr data. 
* `OBS_DIR`: Directory for SOFAR BUFR data.
* `MODEL_COMPARE_DIR`: Base directory for model data.
* `MODEL_PATTERN`: Prefix for the daily model run directories (e.g., `"gfs."` for folders like `gfs.20260101`).
* `MODEL_NAME`: Custom label (e.g., `"GFS_WAVE"`) saved as metadata in the final output files to identify the model.
* `MODEL_FILE_TEMPLATE`: File path template to locate forecast files. Uses datetime formatting (e.g., `{dt:%Y%m%d}` and `{dt:%H}`) to dynamically find files for each cycle.

---
### How to Run
1. Clone the repo to Ursa.
2. Adjust the user inputs in `get_sofar_data.sh` and `make_sofar_matchup.sh`.
3. Submit the data retrieval job as a bash script which will spawn parallel slurm jobs to access data from the HPSS system:
   `bash get_sofar_data.sh`
4. Once data retrieval finishes, submit the matchup job:
   `sbatch make_sofar_matchup.sh`
5. Once the matchup job finishes, edit `plot_sofar_validation_config.yaml` and `CONFIG_FILE` / `WW3TOOLS_DIR` in `plot_sofar_validation.sh`, then submit the plotting job:
   `sbatch plot_sofar_validation.sh`

---
### Pipeline Architecture
1. **Data Check:** Verifies/downloads required SOFAR data from HPSS.
2. **Conversion:** Converts BUFR files to NetCDF.
3. **Zarr Staging:** Consolidates sofar observations for fast parallel reading.
4. **Interpolation:** Maps forecast cycle data to buoy time and positions.
5. **Consolidation:** Aggregates final data into unified netCDF and Zarr store.

---
### Validation Plots
`plot_sofar_validation.sh` runs `ww3tools/scripts/sofar_validation_plots.py`, which groups matchups into forecast-hour windows, computes statistics in each window, and plots one line per model. All inputs come from `plot_sofar_validation_config.yaml`:

* `models`: Model label mapped to its matchup file (`.zarr` or `.nc`). Add one entry per model to compare them on the same plots, e.g. `GFS: /path/to/GFS_WAVE_sofar_matchup.zarr`.
* `variables`: Variables to evaluate: `hs` (significant wave height), `wnd` (wind speed), `mwp` (mean wave period).
* `statistics`: Any of `mse`, `rmse`, `corr`, `error_std`, `mean`, `scatter_index` (RMSE / observation mean, in %), `count` (number of matchups per window).
* `window_hours`: Width of each forecast-hour window (default 1). Windows stop a whole window before the last forecast hour so the final window is never partial.
* `outdir`: Directory for the output PNGs, created if needed.

Only paired matchups (valid model value, valid observation, valid forecast hour) are used, so every statistic in a window is computed from the same points. One PNG is saved per statistic and variable, named `timeseries_{stat}_{var}_{window}hr.png` which will include multiple lines for the multiple models being validated if that is defined by the user.

---
### Environment Setup

**Option 1: Ursa Prebuilt Modules (Default)**
The bash scripts automatically load the required prebuilt modules (`stack-oneapi`, `ursa-env`, `bufr`) and set the `PYTHONPATH`. No manual setup is required but an updated path to the location of ww3tools that you are using is needed.

**Option 2: Custom Conda Environment**
A conda environment is built on Ursa at '/scratch3/NCEPDEV/marine/Edwin.Rainville/envs/sofar' and can be activated using `conda activate path/to/env` after `module load rdhpcs-conda`.
To build the environment manually:
`salloc -A marine-cpu --partition=u1-service -N 1 --mem=8G -t 01:00:00`
`module load rdhpcs-conda`
`conda env create --file environment.yml --prefix /scratch[34]/<LAB>/<PROJECT_ID>/<USER_ID>/envs/sofar`
`conda activate /scratch[34]/<LAB>/<PROJECT_ID>/<USER_ID>/envs/sofar`