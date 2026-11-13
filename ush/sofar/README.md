# SOFAR Spotter Buoy-Model Matchup Pipeline
Interpolates wave model forecasts (e.g., NOAA GFS) with drifting SOFAR Spotter buoy observations.

### Scripts & Inputs 
The workflow is controlled by two bash scripts. Both scripts can be run from anywhere—just update the absolute paths in the user inputs.

1. **`get_sofar_data.sh`**: Fetches SOFAR data from HPSS system on Ursa.
2. **`make_sofar_matchup.sh`**: Runs the parallelized 3D interpolation and generates the final matchup files. Output is converted netCDF files for all SOFAR bufr files and both a netCDF and Zarr matchup file. 

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

---
### Pipeline Architecture
1. **Data Check:** Verifies/downloads required SOFAR data from HPSS.
2. **Conversion:** Converts BUFR files to NetCDF.
3. **Zarr Staging:** Consolidates sofar observations for fast parallel reading.
4. **Interpolation:** Maps forecast cycle data to buoy time and positions.
5. **Consolidation:** Aggregates final data into unified netCDF and Zarr store.

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