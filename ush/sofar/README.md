# SOFAR Spotter Buoy Model Evaluation

## Spotter Buoy-Model Matchup Pipeline 
The SOFAR Spotter buoy observation-model matchup pipeline provides a workflow to collocate operational or experimental wave model forecasts (such as NOAA GFS wave runs) with drifting SOFAR Spotter buoy observations. 

---
### Pipeline Architecture
1. **SOFAR Spotter data check and access:** Check that SOFAR spotter data exists for the full cycle comparison. This includes from the first forecast cycle on the START_DATE user input up to 16 days after the user-defined END_DATE (explained in the inputs section below). If the data does not exist, the pipeline is killed and a new access script is written and executed to access the necessary SOFAR data from the HPSS archive.

2. **BUFR to NetCDF Conversion:** Data is checked to see if all of the SOFAR data exists as netCDF files and if it does not it is automatically converted from BUFR files to netCDFs. This function relies on the [nceplibs-bufr](https://github.com/noaa-emc/nceplibs-bufr) package which is installed into a conda environment via [this](https://anaconda.org/channels/conda-forge/packages/nceplibs-bufr/overview) conda forge host. It is addtionally prebuilt as a module on Ursa that can be loaded. See the [Environment](#setting-up-environment) section below.

3. **Reference SOFAR Zarr Store Staging:** Consolidates observation datasets into a reference store (Zarr) between the given date range for fast parallel reading across cluster nodes.

4. **Parallelized 3D Interpolation:** Maps all forecast cycles (00z, 06z, 12z, 18z) for the given model with initialization times between the user defined start and end dates to parallel Dask workers. Each worker is referencing the SOFAR Zarr store and building a matchup for the corresponding 16 day forecast run. Model fields are linearly interpolated onto the dynamic positions and timestamps of the Spotter buoys.

7. **Matchup Consolidation:** Aggregates interpolated forecast data and observation data into a unified NetCDF and Zarr store for the given start and end dates. The variables are described in the [Matchup Variable Description](#matchup-variable-descriptions) section below.

---
### Pipeline Inputs 
All user inputs are declared in the bash job script `submit_sofar.sh.` This script has user editable inputs to control the data directories, parallelization parameters, and output files.

* **START_DATE="2026-01-01"**

  The start date for the model-observation matchup. 

* **END_DATE="2026-01-31"**
  
  The end date for the model-observation matchup. **NOTE: The start and end dates are relative to the forecast initialization times (i.e. matchup with START_DATE="2026-01-01" and END_DATE="2026-02-02" will have 4 runs per day (00Z, 06Z, 12Z, 18Z) for two days, the matchups will include all forecasts intialized between the given start and end dates.
  

* **WORK_DIR="/scratch3/NCEPDEV/marine/Edwin.Rainville/sofar_work_dir"**

  Directory where the created netCDF files and the matchups will be saved at the end of the pipeline. This directory can be anywhere and all output files will be saved here.

* **OBS_DIR="/scratch3/NCEPDEV/marine/Edwin.Rainville/data/sofar_bufr_data"**
  
  Directory where the observation bufr data is stored.
  
* **MODEL_COMPARE_DIR="/scratch3/NCEPDEV/climate/Jessica.Meixner/Data/retrov17_01/"**

  Directory where the model data to compare with is stored. Default path is for Ursa. Files in this directory are read-only.
  
* **MODEL_PATTERN="gfs."**

  Pattern for the forecast directories that are appended with the date. 
  
* **MODEL_NAME="GFS_WAVE"**

  Model name that is saved in the metadata for the final matchup file.
  
* **MODEL_FILE_TEMPLATE="gfs.{dt:%Y%m%d}/{dt:%H}/products/wave/gridded/global.0p25/gfs.t{dt:%H}z.global.0p25.f\*.grib2"**

  Forecast file path template that is adjusted with the dates and forecast cycles. 
---
### How to run the pipeline
1. Clone repo to Ursa.
2. Edit user input section of `submit_sofar.sh` to ensure that necessary data is available in the observation and model directories.
3. Submit the job script with the command `sbatch submit_sofar.sh`
4. View the log file in real-time with the command `tail -f submit_sofar.out` or view the log at the end and completed matchup data in the user defined work directory.

---
### Matchup Variable Descriptions

| Variable | Description | Units |
| :--- | :--- | :--- |
| `time` | Timestamp of observation | datetime64[ns] |
| `buoy_id` | String of buoy name | str |
| `latitude` | Buoy latitude coordinate | degrees_north |
| `longitude` | Buoy longitude coordinate ($0$ to $360$) | degrees_east |
| `init_time` | Forecast cycle initialization timestamp | datetime64[ns] |
| `fcst_hr` | Forecast lead time since initialization | hours |
| `obs_hs` | Significant wave height from SOFAR Spotter | meters |
| `model_hs` | Significant wave height interpolated from model | meters |
| `obs_wnd` | Wind speed from SOFAR Spotter | m s-1 |
| `model_wnd` | Wind speed interpolated from model | m s-1 |
| `obs_wdir` | Wind direction from SOFAR Spotter | degree |
| `model_wdir` | Wind direction interpolated from model | degree |
| `obs_mwp` | Mean wave period from SOFAR Spotter | seconds |
| `model_mwp` | Mean wave period interpolated from model | seconds |
<!-- | `obs_mwdir` | Dominant peak wave direction from SOFAR Spotter | degree |
| `model_mwdir` | Mean wave direction interpolated from model | degree | -->

The relevant BUFR table where the units of the observation data are accessed is [here](https://www.nco.ncep.noaa.gov/sib/jeff/bufrtab_tableb.html). 

### Setting up environment
**Option 1: Prebuilt module environment**

Ursa includes prebuilt module environments that can be loaded to build the necessary environment. The required modules for this pipeline are the following: 

`module use /scratch4/NCEPDEV/marine/Saeideh.Banihashemi/installs/python-modules/`

`module load stack-oneapi/2024.2.1`

`module load ursa-env`

`module load bufr`

Additionally, the `ursa-env` needs to be prioritized over the `bufr` so the following settings are necessary: 

`export PATH="/scratch4/NCEPDEV/marine/Saeideh.Banihashemi/installs/ursa-env/bin:$PATH"`

The ww3tools functions also need to be added as importable so the ww3tools need to be added to the `PYTHONPATH` environment variable. This is accomplised in the submit script with the following commands:
`cd $SLURM_SUBMIT_DIR`
`export PYTHONPATH="../../:${PYTHONPATH}"`


**Option 2: Conda**

Currently designed to use a conda environment but can expanded to venv in the future. This setup follows the guidelines of NOAA RDHPCS for conda usage stated [here](https://docs.rdhpcs.noaa.gov/software/python/conda_basics.html).

Creating the SOFAR environment on Ursa: 
1. Start an interactive service node: `salloc -A marine-cpu --partition=u1-service -N 1 --mem=8G -t 01:00:00`
2. Load the conda module: `module load rdhpcs-conda`
3. Create the environment: `conda env create --file environment.yml --prefix /scratch[34]/<LAB>/<PROJECT_ID>/<USER_ID>/envs/sofar`
4. Activate the environment: `conda activate /scratch[34]/<LAB>/<PROJECT_ID>/<USER_ID>/envs/sofar`

**Note: to avoid needing to write the path for all future use of this environment make sure it is installed in a directory called `envs` in the location stated above. Then use the command `conda config --append envs_dirs <project_home>/<PROJECT_ID>/<USER_ID>/envs/` to add this directory to be easily discovered by conda. 

<!-- ### Available Functions (`ww3tools.sofar_tools`)

#### `check_model_data`
```python
check_model_data(model_dir, start_date_str, end_date_str, pattern_prefix="gfs.") -->
