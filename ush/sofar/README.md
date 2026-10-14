# SOFAR Spotter Buoy Model Evaluation

## Spotter Buoy-Model Matchup Pipeline 
The SOFAR Spotter buoy observation-model matchup pipeline provides a workflow to collocate operational or experimental wave model forecasts 
(such as NOAA GFS wave runs) with drifting SOFAR Spotter buoy observations. 


---

### Pipeline Architecture

* **Model Data Pre-Check:** Verifies that required forecast run directories and files exist for every day within the specified date range before cluster execution.
* **BUFR to NetCDF Extraction:** Translates NCEP marine BUFR observation files into CF-compliant NetCDF trajectory files, applying an automated forward date buffer to encompass model forecast horizons. Existing NetCDF conversions are detected and skipped automatically.
* **Reference Store Staging:** Consolidates observation datasets into a reference store (Zarr or multi-file NetCDF) for parallel reading across cluster nodes.
* **Distributed 4D Collocation:** Maps forecast initialization dates and six-hourly cycle times (00z, 06z, 12z, 18z) across Dask worker processes. Model fields are linearly interpolated onto the dynamic positions and timestamps of the Spotter buoys.
* **Matchup Consolidation:** Aggregates interpolated forecast records into a unified NetCDF file indexed along the observation dimension, tagged with forecast lead hours and run metadata.

---

### Matchup Variable Descriptions

| Variable | Description | Units |
| :--- | :--- | :--- |
| `time` | Timestamp of observation | datetime64[ns] |
| `latitude` | Buoy latitude coordinate | degrees_north |
| `longitude` | Buoy longitude coordinate ($0$ to $360$) | degrees_east |
| `init_time` | Forecast cycle initialization timestamp | datetime64[ns] |
| `fcst_hr` | Forecast lead time since initialization | hours |
| `obs_hs` | Significant wave height from SOFAR Spotter | meters |
| `model_hs` | Significant wave height interpolated from model | meters |
| `obs_wnd` | Surface wind speed from SOFAR Spotter | m s-1 |
| `model_wnd` | 10-meter wind speed interpolated from model[cite: 1, 2] | m s-1 |
| `obs_wdir` | Surface wind direction from SOFAR Spotter | degree |
| `model_wdir` | Wind direction interpolated from model | degree |
| `obs_mwp` | Mean wave period from SOFAR Spotter | seconds |
| `model_mwp` | Mean wave period interpolated from model | seconds |
| `obs_mwdir` | Dominant peak wave direction from SOFAR Spotter | degree |
| `model_mwdir` | Mean wave direction interpolated from model | degree |

### Available Functions (`ww3tools.sofar_tools`)

#### `check_model_data`
```python
check_model_data(model_dir, start_date_str, end_date_str, pattern_prefix="gfs.")
