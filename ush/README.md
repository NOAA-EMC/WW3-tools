# DATA PROCESSING TOOLS
A collection of Python scripts and shell utilities for processing satellite data (AODN and preBUFR), WAVEWATCH III (WW3) and operational workflow forecast outputs.  
These tools support data conversion, interpolation, and preparation for evaluation workflows.

## Script Descriptions
### - PreBUFR satellite data processing
- `preBUFR_to_NetCDF.sh`: Converts preBUFR JASON3 satellite data to NetCDF format using METplus `pb2nc`.
- `pb2nc_to_ww3tools.py`: Rewrites daily MET pb2nc NetCDF files into ww3tools-compatible NetCDF format.
- `pb2nc_merge_altimeter.py`: Merges daily track-style NetCDF files into monthly NetCDF files.

### - AODN satellite data processing
- `makeprocsatsubmit.py`: Generates Slurm job scripts to process raw AODN satellite data (location-based) into time-based NetCDF files over a user-defined time range.
- `makecombinemonthly.sh`: Merges daily AODN NetCDF files into monthly datasets using CDO.

### - Data processing for evaluaiton
- `makesubmitinterp.py`: Generates Slurm job scripts to interpolate model outputs onto satellite tracks.
- `combineSatInterpOut.py`: Core script that creates combined NetCDF datasets for a specified time period using interpolated model–satellite matchup data.
- `makesubmitcombine.py`: Generates Slurm job scripts to run `combineSatInterpOut.py`. (PR102)

### - Other tools
- `run_all_jobs.sh`: Utility script to submit and run multiple Slurm job scripts.

## Data Processing Steps
All processed data are stored under `WORKDIR`, which contains the `WW3-Tools` repository.

### 1. Create a working directory, _**WORKDIR**_ or other names
```
   mkdir -p WORKDIR
```
### 2. Install WW3-tools in _**WORKDIR**_
```
   cd WORKDIR
   git clone https://github.com/NOAA-EMC/WW3-tools
```
### 3. Process satellite data

#### 3.1. Process AODN satellite data

##### 3.1.1. Generates Slurm job scripts to _**WORKDIR/processsatdata/jobsubs**_
```
    cd WORKDIR/WW3-tools/ush
    module use /scratch3/NCEPDEV/climate/Jessica.Meixner/general/modulefiles
    module load ww3tools
```

modify `makeprocsatsubmit.py` and configuration files in _**WORKDIR/WW3-tools/parm**_ if necessary

```
    python makeprocsatsubmit.py
```
##### 3.1.2. Run jobs and output processed data in _**WORKDIR/processsatdata/out/**_
```
    cd WORKDIR/processsatdata/jobsubs
    ./run_all_jobs.sh
```
##### 3.1.3. Merges daily AODN NetCDF files into monthly datasets and save in _**WORKDIR/processsatdata/combineoutmonthly**_
```
    cd WORKDIR/WW3-tools/ush
```
modify `makecombinemonthly.sh`
```
    sbatch makecombinemonthly.sh
```

#### 3.2. Process preBUFR satellite data

##### 3.2.1. Converts preBUFR JASON3 satellite data to NetCDF format and output in _**WORKDIR/processsatdata/pb2nc_out**_
```
    cd WORKDIR/WW3-tools/ush
```
modify `preBUFR_to_NetCDF.sh`
```
    sbatch preBUFR_to_NetCDF.sh
```    
##### 3.2.2. Rewrites daily MET pb2nc NetCDF files into ww3tools-compatible NetCDF format and save in _**WORKDIR/processsatdata/out/**_
```
    cd WORKDIR/WW3-tools/ush
    module use /scratch3/NCEPDEV/climate/Jessica.Meixner/general/modulefiles
    module load ww3tools
```
modify `pb2nc_to_ww3tools.py`
```
    python pb2nc_to_ww3tools.py
```
##### 3.2.3. Merges daily track-style NetCDF files into monthly NetCDF and save to _**WORKDIR/processsatdata/pb2nc_altimeter_monthly**_
```
    cd WORKDIR/WW3-tools/ush
    module use /scratch3/NCEPDEV/climate/Jessica.Meixner/general/modulefiles
    module load ww3tools
```
modify `pb2nc_merge_altimeter.py`
```
    python pb2nc_merge_altimeter.py
```

### 4. Interpolate model predictions to satellite track
#### 4.1. Generates Slurm job scripts and save to _**WORKDIR/processsatdata/jobinterp**_
```
    cd WORKDIR/WW3-tools/ush
    module use /scratch3/NCEPDEV/climate/Jessica.Meixner/general/modulefiles
    module load ww3tools
```
modify `makesubmitinterp.py`
```
    python makesubmitinterp.py
```
#### 4.2. Run the interpolation jobs and save to _**WORKDIR/processsatdata/outinterp**_
```
    cd WORKDIR/processsatdata/jobinterp/{model}
    ./run_all_jobs.sh
```

### 5. Combine NetCDF datasets for a specified time period using interpolated model–satellite matchup data
#### 5.1. Generates Slurm job scripts and save to _**WORKDIR/processsatdata/jobcombine**_ (PR102)
```
    cd WORKDIR/WW3-tools/ush
    module use /scratch3/NCEPDEV/climate/Jessica.Meixner/general/modulefiles
    module load ww3tools
```
modify `makesubmitcombine.py`
```
    python makesubmitcombine.py
```
#### 5.2. Run the combination jobs and save to _**WORKDIR/processsatdata/outcombine**_
```
    cd WORKDIR/processsatdata/jobcombine/{model}
    ./run_all_jobs.sh
```
