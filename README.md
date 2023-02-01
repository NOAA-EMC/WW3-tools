# WW3-tools
Python tools and utilities for WAVEWATCHIII post-processing and validation.

## Installation
```sh
git clone https://github.com/NOAA-EMC/WW3-tools
cd WW3-tools
pip install .
python3 prep_ww3tools.py
```

## Documentation

### Introduction
&emsp; The third-generation spectral wave model which is actively used for operational wave forecasts, WAVEWATCH III (WW3), is a five-dimensional numerical model covering: space (latitude, longitude), time, and spectrum (direction and frequency). Therefore, the quality of simulations can be evaluated through these dimensions as well as the error signal that can be tracked and sourced based on the model architecture. Additionally, when the model is integrated covering a sufficiently long period of time, and reliable observations are available, a proper validation is possible using error metrics, which can be expanded to analyze the performance in time-domain, probability-domain, and spectral-domain. By considering the dimensional structure of WW3 as well as the characteristics of the outputs through signal processing covering the domains aforementioned, the quality of wave simulations can be assessed. Note that spectral wave models’ extra dimensions require a different type of validation and visualization tools compared to other modeling systems. The existing tools employ conventional methods for statistical analysis, some lack proper skill to quantify spectral wave models. This is an opportunity to deeply investigate the skill and limitations of WW3 simulations using simple visualization, validation, and statistical tools.\
&emsp; In this context, the repository presented here provide a variety of tools to download observations, visualize the wave data (ww3 results and buoys), perform a model validation, and plot statistical results. The descriptions in the header of each program and function as well as the examples provide practical guidance through the aforementioned steps.

### Package structure
&emsp; This repository contains the WW3-tools python package, designed to facilitate the visualization and validation of WW3 results with easy-to-handle modules and functions. The first part is dedicated to the visualization, where codes to plot wave fields and point-outputs (time-series and spectrum) are included. Examples of “wave panels” were also developed, and users can directly visualize in one figure the five dimensions of WW3 results and compare it with buoys – ideal for case studies. Among the wave panels, there are examples to compare two WW3 simulations, side-by-side, accompanied by buoy data – for a straightforward visual comparison.\
&emsp; The second and most important part of the package is the validation, which is ready to read buoy data (NDBC and Copernicus) and satellite altimeter data for WW3 comparisons. The four steps for this task are: (1) pre-processing a grid mask generated with the grid characteristics of WW3 as well as water depth, distance to coast, and information of ocean names and NWS marine zones; (2) cyclone map, generated using IBTracks cyclone information, inserted into the lat/lon grid mask, which allows the option of evaluating the model for cyclonic and non-cyclonic conditions; (3) collocation of WW3 with observations, creating matchups model/buoy or model/satellite data, while also inserting the grid-mask and cyclone information. Step (3) creates a netcdf file with model and observations for the exact time and location, with a suitable data structure that makes easy to sub-sample the later assessments based on water depth (and/or distance to coast), ocean basin, forecast area of interest, and inside/outside the cyclones. The final step (4) contains a function to calculate summary statistics and percentiles of model and observations, another function to calculate error metrics, and several functions for validation plots, such as: QQ-plots, scatter plots, time-series of model and buoy, Taylor diagrams etc.\
&emsp; Scripts to download publicly available observations (altimeters and buoys) and wave forecasts (NOAA GFS and GEFS) are also available so users can, for example, plot and validate the NOAA’s deterministic and ensemble forecasts for any location of interest, using reliable quality-controlled observations.\
&emsp; During the installation the python dependencies required to run the python codes are verified. If you manage to successfully install ww3-tools and run prep_ww3tools.py, then no dependency problems are expected whilst you work with the python scripts.\
&emsp; It is important to mention that all WW3 variables as well as NDBC and Copernicus data hold the metric standard, so wave heights are in meters, periods in seconds etc.

### $\textcolor{darkblue}{ww3tools.downloadobs}$
&emsp; This is a directory with codes to help users to download observations from different sources. Besides ww3 outputs, WW3-tools works with altimeter and buoy data. \
&emsp; The main module wfetchbuoy.py downloads integrated wave parameters and wave spectra from NDBC and Copernicus. \
&emsp; The shell scripts wfetchsatellite_AODN_Altimeter.sh and wfetchsatellite_AODN_Scatterometer.sh download AODN satellite data.

NOAA’s National Data Buoy Center:\
https://www.ndbc.noaa.gov/rsa.shtml \
https://www.ndbc.noaa.gov/measdes.shtml \
https://www.ndbc.noaa.gov/stndesc.shtml \
Copernicus buoy data (require registration and username/password): \
https://marine.copernicus.eu/ \
Integrated Marine Observing System (IMOS) Australian Ocean Data Network (AODN): \
https://portal.aodn.org.au/ \
http://thredds.aodn.org.au/thredds/catalog/IMOS/SRS/Surface-Waves/Wave-Wind-Altimetry-DM00/catalog.html \
https://doi.org/10.1038/s41597-019-0083-9

&emsp; The wind speed data measured by the buoys have multiple different heights. However, when loading and reading buoy data (both NDBC and Copernicus) with wread.py, it converts the wind speeds to 10-meter height. Therefore, wind speed from ww3 outputs, satellite, and buoy data are consistent and fixed to 10 meters.

### $\textcolor{darkblue}{ww3tools}$
&emsp; This directory is the main location where modules and functions are saved. \
&emsp; wread.py is a key python module containing several functions to read WW3 data (field outputs, table of point outputs, and spectrum), NDBC and Copernicus buoy data (table of integrated parameters and spectrum). The header of wread.py shows all the details, inputs and outputs, and examples. More specific information can be found in the header of each function as well.\
&emsp; The WW3 model provides three output types that are plotted by three python scripts, associated with: field outputs, time-series of point-output, and spectral point-output. These three python scripts are respectively ww3fields.py, ww3pointimeseries.py, ww3pointspec.py.\
&emsp; The program ww3fields.py is designed to read both netcdf and grib2 formats, and there are options to customize the area and time steps to be plotted. The code is executed for a specific variable defined by the user as an input argument, so entering the file name and variable is mandatory. By using a structure with arguments, common to all the visualization codes, it is simple to embed it in another script, for example an operational shell script. The next codes, ww3pointimeseries.py and ww3pointspec.py are site-specific, associated with a point-output (single lat/lon). ww3pointimeseries.py plots time-series of integrated parameters for a specific point (input argument) while ww3pointspec.py reads the spectral output for a specific point (input argument) and plots the directional spectra, using a polar plot, and the power spectrum.\
&emsp; The program ndbcpointspec.py is similar to ww3pointspec.py but applied to NDBC spectrum. It calls the function wread.spec_ndbc that reads and mounts the buoy’s directional wave spectrum, and then plot it. All polar spectra use a square-root based scale for the levels, colors, and colorbar.\
&emsp; The results of the aforementioned programs are png figures, saved in the directory where the user is running the code. Since `matplotlib.use('Agg')` has been added to the scripts, no figure is expected to pop-up on the screen.

&emsp; Although the general validation principle looks intuitively simple, there are important scientific questions hidden in this process as well as critical complex steps. Above all, the validation must allow analysts to deeply investigate the model’s performance and its deficiencies, without being restricted to oversimplified results, such as for example associated with the RMSE of Hs alone. The challenge of having a complete validation package while preserving the practicability has been obtained through the structure shown in the next figures. The two figures present the validation steps using altimeter and buoy data, which are very similar. The validation scheme using altimeters has one more step in the process, related to the satellite collocation into a given regular grid, i.e, first the satellite data is collocated and then the matchups ww3/satellite are built. Whereas, for the buoy data, fixed in the same position (moored buoys), the matchups ww3/buoy is straightforward.\
&emsp; The process starts with the download of observations, using the given codes described in the figure, or users can adapt to enter their own measurements. The second step is very important because it builds the foundation of where and how the evaluation will occur, through the lat/lon arrays, water depth, and distance to coast. Including marine zone and offshore zone names are optional. This information is carried through the whole validation process, being very useful at the last steps when users want to determine the errors at coastal zones, or, inversely, to ensure the validation is not affected by the continental shelf or islands, or even to specify the model performance at each ocean and offshore/marine areas. Step (3), associated with the cyclones, is also optional, but it is a great tool that allow users to separate the validation considering cyclonic and non-cyclonic conditions – which used to be a time demanding task that is now simple by using procyclmap.py and IBtracks data.\
&emsp; The satellite collocation with gridSatGlobal_Altimeter.py involves some criteria for the weighted average using pyresample.kd_tree, associated with maximum distance and time to the centered grid point to select from the altimeter tracks. Default values, following the literature, are included; however, depending on the model domain and resolution, users should customize these criteria. The following step, with the model/observation collocation, generates an array with simple structure and dimension. For model/satellite, it is a two-dimensional array, where each column is one feature (model_hs, obs_hs, time, lat, lon, cycloneInfo etc) and the lines are the records of each matchup that has been reshaped. Therefore, if users want to subsample the validation for an specific time interval and location it can be easily done in python by selecting the indexes; for example ind=np.where(cyclone>2) to select model_hs[ind] and obs_hs[ind] for an analysis within the cyclones. The same is valid to select specific oceanic areas, marine/offshore zones, water depths intervals, latitude ranges etc. The main difference between buoy/ww3 matchups and altimeter/ww3 matchups, that must be considered at the final validation steps, is that the buoys/ww3 matchups have three dimensions instead of two, being this extra one associated with buoyID.\
&emsp; Having the arrays of model and observations for the same time and location is the basis for the last part of the process, which includes two scripts with functions for the validation statistics and plots. These functions can be included in the main script, designed by the user, to finally obtain the summary statistics comparisons, error metrics, and plots. There is a WW3-tools directory with examples.

![Screenshot](https://github.com/NOAA-EMC/WW3-tools/blob/develop/docs/BuoyModel_scheme.png)
**<div align="center">Figure 1 - Flowchart of the wave validation process using buoy data.</div>**

![Screenshot](https://github.com/NOAA-EMC/WW3-tools/blob/develop/docs/SatModel_scheme.png)
**<div align="center">Figure 2 - Flowchart of the wave validation process using altimeter data.</div>**

&emsp; The validation package was designed using the methodology described in the following papers, which can be of great support for users to improve and to better interpret the assessments.\
&emsp; Mentaschi L, Besio G, Cassola F, Mazzino, A., 2013. Problems in RMSE- based wave model validations. Ocean Model 72:53–58. https://doi.org/10.1016/j.ocemod.2013.08.003 \
&emsp; Mentaschi, L., Besio, G., Cassola, F., Mazzino, A., 2013. Developing and validating a forecast/hindcast system for the Mediterranean Sea. J. Coastal Res. SI 65, 1551– 1556. http:/dx.doi.org/102112/SI65-262.1. \
&emsp; Jolliff, J.K., Kindle, J.C., Shulman, I., Penta, B., Friedrichs, M.A.M., Helber, R., Arnone, R.A., 2009. Summary diagrams for coupled hydrodynamic-ecosystem model skill assessment. Journal of Marine Systems, 76, 64–82. https://doi.org/10.1016/j.jmarsys.2008.05.014 \
&emsp; Campos, R.M., Alves, J.H.G.M., Penny, S.G., Krasnopolsky, V., 2020. Global assessments of the NCEP Ensemble Forecast System using altimeter data. Ocean Dynamics, 70, 405-419, https://doi.org/10.1007/s10236-019-01329-4 \
&emsp; Willmott C, Matsuura, K., 2005. Advantages of the mean absolute error (MAE) over the root mean square error (RMSE) in assessing average model performance. Clim Res 30(79–82). \
&emsp; Hanna, S., Heinold, D. 1985. Development and application of a simple method for evaluating air quality. In: API Pub. No. 4409, Washington, DC, Washington, USA.

### $\textcolor{darkblue}{ww3tools.opforecast}$
&emsp; Similar to download_observations directory, this place contains programs and functions to download publicly available data, namely operational wave forecast data. Two options are included, for NOAA’s deterministic (GFS) and ensemble (GEFS) wave forecasts: \
https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod \
https://ftpprd.ncep.noaa.gov/data/nccf/com/gens/prod \
&emsp; The operational forecast download and post-processing is restricted to linux environment and require the programs: wget, perl, NCO, and wgrib2. Users can edit the shell scripts to select specific variables of interest and level of netcdf4 compression.

### $\textcolor{darkblue}{examples}$
&emsp; This final directory contains practical examples of visualization codes developed to evaluate and compare ww3 simulations. The wave panels included gather many tools described above and can be used as a suggestion of visualization and model evaluation. The directory is under constant development, including small practical examples of ww3 simulations and validations.

## Final Remarks
&emsp; This is an Open Science initiative. See the webinar \
&emsp; "Small Steps, Big Impact: Supporting Open Science Through Open Access & Data Initiatives": \
&emsp; https://www.youtube.com/watch?v=JC1CCpuhlqw&ab_channel=NOAACentralLibrary
