# RETROTESTS EVALUATION TOOLS

A collection of Python scripts for wave evaluation and visualization.  
These tools generate QQ plots, Taylor diagrams, colored scatter plots, and global maps of Bias, RMSE, and sample counts for each forecast day.  
In addition, they provide statistical analysis and plots of key metrics as a function of forecast lead time.

## Script Descriptions
### - Configuration setting
- `evalsumconfig.json`: Defines evaluation configuration, including model directories, model names, filenames for combined datasets (both full-period and per-forecast-day), satellite name, season, and output directory.

### - Evaluation and visualization for each forecast day
- `eval.py`: Core script for loading and processing per-forecast-day datasets. Generates evaluation figures including QQ plots, Taylor diagrams, colored scatter plots, and global maps of Bias, RMSE, and sample counts. Plotting is performed using classes from `pvalstats.py`.
- `job_eval.sh`: Slurm job script to run `eval.py`.

### - Statistical analysis and plots
- `eval_stat.py`: Core script for loading the all-in-one combined NetCDF dataset, performing statistical analysis, and generating plots of Bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, and sample count (N) as functions of forecast lead time. The statistical computation is performed using `mvalstats.py`.
- `job_eval_stat.sh`: Slurm job script to run `eval_stat.py`.

## Usage
### - Evaluation for each forecast day
#### 1. Edit `evalsumconfig.json` to define:
- directories     : input directories for each model
- filenames       : list of combined file prefixes to process
- satellite_name  : satellite name used in the filenames and plot labels
- season          : season label used in the filenames
- output_dir      : directory for output plots

#### 2. Modify the global map settings in this script as needed:
- DLAT:          latitude bin size
- DLON:          longitude bin size
- MIN_COUNT:     minimum data required per bin, MIN_COUNT = 1 means no filtering
- LATMIN:        minimum latitude included in analysis
- LATMAX:        maximum latitude included
- HS_BIAS_VMAX:  maximum absolute colorbar range for Hs bias
- HS_RMSE_VMAX:  maximum colorbar range for Hs RMSE
- WND_BIAS_VMAX: maximum absolute colorbar range for wind bias
- WND_RMSE_VMAX: maximum colorbar range for wind RMSE
- QC_HS:         quality-control thresholds for Hs
- QC_WND:        quality-control thresholds for wind speed

#### 3. OUTPUT:
The script generates the following plot types for each configured input file:
1. Hs QQ plot comparing all models
2. Hs Taylor diagram comparing all models
3. WND QQ plot comparing all models
4. WND Taylor diagram comparing all models
5. Hs scatter plot for each model
6. WND scatter plot for each model
7. Hs global bias map for each model
8. Hs global RMSE map for each model
9. Hs global sample counts map for each model (PR104)
10. WND global bias map for each model
11. WND global RMSE map for each model
12. WND global sample counts map for each model (PR104)

#### 4. Output filename formats:
- Hs QQ / Taylor plots: plot_HS_{filename}_{satellite_name}_{season}*
- WND QQ / Taylor plots: plot_WND_{filename}_{satellite_name}_{season}*
- Hs scatter plots: plot_HS_scatter_{filename}_{satellite_name}_{season}_{model_label}_*
- WND scatter plots: plot_WND_scatter_{filename}_{satellite_name}_{season}_{model_label}_*
- Hs global bias map: plot_Hs_{model_label}_{filename}_{satellite_name}_global_Bias.png
- Hs global RMSE map: plot_Hs_{model_label}_{filename}_{satellite_name}_global_RMSE.png
- WND global bias map: plot_WND_{model_label}_{filename}_{satellite_name}_global_Bias.png
- WND global RMSE map: plot_WND_{model_label}_{filename}_{satellite_name}_global_RMSE.png

#### NOTE
- The script currently assumes two models when assigning suffixes and labels: retrov17_01 and gfsv16
- Hs uses `obs_hs` as the observation field.
- WND uses `obs_wnd_cal` as the observation field.
- Only files that exist for all configured models are processed.
- Rows with missing values are removed before plotting.

### - Statistical analysis and plots for all forecast days
#### 1. Edit `evalsumconfig.json` to define:
- directories    : input directories for each model
- filenames_all  : all-in-one combined file prefix
- satellite_name : satellite name used in the filenames and plot titles
- season         : season label used in the filenames and plot titles
- output_dir     : directory for output figures

#### 2. Modify the script as needed for:
- endday        : maximum forecast day to include
- n_filters     : number of filters to apply

#### 3. OUTPUT:
One PNG figure is created for each verification metric. Each figure contains two panels:
- top panel    : Hs metric versus forecast hour
- bottom panel : WND metric versus forecast hour

Output filename format: fig_{stat_name}_{satellite_name}_{season}.png

#### NOTE:
- Statistics: The script computes the following metrics from `mvalstats.metrics`: bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N

- Forecast bins: Statistics are computed in cumulative 24-hour forecast bins:
```
            Day 1  :   0 < fcst_hr <=  24
            Day 2  :  24 < fcst_hr <=  48
            ...
            Day 16 : 360 < fcst_hr <= 384
```

- Filters: The script supports up to three filter levels:
   - filter 0 : all samples
   - filter 1 : obs_hs >= 4 m
   - filter 2 : obs_hs >= 7 m

- Others:
   - Hs statistics are computed from: model_hs vs obs_hs
   - WND statistics are plotted from: model_wnd vs obs_wnd_cal
