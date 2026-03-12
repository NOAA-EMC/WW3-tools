import netCDF4 as nc
import numpy as np
import pandas as pd
import os
import json
import matplotlib.pyplot as plt  # Import matplotlib for saving figures
from pvalstats import ModelObsPlot
from pvalstats import GlobalSkillMap

# Load configuration settings from JSON file
with open('evalsumconfig.json') as config_file:
    config = json.load(config_file)

directories = config["directories"]
filenames = config["filenames"]
satellite_name = config["satellite_name"]
season = config["season"]
output_dir = config["output_dir"]

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# -------------------- global plot setting ------------------------------
DLAT = 1.0             # latitude bin size
DLON = 1.0             # longitude bin size
MIN_COUNT = 10         # minimum data required per bin, MIN_COUNT = 1 means no filtering
LATMIN = -80.0         # minimum latitude included in analysis
LATMAX = 80.0          # maximum latitude included

HS_BIAS_VMAX = 0.8     # color range for SWH, bias range: [-0.5, +0.5] m
HS_RMSE_VMAX = 1.5     # color range, rmse range: [0, 1.0] m
WND_BIAS_VMAX = 2.0    # example m/s (edit to yours)
WND_RMSE_VMAX = 5.0    # example m/s (edit to yours)

# Example for Hs: between 0 and 30 m
QC_HS = {"model_min": 0.0, "model_max": 30.0, "obs_min": 0.0, "obs_max": 30.0}
# Example for WND: between 0 and 60 m/s
QC_WND = {"model_min": 0.0, "model_max": 60.0, "obs_min": 0.0, "obs_max": 60.0}

# -------------------------------------------------------------------------

for filename in filenames:
    # Construct the full file paths for each model
    file_paths = [os.path.join(directories[key], f"{filename}_{key}_{season}_{satellite_name}.nc") for key in directories]

    # Check if all files exist and then process
    if all(os.path.exists(fp) for fp in file_paths):
        datasets = [nc.Dataset(fp, 'r') for fp in file_paths]
        dfs_hs = []
        dfs_wnd = []
        suffixes = ['_retrov17_01', '_gfsv16']
        day = filename.replace("combined_", "")

        for ds, suffix in zip(datasets, suffixes):
            df_hs = pd.DataFrame({
                    "time": ds.variables["time"][:],
                    "latitude": ds.variables["latitude"][:],
                    "longitude": ds.variables["longitude"][:],
                    "obs_hs": ds.variables["obs_hs"][:],
                    "hs" + suffix: ds.variables["model_hs"][:]
            })
            df_hs = df_hs.reset_index().rename(columns={'index': 'row_id'})
            dfs_hs.append(df_hs)

            df_wnd = pd.DataFrame({
                    "time": ds.variables["time"][:],
                    "latitude": ds.variables["latitude"][:],
                    "longitude": ds.variables["longitude"][:],
                    "obs_wnd": ds.variables["obs_wnd_cal"][:],
                    "wnd" + suffix: ds.variables["model_wnd"][:]
            })
            df_wnd = df_wnd.reset_index().rename(columns={'index': 'row_id'})
            dfs_wnd.append(df_wnd)

        # Merging HS dataframes
        merged_hs = dfs_hs[0]
        for df in dfs_hs[1:]:
            merged_hs = pd.merge(merged_hs, df, on=['time', 'row_id', 'latitude', 'longitude', 'obs_hs'], how='inner')

        # Merging WND dataframes
        merged_wnd = dfs_wnd[0]
        for df in dfs_wnd[1:]:
            merged_wnd = pd.merge(merged_wnd, df, on=['time', 'row_id', 'latitude', 'longitude', 'obs_wnd'], how='inner')

        merged_hs = merged_hs.drop(columns=['row_id'])
        merged_wnd = merged_wnd.drop(columns=['row_id'])

        merged_hs = merged_hs.dropna()
        merged_wnd = merged_wnd.dropna()

        # Close datasets
        for ds in datasets:
            ds.close()

        # Create ModelObsPlot for HS
        mop_hs = ModelObsPlot(
            model=np.c_[merged_hs['hs_retrov17_01'], merged_hs['hs_gfsv16']],
            obs=merged_hs['obs_hs'],
            axisnames=["Models", "Satellite"],
            mlabels=["Retrov17_01", "GFSv16"],
            ftag=os.path.join(output_dir, f"plot_HS_{filename}_{satellite_name}_{season}")
        )
        mop_hs.qqplot()
        mop_hs.taylordiagram()

        # Create ModelObsPlot for WND
        mop_wnd = ModelObsPlot(
            model=np.c_[merged_wnd['wnd_retrov17_01'], merged_wnd['wnd_gfsv16']],
            obs=merged_wnd['obs_wnd'],
            axisnames=["Models", "Satellite"],
            mlabels=["Retrov17_01", "GFSv16"],
            ftag=os.path.join(output_dir, f"plot_WND_{filename}_{satellite_name}_{season}")
        )
        mop_wnd.qqplot()
        mop_wnd.taylordiagram()

        model_columns_hs = ['hs_retrov17_01', 'hs_gfsv16']
        model_labels = ['Retrov17_01', 'GFSv16']
        for i, (col, label) in enumerate(zip(model_columns_hs, model_labels)):
             mop_hs_sc = ModelObsPlot(
                model=merged_hs[col].values.reshape(-1, 1),
                obs=merged_hs['obs_hs'].values,
                linreg=True,  # regression line + (your updated) R^2 + equation + corner stats
                axisnames=[f"{label} Hs (m)", f"{satellite_name} Hs (m)"],
                mlabels=[''],
                mtitle=f"Hs (m) {day} {satellite_name} {season}",
                ftag=os.path.join(output_dir, f"plot_HS_scatter_{filename}_{satellite_name}_{season}_{label}_")
             )
             mop_hs_sc.scatterplot(dwscl='yes')

             gsm = GlobalSkillMap(
                lat=dfs_hs[i]['latitude'],
                lon=dfs_hs[i]['longitude'],
                model=dfs_hs[i][col],
                obs=dfs_hs[i]['obs_hs'],
                mlabels=[label]
             )

             gsm.plot_global_like_global_plot(
                metric="bias",
                model_index=0,
                dlat=DLAT,
                dlon=DLON,
                lon_0_360=True,
                min_count=MIN_COUNT,
                latmin=LATMIN,
                latmax=LATMAX,
                qc_kwargs=QC_HS,
                vmax=HS_BIAS_VMAX,
                title=f"Hs Bias – {day} ({label} vs {satellite_name})",
                outfile=os.path.join(output_dir, f"plot_Hs_{label}_{filename}_{satellite_name}_global_Bias.png"),
             )
             gsm.plot_global_like_global_plot(
                metric="rmse",
                model_index=0,
                dlat=DLAT,
                dlon=DLON,
                lon_0_360=True,
                min_count=MIN_COUNT,
                latmin=LATMIN,
                latmax=LATMAX,
                qc_kwargs=QC_HS,
                vmax=HS_RMSE_VMAX,
                title=f"Hs RMSE – {day} ({label} vs {satellite_name})",
                outfile=os.path.join(output_dir, f"plot_Hs_{label}_{filename}_{satellite_name}_global_RMSE.png"),
             )

        model_columns_wnd = ['wnd_retrov17_01', 'wnd_gfsv16']
        for i, (col, label) in enumerate(zip(model_columns_wnd, model_labels)):
             mop_wnd_sc = ModelObsPlot(
                model=merged_wnd[col].values.reshape(-1, 1),
                obs=merged_wnd['obs_wnd'].values,
                linreg=True,
                axisnames=[f"{label} WND (m/s)", f"{satellite_name} WND (m/s)"],
                mlabels=[''],
                mtitle=f"Wind speed (m/s) {day} {satellite_name} {season}",
                ftag=os.path.join(output_dir, f"plot_WND_scatter_{filename}_{satellite_name}_{season}_{label}_")
             )
             mop_wnd_sc.scatterplot(dwscl='yes')

             gsm = GlobalSkillMap(
                lat=dfs_wnd[i]['latitude'],
                lon=dfs_wnd[i]['longitude'],
                model=dfs_wnd[i][col],
                obs=dfs_wnd[i]['obs_wnd'],
                mlabels=[label]
             )

             gsm.plot_global_like_global_plot(
                metric="bias",
                model_index=0,
                dlat=DLAT,
                dlon=DLON,
                lon_0_360=True,
                min_count=MIN_COUNT,
                latmin=LATMIN,
                latmax=LATMAX,
                qc_kwargs=QC_WND,
                vmax=WND_BIAS_VMAX,
                title=f"WND Bias – {day} ({label} vs {satellite_name})",
                outfile=os.path.join(output_dir, f"plot_WND_{label}_{filename}_{satellite_name}_global_Bias.png"),
             )

             gsm.plot_global_like_global_plot(
                metric="rmse",
                model_index=0,
                dlat=DLAT,
                dlon=DLON,
                lon_0_360=True,
                min_count=MIN_COUNT,
                latmin=LATMIN,
                latmax=LATMAX,
                qc_kwargs=QC_WND,
                vmax=WND_RMSE_VMAX,
                title=f"WND RMSE – {day} ({label} vs {satellite_name})",
                outfile=os.path.join(output_dir, f"plot_WND_{label}_{filename}_{satellite_name}_global_RMSE.png"),
             )

    else:
        # Not all file paths exist message
        print(f"Some files for {filename} do not exist.")

