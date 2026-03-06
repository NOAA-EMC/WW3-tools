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
# -----------------------------------------------------------------------
def _get_var(ds, candidates):
    """Return the first matching variable array from candidates."""
    for v in candidates:
        if v in ds.variables:
            return ds.variables[v][:]
    raise KeyError(f"Missing variable. Tried: {candidates}")

DLAT = 1.0             # latitude bin size
DLON = 1.0             # longitude bin size
MIN_COUNT = 10         # minimum data required per bin, MIN_COUNT = 1 means no filtering
LATMIN = -60.0         # minimum latitude included in analysis
LATMAX = 60.0          # maximum latitude included

HS_BIAS_VMAX = 0.5     # color range for SWH, bias range: [-0.5, +0.5] m
HS_RMSE_VMAX = 1.0     # color range, rmse range: [0, 1.0] m
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
                    "hs" + suffix: ds.variables["model_hs"][:]
            })
            df_hs = df_hs.reset_index().rename(columns={'index': 'row_id'})
            dfs_hs.append(df_hs)

            df_wnd = pd.DataFrame({
                    "time": ds.variables["time"][:],
                    "wnd" + suffix: ds.variables["model_wnd"][:]
            })
            df_wnd = df_wnd.reset_index().rename(columns={'index': 'row_id'})
            dfs_wnd.append(df_wnd)

        # Merging HS dataframes
        merged_hs = dfs_hs[0]
        for df in dfs_hs[1:]:
            merged_hs = pd.merge(merged_hs, df, on=['time', 'row_id'], how='inner')

        # Merging WND dataframes
        merged_wnd = dfs_wnd[0]
        for df in dfs_wnd[1:]:
            merged_wnd = pd.merge(merged_wnd, df, on=['time', 'row_id'], how='inner')


        # Adding observation data

        df_obs_hs = pd.DataFrame({
            'time': datasets[0].variables['time'][:],
            'obs_hs': datasets[0].variables['obs_hs_cal'][:]
        }).reset_index().rename(columns={'index': 'row_id'})

        df_obs_wnd = pd.DataFrame({
            'time': datasets[0].variables['time'][:],
            'obs_wnd': datasets[0].variables['obs_wnd_cal'][:]
        }).reset_index().rename(columns={'index': 'row_id'})

        merged_hs = pd.merge(merged_hs, df_obs_hs, on=['time', 'row_id'], how='inner')
        merged_wnd = pd.merge(merged_wnd, df_obs_wnd, on=['time', 'row_id'], how='inner')

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

        model_labels = ['Retrov17_01', 'GFSv16']

        model_columns_hs = ['hs_retrov17_01', 'hs_gfsv16']
        for col, label in zip(model_columns_hs, model_labels):
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

        model_columns_wnd = ['wnd_retrov17_01', 'wnd_gfsv16']
        for col, label in zip(model_columns_wnd, model_labels):
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

# ------------------------------ global plot ----------------------------------------------------
        for model_key, model_dir in directories.items():
            fp = os.path.join(model_dir, f"{filename}_{model_key}_{season}_{satellite_name}.nc")
            if not os.path.exists(fp):
                print(f"Missing file for global plots: {fp}")
                continue

            with nc.Dataset(fp, "r") as ds:
                lat = _get_var(ds, ["latitude", "lat"])
                lon = _get_var(ds, ["longitude", "lon"])

                model_hs = _get_var(ds, ["model_hs"])
                obs_hs = _get_var(ds, ["obs_hs_cal", "obs_hs"])

                model_wnd = _get_var(ds, ["model_wnd"])
                obs_wnd = _get_var(ds, ["obs_wnd_cal", "obs_wnd"])

            # --- Hs bias + rmse ---
            gsm = GlobalSkillMap(lat=lat, lon=lon, model=model_hs, obs=obs_hs, mlabels=[model_key])

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
                title=f"Hs Bias – {filename} ({model_key} vs {satellite_name})",
                outfile=os.path.join(output_dir, f"plot_Hs_{model_key}_{filename}_{satellite_name}_global_Bias.png"),
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
                title=f"Hs RMSE – {filename} ({model_key} vs {satellite_name})",
                outfile=os.path.join(output_dir, f"plot_Hs_{model_key}_{filename}_{satellite_name}_global_RMSE.png"),
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
                title=f"WND Bias – {filename} ({model_key} vs {satellite_name})",
                outfile=os.path.join(output_dir, f"plot_WND_{model_key}_{filename}_{satellite_name}_global_Bias.png"),
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
                title=f"WND RMSE – {filename} ({model_key} vs {satellite_name})",
                outfile=os.path.join(output_dir, f"plot_WND_{model_key}_{filename}_{satellite_name}_global_RMSE.png"),
            )

    else:
        # Not all file paths exist message
        print(f"Some files for {filename} do not exist.")

