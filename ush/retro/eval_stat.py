"""
eval_stat.py

PURPOSE:
    Load one all-in-one combined NetCDF file for each configured model and
compute verification statistics by forecast-day range for significant
wave height (Hs) and wind speed (WND).

    The script:
      1. reads configuration settings from `evalsumconfig.json`
      2. loads one combined NetCDF file for each model
      3. converts the dataset to a pandas DataFrame
      4. removes rows containing missing values
      5. converts `fcst_hr` from timedelta to hours when needed
      6. groups data into cumulative 24-hour forecast ranges:
             0–24 h, 24–48 h, ..., 360–384 h
      7. computes verification metrics using `mvalstats.metrics`
      8. plots each metric versus forecast hour for Hs and WND

    The script computes statistics using 'mvalstats' for:
      - Hs using `obs_hs`
      - WND using calibrated wind observations `obs_wnd_cal`

USAGE:
    Edit `evalsumconfig.json` to define:
      - directories    : input directories for each model
      - filenames_all  : all-in-one combined file prefix
      - satellite_name : satellite name used in the filenames and plot titles
      - season         : season label used in the filenames and plot titles
      - output_dir     : directory for output figures

    Modify the script as needed for:
      - endday        : maximum forecast day to include
      - n_filters     : number of filters to apply

OUTPUT:
    One PNG figure is created for each verification metric. Each figure
    contains two panels:
      - top panel    : Hs metric versus forecast hour
      - bottom panel : WND metric versus forecast hour

    Output filename format:

        fig_{stat_name}_{satellite_name}_{season}.png

NOTE:
    Statistics:
        The script computes the following metrics from `mvalstats.metrics`:
            bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N

    Forecast bins:
        Statistics are computed in cumulative 24-hour forecast bins:
            Day 1  :   0 < fcst_hr <=  24
            Day 2  :  24 < fcst_hr <=  48
            ...
            Day 16 : 360 < fcst_hr <= 384

    Filters:
        The script supports up to three filter levels:
            - filter 0 : all samples
            - filter 1 : obs_hs >= 4 m
            - filter 2 : obs_hs >= 7 m

    Others:
        - Hs statistics are computed from:
            model_hs vs obs_hs
        - WND statistics are plotted from:
            model_wnd vs obs_wnd_cal

AUTOR and DATE:
    03/12/2026: Ming Chen, first version

"""

import netCDF4 as nc
import numpy as np
import pandas as pd
import os
import json
import xarray as xr
import matplotlib.pyplot as plt
import mvalstats

endday     = 16
n_filters  = 1

with open('evalsumconfig.json') as config_file:
    config = json.load(config_file)

directories = config["directories"]
filename = config["filenames_all"]
satellite = config["satellite_name"]
season = config["season"]
output_dir = config["output_dir"]

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

file_paths = [os.path.join(directories[key], f"{filename}_{key}_{season}_{satellite}.nc") for key in directories]
models     = list(directories.keys())

allstats_hs      = np.zeros([len(models),n_filters,endday,9])*np.nan
allstats_wnd     = np.zeros([len(models),n_filters,endday,9])*np.nan
allstats_hs_cal  = np.zeros([len(models),n_filters,endday,9])*np.nan
allstats_wnd_cal = np.zeros([len(models),n_filters,endday,9])*np.nan

for m, path in enumerate(file_paths):
    ds = xr.open_dataset(path)
    df = ds.to_dataframe().reset_index()
    df = df.dropna()

    if df["fcst_hr"].dtype == 'timedelta64[ns]':
        df["fcst_hr"] = df["fcst_hr"].dt.total_seconds() / 3600.0

    day0 = 0
    day  = 1

    while day <= endday:
        f0 = day0*24
        f1 = day*24
        df_day = df[(df["fcst_hr"] <= f1) & (df["fcst_hr"] > f0)]

        # all Hs
        if n_filters >= 1:
            allstats_hs[m,0,day0,:] = mvalstats.metrics(df_day["model_hs"].values, df_day["obs_hs"].values)
            allstats_wnd[m,0,day0,:] = mvalstats.metrics(df_day["model_wnd"].values, df_day["obs_wnd"].values)
            allstats_hs_cal[m,0,day0,:] = mvalstats.metrics(df_day["model_hs"].values, df_day["obs_hs_cal"].values)
            allstats_wnd_cal[m,0,day0,:] = mvalstats.metrics(df_day["model_wnd"].values, df_day["obs_wnd_cal"].values)

        # Hs above 4m
        if n_filters >=2:
            df_4m = df_day[df_day["obs_hs"] >= 4]
            allstats_hs[m,1,day0,:] = mvalstats.metrics(df_4m["model_hs"].values, df_4m["obs_hs"].values)
            allstats_wnd[m,1,day0,:] = mvalstats.metrics(df_4m["model_wnd"].values, df_4m["obs_wnd"].values)

            df_4m_cal = df_day[df_day["obs_hs_cal"] >= 4]
            allstats_hs_cal[m,1,day0,:] = mvalstats.metrics(df_4m_cal["model_hs"].values, df_4m_cal["obs_hs_cal"].values)
            allstats_wnd_cal[m,1,day0,:] = mvalstats.metrics(df_4m_cal["model_wnd"].values, df_4m_cal["obs_wnd_cal"].values)

        # Hs above 7m
        if n_filters >= 3:
            df_7m = df_day[df_day["obs_hs"] >= 7]
            allstats_hs[m, 2, day0, :] = mvalstats.metrics(df_7m["model_hs"].values, df_7m["obs_hs"].values)
            allstats_wnd[m, 2, day0, :] = mvalstats.metrics(df_7m["model_wnd"].values, df_7m["obs_wnd"].values)

            df_7m_cal = df_day[df_day["obs_hs_cal"] >= 7]
            allstats_hs_cal[m, 2, day0, :] = mvalstats.metrics(df_7m_cal["model_hs"].values, df_7m_cal["obs_hs_cal"].values)
            allstats_wnd_cal[m, 2, day0, :] = mvalstats.metrics(df_7m_cal["model_wnd"].values, df_7m_cal["obs_wnd_cal"].values)

        day0 = day
        day  = day + 1

# Plot stats vs forecast hours
xday = np.arange(1, endday + 1) * 24.0
yday = np.arange(0, endday)
stats_names = ['bias', 'RMSE', 'NBias', 'NRMSE', 'SCrmse', 'SI', 'HH', 'CC', 'N']

linestyles = {0: "-", 1: "--", 2: "-."}
filter_labels = {0: "all", 1: "Hs>=4m", 2: "Hs>=7m"}
model_colors = {0: "red", 1: "black"}

for s, stat_name in enumerate(stats_names):
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9, 7), sharex=True)
    ax_hs, ax_wnd = axes

    # --- Hs ---
    for mi, model_name in enumerate(models):
        for f in range(n_filters):
            y = allstats_hs[mi, f, yday, s]
            ax_hs.plot(
                xday, y,
                color=model_colors.get(mi, "black"),
                linestyle=linestyles.get(f, "-"),
                marker="o",
                markersize=4,
                linewidth=1.5,
                label=f"{model_name} ({filter_labels.get(f, f'filter{f}')})"
            )
    ax_hs.set_title("Significant Wave Height")
    ax_hs.set_ylabel(stat_name)
    ax_hs.grid(True)
    ax_hs.set_xlim(0, endday * 24)

    # --- Wind ---
    for mi, model_name in enumerate(models):
        for f in range(n_filters):
            y = allstats_wnd_cal[mi, f, yday, s]
            ax_wnd.plot(
                xday, y,
                color=model_colors.get(mi, "black"),
                linestyle=linestyles.get(f, "-"),
                marker="o",
                markersize=4,
                linewidth=1.5,
                label=f"{model_name} ({filter_labels.get(f, f'filter{f}')})"
            )

    ax_wnd.set_title("Wind Speed (calibrated obs)")
    ax_wnd.set_xlabel("Forecast Hour")
    ax_wnd.set_ylabel(stat_name)
    ax_wnd.grid(True)
    ax_wnd.set_xlim(0, endday * 24)
    ax_wnd.set_xticks(np.arange(0, endday + 1) * 24)

    # Legends
    ax_hs.legend(loc="best", fontsize=8, framealpha=0.8)
    ax_wnd.legend(loc="best", fontsize=8, framealpha=0.8)

    supertitle = f"fig_{stat_name}_{satellite}_{season}"
    fig.suptitle(supertitle)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    outpng = os.path.join(output_dir, f"{supertitle}.png")
    fig.savefig(outpng, dpi=150)
    plt.close(fig)
