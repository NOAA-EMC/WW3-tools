import netCDF4 as nc
import numpy as np
import os
import json
import matplotlib.pyplot as plt
import mvalstats
from pvalstats import GlobalSkillMap

# -------------------- global map configuration --------------------
DLAT = 1.0
DLON = 1.0
MIN_COUNT = 10
LATMIN = -60.0
LATMAX = 60.0

# Color scale limits
HS_BIAS_VMAX = 0.5
HS_RMSE_VMAX = 1.0
WND_BIAS_VMAX = 2.0
WND_RMSE_VMAX = 5.0

# QC ranges
QC_HS  = {"model_min": 0.0, "model_max": 30.0, "obs_min": 0.0, "obs_max": 30.0}
QC_WND = {"model_min": 0.0, "model_max": 60.0, "obs_min": 0.0, "obs_max": 60.0}
# ------------------------------------------------------------------
def _get_var(ds, candidates):
    for v in candidates:
        if v in ds.variables:
            return ds.variables[v][:]
    raise KeyError(f"Missing variable. Tried: {candidates}")

def _to_hours(fcst_hr):
    """
    Try to convert fcst_hr to numeric hours for comparisons.
    Handles common cases:
      - numeric hours already
      - timedelta64 arrays (convert to hours)
      - netCDF masked arrays
    """
    fcst_hr = np.asarray(fcst_hr)

    # netCDF4 sometimes returns masked arrays
    if np.ma.isMaskedArray(fcst_hr):
        fcst_hr = fcst_hr.filled(np.nan)

    # timedelta64 -> hours
    if np.issubdtype(fcst_hr.dtype, np.timedelta64):
        return fcst_hr / np.timedelta64(1, "h")

    # Otherwise assume numeric already (hours)
    return fcst_hr.astype(float)

def _metrics_safe(mod, obs, n_metrics=9):
    """
    NaN-safe wrapper:
      - mask non-finite pairs
      - if empty, return all-NaNs (length n_metrics)
      - else return mvalstats.metrics(...)
    """
    mod = np.asarray(mod)
    obs = np.asarray(obs)

    if np.ma.isMaskedArray(mod):
        mod = mod.filled(np.nan)
    if np.ma.isMaskedArray(obs):
        obs = obs.filled(np.nan)

    mask = np.isfinite(mod) & np.isfinite(obs)
    if mask.sum() == 0:
        return np.full(n_metrics, np.nan, dtype=float)

    out = mvalstats.metrics(mod[mask], obs[mask])
    out = np.asarray(out, dtype=float)

    # Guard in case metrics returns list/tuple of different length
    if out.size != n_metrics:
        raise ValueError(f"metrics() returned length {out.size}, expected {n_metrics}")
    return out

def main(config_path="evalsumconfig.json"):

    with open(config_path) as config_file:
        config = json.load(config_file)

    directories = config["directories"]
    filenames_all = config["filenames_all"]
    satellite_name = config["satellite_name"]
    season = config["season"]
    output_dir = config["output_dir"]

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    models = list(directories.keys())
    n_models = len(models)

    endday = 16
    n_filters = 3
    n_metrics = 9   # bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N
    stats_names = ['bias', 'RMSE', 'NBias', 'NRMSE', 'SCrmse', 'SI', 'HH', 'CC', 'N']

    allstats_hs      = np.full((n_models, n_filters, endday, n_metrics), np.nan, dtype=float)
    allstats_wnd     = np.full((n_models, n_filters, endday, n_metrics), np.nan, dtype=float)
    allstats_hs_cal  = np.full((n_models, n_filters, endday, n_metrics), np.nan, dtype=float)
    allstats_wnd_cal = np.full((n_models, n_filters, endday, n_metrics), np.nan, dtype=float)

    file_paths = [
        os.path.join(directories[key], f"{filenames_all}_{key}_{season}_{satellite_name}.nc")
        for key in models
    ]

    # Compute stats arrays
    for mi, (model_name, fp) in enumerate(zip(models, file_paths)):
        if not os.path.exists(fp):
            print(f"[WARN] Missing input: {fp}")
            continue

        with nc.Dataset(fp, "r") as ds:
            fcst_hr = _to_hours(ds.variables["fcst_hr"][:])

            model_hs  = ds.variables["model_hs"][:]
            model_wnd = ds.variables["model_wnd"][:]

            obs_hs      = ds.variables["obs_hs"][:]
            obs_wnd     = ds.variables["obs_wnd"][:]
            obs_hs_cal  = ds.variables["obs_hs_cal"][:]
            obs_wnd_cal = ds.variables["obs_wnd_cal"][:]

        for day0 in range(endday):
            f0 = day0 * 24.0
            f1 = (day0 + 1) * 24.0

            # Forecast-hour window (same rule as your original code: (f0, f1])
            idx = (fcst_hr > f0) & (fcst_hr <= f1)
            if np.count_nonzero(idx) == 0:
                continue

            # filter 0: all
            allstats_hs[mi, 0, day0, :]      = _metrics_safe(model_hs[idx],  obs_hs[idx], n_metrics)
            allstats_wnd[mi, 0, day0, :]     = _metrics_safe(model_wnd[idx], obs_wnd[idx], n_metrics)
            allstats_hs_cal[mi, 0, day0, :]  = _metrics_safe(model_hs[idx],  obs_hs_cal[idx], n_metrics)
            allstats_wnd_cal[mi, 0, day0, :] = _metrics_safe(model_wnd[idx], obs_wnd_cal[idx], n_metrics)

            # filter 1: Hs >= 4m (use obs thresholds separately for raw/cal)
            idx4  = idx & (np.asarray(obs_hs) >= 4.0)
            idx4c = idx & (np.asarray(obs_hs_cal) >= 4.0)

            allstats_hs[mi, 1, day0, :]      = _metrics_safe(model_hs[idx4],   obs_hs[idx4], n_metrics)
            allstats_wnd[mi, 1, day0, :]     = _metrics_safe(model_wnd[idx4],  obs_wnd[idx4], n_metrics)
            allstats_hs_cal[mi, 1, day0, :]  = _metrics_safe(model_hs[idx4c],  obs_hs_cal[idx4c], n_metrics)
            allstats_wnd_cal[mi, 1, day0, :] = _metrics_safe(model_wnd[idx4c], obs_wnd_cal[idx4c], n_metrics)

            # filter 2: Hs >= 7m
            idx7  = idx & (np.asarray(obs_hs) >= 7.0)
            idx7c = idx & (np.asarray(obs_hs_cal) >= 7.0)

            allstats_hs[mi, 2, day0, :]      = _metrics_safe(model_hs[idx7],   obs_hs[idx7], n_metrics)
            allstats_wnd[mi, 2, day0, :]     = _metrics_safe(model_wnd[idx7],  obs_wnd[idx7], n_metrics)
            allstats_hs_cal[mi, 2, day0, :]  = _metrics_safe(model_hs[idx7c],  obs_hs_cal[idx7c], n_metrics)
            allstats_wnd_cal[mi, 2, day0, :] = _metrics_safe(model_wnd[idx7c], obs_wnd_cal[idx7c], n_metrics)

    print("Done computing allstats arrays.")

    # Plot stats vs forecast hours
    xday = np.arange(1, endday + 1) * 24.0   # 24..384
    yday = np.arange(0, endday)

    # Line styles by filter (solid for all, dashed for >=4, dashdot for >=7)
    linestyles = {0: "-", 1: "--", 2: "-."}
    filter_labels = {0: "all", 1: "Hs>=4m", 2: "Hs>=7m"}
    model_colors = {0: "blue", 1: "red"}

    for s, stat_name in enumerate(stats_names):
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(9, 7), sharex=True)
        ax_hs, ax_wnd = axes

        # --- Hs ---
        for mi, model_name in enumerate(models):
            for f in range(n_filters):
                y = allstats_hs_cal[mi, f, yday, s]
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

        supertitle = f"fig_{stat_name}_{satellite_name}_{season}"
        fig.suptitle(supertitle)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])

        outpng = os.path.join(output_dir, f"{supertitle}.png")
        fig.savefig(outpng, dpi=150)
        plt.close(fig)

    print(f"[INFO] Plots saved under: {output_dir}")

    # ------------------------------ global plots (Bias/RMSE maps) ------------------------------
    for mi, model_key in enumerate(models):
        fp = os.path.join(directories[model_key],
                          f"{filenames_all}_{model_key}_{season}_{satellite_name}.nc")

        if not os.path.exists(fp):
            print(f"[WARN] Missing file for global plots: {fp}")
            continue

        with nc.Dataset(fp, "r") as ds:
            lat = _get_var(ds, ["latitude", "lat"])
            lon = _get_var(ds, ["longitude", "lon"])

            model_hs  = _get_var(ds, ["model_hs"])
            model_wnd = _get_var(ds, ["model_wnd"])

            obs_hs  = _get_var(ds, ["obs_hs_cal", "obs_hs"])
            obs_wnd = _get_var(ds, ["obs_wnd_cal", "obs_wnd"])

        # ---- Hs maps ----
        gsm_hs = GlobalSkillMap(lat=lat, lon=lon, model=model_hs, obs=obs_hs, mlabels=[model_key])

        gsm_hs.plot_global_like_global_plot(
            metric="bias",
            model_index=0,
            dlat=DLAT, dlon=DLON,
            lon_0_360=True,
            min_count=MIN_COUNT,
            latmin=LATMIN, latmax=LATMAX,
            qc_kwargs=QC_HS,
            vmax=HS_BIAS_VMAX,
            title=f"Hs Bias – {model_key} vs {satellite_name}",
            outfile=os.path.join(output_dir,
                f"plot_HS_{model_key}_{filenames_all}_{season}_{satellite_name}_global_Bias.png"),
        )

        gsm_hs.plot_global_like_global_plot(
            metric="rmse",
            model_index=0,
            dlat=DLAT, dlon=DLON,
            lon_0_360=True,
            min_count=MIN_COUNT,
            latmin=LATMIN, latmax=LATMAX,
            qc_kwargs=QC_HS,
            vmax=HS_RMSE_VMAX,
            title=f"Hs RMSE – {model_key} vs {satellite_name}",
            outfile=os.path.join(output_dir,
                f"plot_HS_{model_key}_{filenames_all}_{season}_{satellite_name}_global_RMSE.png"),
        )

        # ---- Wind maps ----
        gsm_wnd = GlobalSkillMap(lat=lat, lon=lon, model=model_wnd, obs=obs_wnd, mlabels=[model_key])

        gsm_wnd.plot_global_like_global_plot(
            metric="bias",
            model_index=0,
            dlat=DLAT, dlon=DLON,
            lon_0_360=True,
            min_count=MIN_COUNT,
            latmin=LATMIN, latmax=LATMAX,
            qc_kwargs=QC_WND,
            vmax=WND_BIAS_VMAX,
            title=f"WND Bias – {model_key} vs {satellite_name}",
            outfile=os.path.join(output_dir,
                f"plot_WND_{model_key}_{filenames_all}_{season}_{satellite_name}_global_Bias.png"),
        )

        gsm_wnd.plot_global_like_global_plot(
            metric="rmse",
            model_index=0,
            dlat=DLAT, dlon=DLON,
            lon_0_360=True,
            min_count=MIN_COUNT,
            latmin=LATMIN, latmax=LATMAX,
            qc_kwargs=QC_WND,
            vmax=WND_RMSE_VMAX,
            title=f"WND RMSE – {model_key} vs {satellite_name}",
            outfile=os.path.join(output_dir,
                f"plot_WND_{model_key}_{filenames_all}_{season}_{satellite_name}_global_RMSE.png"),
        )

    print("[INFO] Global maps saved under {output_dir}.")

if __name__ == "__main__":
    main("evalsumconfig.json")
