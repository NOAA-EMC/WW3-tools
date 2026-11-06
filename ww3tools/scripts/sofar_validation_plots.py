import argparse
import matplotlib.pyplot as plt
import yaml
import numpy as np
import xarray as xr

from ww3tools import sofar_tools as sf

def main():
    parser = argparse.ArgumentParser(description="Generate SOFAR matchup plots from YAML config.")
    parser.add_argument('--config', required=True, help="Path to the YAML configuration file")
    args = parser.parse_args()

    # Load the yaml configuration inputs
    with open(args.config, 'r') as file:
        config = yaml.safe_load(file)

    # Get input values from the cofiguration
    labels = list(config['models'].keys())
    files = list(config['models'].values())
    variables = config['variables']
    stats = config['statistics']
    window = config.get('window_hours', 1.0) # Defaults to 1.0 if missing
    outdir = config['outdir']

    # Print out configuration
    print(f"Loaded config: Comparing {labels}")
    print(f"matchup files: {files}")
    print(f"variables: {variables}")
    print(f"stats: {stats}")
    print(f"window time: {window}")
    print(f"output directory: {outdir}")
    
    # Load the matchup datasets
    datasets = []
    for file_path in files:
        if file_path.endswith('.zarr'):
            ds = xr.open_zarr(file_path, decode_timedelta=False)
        else:
            ds = xr.open_dataset(file_path)
        datasets.append(ds)

    # Compute the windowed statistics
    print("Computing statistics...")
    computed_stats_ds = sf.compute_windowed_statistics(
        datasets=datasets,
        variables=variables,
        labels=labels,
        statistics=stats,
        window_hours=window
    )

    # Plot the windowed statistics
    print(f"Generating plots and saving to {outdir}...")
    sf.plot_statistics(
        computed_ds=computed_stats_ds,
        variables=variables,
        labels=labels,
        statistics=stats,
        window_hours=window,
        save_dir=outdir
    )
    
    print("Done!")

if __name__ == "__main__":
    main()