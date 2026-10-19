import argparse
import os
import xarray as xr
import dask
import ww3tools.sofar_tools as sf

def main():
    parser = argparse.ArgumentParser(description="Build WW3 vs Sofar Buoy Matchup Files")
    parser.add_argument("--start-date", type=str, required=True, help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end-date", type=str, required=True, help="End date (YYYY-MM-DD)")
    parser.add_argument("--work-dir", type=str, required=True, help="Working directory for data output")
    parser.add_argument("--model-dir", type=str, required=True, help="Directory where model data to matchup is stored")
    parser.add_argument("--model-pattern", type=str, required=True, help="Pattern for model directories")
    parser.add_argument("--model-name", type=str, required=True, help="Name of the model to be saved as metadata in the matchup file")
    parser.add_argument("--model-file-template", type=str, required=True, help="relative subpath inside each cycle hour directory")
    args = parser.parse_args()

    print('\n')
    print('==========================')
    print(f"Starting pipeline for dates: {args.start_date} to {args.end_date}")

    # Ensure all SOFAR data is present or kill job
    sf.check_and_fetch_sofar_hpss(
        start_date_str=args.start_date,
        end_date_str=args.end_date,
        work_dir=args.work_dir,
        data_dir="bufr_data",
        slurm_params=None
    )

    # Ensure all the model data is present between the start and end dates
    sf.check_model_data(args.model_dir, args.start_date,
                        args.end_date, args.model_pattern)

    # Check if the bufr files have been converted to netCDF's and are all present, otherwise convert to netCDF
    sf.convert_sofar_bufr_to_nc(
                    args.work_dir,
                    bufr_subdir="bufr_data",
                    nc_subdir="nc_data",
                    start_date=args.start_date,
                    end_date=args.end_date,
                    pattern="sofar_b001_xx121_*")

    # Create a zarr store of the combined SOFAR data that can be accessed by each parallel worker
    nc_dir = os.path.join(args.work_dir, "nc_data")
    sofar_zarr_path = sf.consolidate_sofar_to_zarr(nc_dir=nc_dir, 
                                                    work_dir=args.work_dir, 
                                                    start_date_str=args.start_date, 
                                                    end_date_str=args.end_date)

    # Initialize the parallel Dask cluster
    client = sf.setup_dask_cluster()

    # Broadcast static variables to all workers to save bandwidth
    sofar_path_future = client.scatter(sofar_zarr_path, broadcast=True)
    model_name_future = client.scatter(args.model_name, broadcast=True)

    # Get initialization times and directories
    init_times, cycle_files = sf.build_cycle_inputs(start_date_str=args.start_date,
                                                    end_date_str=args.end_date,
                                                    model_base_dir=args.model_dir,
                                                    rel_path_template=args.model_file_template)

    # Test process a single cylce
    # ds_combined = sf.process_single_cycle(cycle_files[0], init_times[0],
    #                                     sofar_zarr_path, args.model_name)
    # print(ds_combined)

    print("Submitting tasks to Dask workers...")
    # Map the function across the lists of varying arguments
    futures = client.map(sf.process_single_cycle,
                         cycle_files,
                         init_times,
                         sofar_zarr_path=sofar_path_future,
                         model_name=model_name_future, 
                         parallel=True)

    # print("Computing interpolations...")

    print(f"Computing interpolations for {len(futures)} cycles...")
    
    valid_datasets = []
    completed_count = 0
    
    # as_completed yields each future the moment a worker finishes it
    for future in dask.distributed.as_completed(futures):
        completed_count += 1
        
        # Print a clean, sequential progress update to your log file
        print(f"Progress: {completed_count}/{len(futures)} cycles completed.", flush=True)
        
        # Get the actual data from the completed future
        res = future.result()
        if res is not None:
            valid_datasets.append(res)

    if not valid_datasets:
        print("Error: All forecast cycles were skipped. No data to save.")
        client.close()
        return None

    # # Gather pulls the actual datasets back from the workers to the main script
    # print("\nAll workers finished. Gathering data...")
    # results = client.gather(futures)

    # # Filter out empty results from skipped cycles
    # valid_datasets = [res for res in results if res is not None]

    # if not valid_datasets:
    #     print("Error: All forecast cycles were skipped. No data to save.")
    #     client.close()
    #     return None

    print(f"Concatenating {len(valid_datasets)} forecast cycles...")
    final_dataset = xr.concat(valid_datasets, dim="obs")
    final_dataset = final_dataset.sortby("time")
    
    matchup_dir = os.path.join(args.work_dir, "matchup_data")
    os.makedirs(matchup_dir, exist_ok=True)
    output_filename = os.path.join(matchup_dir, f"{args.model_name}_sofar_matchup_{args.start_date}.nc")
    
    print(f"Saving final matchup file to {output_filename}...")
    final_dataset.to_netcdf(output_filename)
    
    print("Pipeline successfully completed.")
    client.close()

if __name__ == "__main__":
    main()