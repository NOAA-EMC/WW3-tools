import argparse
import dask
import ww3tools.sofar_tools as sf # Import your other functions here too

def main():
    parser = argparse.ArgumentParser(description="Build WW3 vs Sofar Buoy Matchup Files")
    parser.add_argument("--start-date", type=str, required=True, help="Start date (YYYY-MM-DD)")
    parser.add_argument("--end-date", type=str, required=True, help="End date (YYYY-MM-DD)")
    parser.add_argument("--work-dir", type=str, required=True, help="Working directory for data output")
    
    args = parser.parse_args()

    print(f"Starting pipeline for dates: {args.start_date} to {args.end_date}")

    # Initialize the parallel Dask cluster using your library function
    client = sf.setup_dask_cluster()

    print("yay it is me the parallelized section")

    # Gracefully close the Dask cluster
    client.close()

if __name__ == "__main__":
    main()