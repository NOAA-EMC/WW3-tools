import os
import shutil
import subprocess
from datetime import datetime, timedelta
"""
Module for fetching Sofar Spotter buoy data from the NOAA HPSS system.
"""

def get_sofar_hpss(start_date: str, end_date: str, working_dir: str) -> None:
    """
    Retrieves Sofar HPSS data for a given date range.
    
    Args:
        start_date (str): The start date in 'YYYY-MM-DD' format.
        end_date (str): The end date in 'YYYY-MM-DD' format.
        working_dir (str): The target directory for data extraction.
    """
    # Ensure the working directory exists and navigate to it
    os.makedirs(working_dir, exist_ok=True)
    os.chdir(working_dir)

    # Convert string inputs to datetime objects
    start_dt = datetime.strptime(start_date, "%Y-%m-%d")
    end_dt = datetime.strptime(end_date, "%Y-%m-%d")

    print(f"start date: {start_dt.strftime('%Y-%m-%d %H')}")
    print(f"end date: {end_dt.strftime('%Y-%m-%d %H')}")
    print("retrieving from hpss/dcom")

    print("*** This is being built to make sure it connects with hpss through the correct node")
    # # Iterate through the dates day-by-day
    # current_date = start_dt
    # while current_date <= end_dt:
    #     yy = current_date.strftime("%Y")
    #     mm = current_date.strftime("%m")
    #     dd = current_date.strftime("%d")
        
    #     # Define the archive path and the local extraction directory
    #     file_atm = f"/NCEPPROD/hpssprod/runhistory/rh{yy}/{yy}{mm}/{yy}{mm}{dd}/dcom_{yy}{mm}{dd}.tar"
    #     atmos = "./b001"
        
    #     try:
    #         # 1. Run the htar extraction command
    #         subprocess.run(["htar", "-P", "-xvf", file_atm, atmos], check=True)
            
    #         # 2. Move and rename the target file
    #         src_file = os.path.join(atmos, "xx121")
    #         dest_file = f"sofar_b001_xx121_{yy}{mm}{dd}"
            
    #         if os.path.exists(src_file):
    #             shutil.move(src_file, dest_file)
    #             print(f"Successfully retrieved and renamed: {dest_file}")
    #         else:
    #             print(f"Warning: Extracted file {src_file} not found for {yy}-{mm}-{dd}")
                
    #     except subprocess.CalledProcessError as e:
    #         print(f"Error executing htar for {file_atm}: {e}")
            
    #     finally:
    #         # 3. Clean up the extraction directory (equivalent to rm -rf)
    #         if os.path.exists(atmos):
    #             shutil.rmtree(atmos)
                
    #     # Increment by 1 day
    #     current_date += timedelta(days=1)