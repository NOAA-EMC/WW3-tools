import os

def plot_nc_files(directory_path, variable_name):
    # Check if directory exists
    if not os.path.exists(directory_path):
        print(f"Directory {directory_path} does not exist.")
        return

    # Iterate through files in directory
    for filename in os.listdir(directory_path):
        if filename.endswith(".nc"):
            file_path = os.path.join(directory_path, filename)
            # Generate plot for the file
            command = f"python3 ww3fields.py {file_path} {variable_name}"
            os.system(command)

# Example usage:
directory_path = "/work/noaa/marine/jmeixner/Data/HR1/summer/gfs.20200601/00/wave/gridded"
variable_name = "HTSGW_surface"  # Example variable name
plot_nc_files(directory_path, variable_name)

