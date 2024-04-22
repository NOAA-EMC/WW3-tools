import mvalstats
import netCDF4 as nc

# Open the .nc file
nc_file = nc.Dataset('/work/noaa/marine/Ghazal.Mohammadpour/Ghazal.Mohammadpour1/hr3_eva/WW3-tools/ww3tools/output/HR3/winter/day7/HR1_filtered_day_winter.nc', 'r')

# Load model and observation data
model_data = nc_file.variables['model_hs'][:]
obs_data = nc_file.variables['obs_hs'][:]

# Close the .nc file
nc_file.close()

# Call metrics function
rmse_result = mvalstats.metrics(model_data, obs_data)

# Print error metrics
print("Error Metrics:")
print("Bias:", rmse_result[0])
print("RMSE:", rmse_result[1])
print("Normalized Bias:", rmse_result[2])
print("Normalized RMSE:", rmse_result[3])
print("SCrmse:", rmse_result[4])
print("SI:", rmse_result[5])
print("HH:", rmse_result[6])
print("CC:", rmse_result[7])

# Call smrstat function
summary_stats = mvalstats.smrstat(obs_data)

# Print summary statistics
print("\nSummary Statistics:")
print("Mean:", summary_stats[0])
print("Variance:", summary_stats[1])
print("Skewness:", summary_stats[2])
print("Kurtosis:", summary_stats[3])
print("Min:", summary_stats[4])
print("Max:", summary_stats[5])
print("80th Percentile:", summary_stats[6])
print("90th Percentile:", summary_stats[7])
print("95th Percentile:", summary_stats[8])
print("99th Percentile:", summary_stats[9])
print("99.9th Percentile:", summary_stats[10])

