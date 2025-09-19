import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse

import datetime as dt
from dateutil.relativedelta import relativedelta
import os
import xarray as xr


'''
Create combined NetCDF files for easier post processing. 
'''

def main():

  ap = argparse.ArgumentParser()
  ap.add_argument('-m', '--model', help="String Identifier of Model 'multi1', 'GFSv16', 'HR1', 'HR2', 'HR3a', 'HR3b'", required=True)
  ap.add_argument('-o', '--outdir', help="Output directory for files", default='./')
  MyArgs = ap.parse_args()

  model = MyArgs.model
  OUTDIR = MyArgs.outdir


  #Check output directory exists: 
  INPUTDIR=f"/work2/noaa/marine/jmeixner/processsatdata/outinterp/{model}"
  if not os.path.isdir(INPUTDIR):
    INPUTDIR=f"/scratch1/NCEPDEV/climate/Jessica.Meixner/processsatdata/outinterp/{model}"
    if not os.path.isdir(INPUTDIR):
      print('INPUTDIR ({INPUTDIR}) does not exist!!!!') 
      exit(1) 

  #create OUTDIR directory if it does not exist: 
  if not os.path.isdir(OUTDIR):
    os.makedirs(OUTDIR)



  satelites=['JASON3', 'CRYOSAT2', 'SARAL', 'SENTINEL3A'] #JASON3,JASON2,CRYOSAT2,JASON1,HY2,SARAL,SENTINEL3A,ENVISAT,ERS1,ERS2,GEOSAT,GFO,TOPEX,SENTINEL3B,CFOSAT

  if model == "GFSv16": 
    season=['summer', 'hurricane']
  else: 
    season=['winter', 'summer', 'hurricane']

  for k in range(len(season)):
    if season[k] == "winter":
       startdate = dt.datetime(2019,12,3)
       enddate = dt.datetime(2020,2,26)
       datestride = 3 
       endday = 16
    elif season[k] == "summer":
       startdate = dt.datetime(2020,6,1)
       enddate = dt.datetime(2020,8,30)
       datestride = 3
       endday = 16
    elif season[k] == "hurricane":
       startdate = dt.datetime(2020,7,20)
       enddate = dt.datetime(2020,11,20)
       datestride = 1
       endday = 7

    #if multi-1 end at day 7
    if model == "multi1":
      endday = 7


    nowdate = startdate
    dates1 = []
    while nowdate <= enddate:
       dates1.append(nowdate.strftime('%Y%m%d%H'))
       nowdate = nowdate + dt.timedelta(days=datestride)

    for j in range(len(satelites)): 
      time = []; lats = []; lons = []
      fhrs = []; fhrsall = [];
      obs_hs = []; obs_wnd = []
      model_hs = []; model_wnd = []
      obs_hs_cal = []; obs_wnd_cal = []
      for i in range(len(dates1)):
         #list of grids for each model.  First should be "global" or the base, followed by high resolution inserts in the order 
         #of lower(global) to higher(regional) resolution. 
         if model == "multi1":
            grids=['global.0p50', 'alaska.0p16', 'atlocn.0p16', 'epacif.0p16', 'wcoast.0p16', 'alaska.0p06', 'atlocn.0p06', 'wcoast.0p06']
         elif model == "GFSv16":
            grids=['global.0p25', 'global.0p16']
         else:
            grids=['global.0p25']

         for g in range(len(grids)): 
            
            if model == "multi1": 
               INPUT_FILE=f"{model}_{grids[g]}_{dates1[i]}_{satelites[j]}.nc"
            elif model == "GFSv16":
               INPUT_FILE=f"{model}_{grids[g]}_{dates1[i]}_{satelites[j]}.nc"
            else: 
               INPUT_FILE=f"{model}_{grids[g]}_{season[k]}_{dates1[i]}_{satelites[j]}.nc"

            datapath = INPUTDIR + "/" + INPUT_FILE
            datanc  = nc.Dataset(datapath)
            if g == 0:
               #this is the global/base grid:  
               time_tmpbase = np.array(datanc.variables['time'][:])
               fhrs_tmpbase = np.array(datanc.variables['fcst_hr'][:])
               lats_tmpbase = np.array(datanc.variables['latitude'][:])
               lons_tmpbase = np.array(datanc.variables['longitude'][:])
               obs_hs_tmpbase = np.array(datanc.variables['obs_hs'][:])
               obs_hs_cal_tmpbase = np.array(datanc.variables['obs_hs_cal'][:])
               obs_wnd_tmpbase = np.array(datanc.variables['obs_wnd'][:])
               obs_wnd_cal_tmpbase = np.array(datanc.variables['obs_wnd_cal'][:])
               model_hs_tmpbase = np.array(datanc.variables['model_hs'][:])
               model_wnd_tmpbase = np.array(datanc.variables['model_wnd'][:]) 
               initial_condition_time = datanc.getncattr('initial_condition_time')
               fhrsall_tmpbase = (time_tmpbase - initial_condition_time)/3600 
            else: 
               #this is s a higher resolution sub-grid 
               time_tmphigh = np.array(datanc.variables['time'][:])
               fhrs_tmphigh = np.array(datanc.variables['fcst_hr'][:])
               lats_tmphigh = np.array(datanc.variables['latitude'][:])
               lons_tmphigh = np.array(datanc.variables['longitude'][:])
               obs_hs_tmphigh = np.array(datanc.variables['obs_hs'][:])
               obs_hs_cal_tmphigh = np.array(datanc.variables['obs_hs_cal'][:])
               obs_wnd_tmphigh = np.array(datanc.variables['obs_wnd'][:])
               obs_wnd_cal_tmphigh = np.array(datanc.variables['obs_wnd_cal'][:])
               model_hs_tmphigh = np.array(datanc.variables['model_hs'][:])
               model_wnd_tmphigh = np.array(datanc.variables['model_wnd'][:])
               #Check that obs values are the same for sanity check and if so, 
               #replace model values with high res inserts where HS is not nan 
               if ((obs_hs_tmphigh == obs_hs_tmpbase).all()): 
                  np.where(~np.isnan(model_hs_tmphigh), model_hs_tmpbase, model_hs_tmphigh)
                  np.where(~np.isnan(model_hs_tmphigh), model_wnd_tmpbase, model_wnd_tmphigh)


         time = np.append(time, time_tmpbase) 
         lats = np.append(lats, lats_tmpbase)
         lons = np.append(lons, lons_tmpbase)
         fhrs = np.append(fhrs, fhrs_tmpbase)
         fhrsall = np.append(fhrsall, fhrsall_tmpbase)

         obs_hs = np.append(obs_hs, obs_hs_tmpbase)
         obs_wnd = np.append(obs_wnd, obs_wnd_tmpbase)
         obs_hs_cal = np.append(obs_hs_cal, obs_hs_cal_tmpbase)
         obs_wnd_cal = np.append(obs_wnd_cal, obs_wnd_cal_tmpbase)

         model_hs = np.append(model_hs, model_hs_tmpbase)
         model_wnd = np.append(model_wnd, model_wnd_tmpbase) 
  
      #remove values we should not use for wind due to zero at fhr0
      if model == "HR1":
         model_wnd[fhrs<3]=np.nan 
      elif model == "HR2":
         model_wnd[fhrs<3]=np.nan 
      elif model == "HR3a":
         model_wnd[fhrs<1]=np.nan 
      elif model == "HR3b":
         model_wnd[fhrs<1]=np.nan 

      #Call function to write out netcdf file with all forecast hours
      outfilename=f"combined_{model}_{season[k]}_{satelites[j]}.nc"
      OUTFILE = OUTDIR + '/' + outfilename 
      write_netcdf_file(OUTFILE, model, satelites[j], time,lats, lons, fhrs, obs_hs, obs_hs_cal, obs_wnd, obs_wnd_cal, model_hs, model_wnd)
   
      day0=0
      day=1 
      while day <= endday:
        f0 = day0*24 
        f1 = day*24
        #it will likely be easier to match all models up if we don't filter nans out here... 
        #indx=np.where(( fhrs < f1 ) & ( fhrs > f0 ) & (~np.isnan(model_hs))) 
        indx=np.where(( fhrsall <= f1 ) & ( fhrsall > f0 )) 
        time_day = time[indx]
        lats_day = lats[indx] 
        lons_day = lons[indx] 
        fhrs_day = fhrs[indx]
        obs_hs_day = obs_hs[indx]
        obs_wnd_day = obs_wnd[indx]
        obs_hs_cal_day = obs_hs_cal[indx]
        obs_wnd_cal_day = obs_wnd_cal[indx]
        model_hs_day = model_hs[indx]
        model_wnd_day = model_wnd[indx]
      
        #Call function to write out netcdf file for each day
        outfilename=f"combined_day{day:02d}_{model}_{season[k]}_{satelites[j]}.nc"
        OUTFILE = OUTDIR + '/' + outfilename
        write_netcdf_file(OUTFILE, model, satelites[j], time_day,lats_day, lons_day, fhrs_day, obs_hs_day, obs_hs_cal_day, obs_wnd_day, obs_wnd_cal_day, model_hs_day, model_wnd_day)

        day0 = day
        day = day + 1

def write_netcdf_file(nameoffile, nameofmodel,nameofsat, val_time, val_lats, val_lons, val_fhrs, val_obs_hs, val_obs_hs_cal, val_obs_wnd, val_obs_wnd_cal, val_model_hs, val_model_wnd): 

        time_dataarray = xr.DataArray(val_time, dims=['time'], name='time', attrs={
           'standard_name': 'time',
           'units': 'seconds since 1970-01-01 00:00:00',
           'calendar': 'standard',
           'axis': 'T'
        })

        interpolated_dataset = xr.Dataset({
           'time': time_dataarray,
           'latitude': xr.DataArray(val_lats, coords={'time': val_time}, dims=['time'], name='latitude').assign_attrs(units='degree_north'),
           'longitude': xr.DataArray(val_lons, coords={'time': val_time}, dims=['time'], name='longitude').assign_attrs(units='degree_east'),
           'model_hs': xr.DataArray(val_model_hs, coords={'time': val_time}, dims=['time'], name='model_hs').assign_attrs(units='m'),
           'model_wnd': xr.DataArray(val_model_wnd, coords={'time': val_time}, dims=['time'], name='model_wnd').assign_attrs(units='m'),
           'obs_hs': xr.DataArray(val_obs_hs, coords={'time': val_time}, dims=['time'], name='obs_hs').assign_attrs(units='m'),
           'obs_hs_cal': xr.DataArray(val_obs_hs_cal, coords={'time': val_time}, dims=['time'], name='obs_hs_cal').assign_attrs(units='m'),
           'obs_wnd': xr.DataArray(val_obs_wnd, coords={'time': val_time}, dims=['time'], name='obs_wnd').assign_attrs(units='m/s'),
           'obs_wnd_cal': xr.DataArray(val_obs_wnd_cal, coords={'time': val_time}, dims=['time'], name='obs_wnd_cal').assign_attrs(units='m/s'),
           'fcst_hr': xr.DataArray(val_fhrs,coords={'time': val_time}, dims=['time'], name='fcst_hr').assign_attrs(description="Forecast hour relative to initial condition time", units='hours')
        })

        interpolated_dataset.attrs['satellite_name'] = f"{nameofsat}"
        interpolated_dataset.attrs['model_name'] = f"{nameofmodel}"
        interpolated_dataset.to_netcdf(nameoffile, format='NETCDF4')


if __name__ == '__main__':
    main()
