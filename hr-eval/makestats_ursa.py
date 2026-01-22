import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse

import datetime as dt
from dateutil.relativedelta import relativedelta
import os
import xarray as xr

import mvalstats



'''
Calculate stats and plot
'''

def main():

  #ap = argparse.ArgumentParser()
  #ap.add_argument('-m', '--model', help="String Identifier of Model 'multi1', 'GFSv16', 'HR1', 'HR2', 'HR3a', 'HR3b'", required=True)
  #ap.add_argument('-o', '--outdir', help="Output directory for files", default='./')
  #MyArgs = ap.parse_args()

  #model = MyArgs.model
  #OUTDIR = MyArgs.outdir
  OUTDIR = f"/scratch4/NCEPDEV/marine/Ming.Chen/wave_eval/processsatdata_ursa/temp"

  #Check output directory exists:
  INPUTDIR=f"/scratch4/NCEPDEV/marine/Ming.Chen/wave_eval/processsatdata_ursa/outcombine"
  if not os.path.isdir(INPUTDIR):
    print(f"INPUTDIR ({INPUTDIR}) does not exist!!!!")
    exit(1)

  #create OUTDIR directory if it does not exist:
  if not os.path.isdir(OUTDIR):
    os.makedirs(OUTDIR)

  satelites=['JASON3', 'CRYOSAT2', 'SARAL', 'SENTINEL3A'] #JASON3,JASON2,CRYOSAT2,JASON1,HY2,SARAL,SENTINEL3A,ENVISAT,ERS1,ERS2,GEOSAT,GFO,TOPEX,SENTINEL3B,CFOSAT

  #seasons[k],satelites,model, (all, above 4hs, above 7hs),day, stats:bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N
  allstats_hs=np.zeros([3,4,6,3,16,9])*np.nan
  allstats_wnd=np.zeros([3,4,6,3,16,9])*np.nan
  allstats_hs_cal=np.zeros([3,4,6,3,16,9])*np.nan
  allstats_wnd_cal=np.zeros([3,4,6,3,16,9])*np.nan

  season=['winter']

  for k in range(len(season)):
    if season[k] == "winter":
       startdate = dt.datetime(2024,12,1,0)
       enddate = dt.datetime(2024,12,17,18)
       datestride = 6
       endday = 16
       model=['Retrotests', 'GFSv16']
    elif season[k] == "summer":
       startdate = dt.datetime(2020,6,1)
       enddate = dt.datetime(2020,8,30)
       datestride = 3
       endday = 16
       model=['multi1','HR1', 'HR2', 'HR3a', 'HR3b', 'GFSv16']
    elif season[k] == "hurricane":
       startdate = dt.datetime(2020,7,20)
       enddate = dt.datetime(2020,11,20)
       datestride = 1
       endday = 7
       model=['multi1', 'HR1', 'HR2', 'HR3a', 'HR3b', 'GFSv16']

    for j in range(len(satelites)):
      for m in range(len(model)):
        time = []; lats = []; lons = []
        fhrs = []; fhrsall = [];
        obs_hs = []; obs_wnd = []
        model_hs = []; model_wnd = []
        obs_hs_cal = []; obs_wnd_cal = []

        INPUT_FILE=f"combined_{model[m]}_{season[k]}_{satelites[j]}.nc"

        datapath = INPUTDIR + "/" + INPUT_FILE
        datanc  = nc.Dataset(datapath)

        time = np.array(datanc.variables['time'][:])
        fhrs = np.array(datanc.variables['fcst_hr'][:])
        lats = np.array(datanc.variables['latitude'][:])
        lons = np.array(datanc.variables['longitude'][:])
        obs_hs = np.array(datanc.variables['obs_hs'][:])
        obs_hs_cal = np.array(datanc.variables['obs_hs_cal'][:])
        obs_wnd = np.array(datanc.variables['obs_wnd'][:])
        obs_wnd_cal = np.array(datanc.variables['obs_wnd_cal'][:])
        model_hs = np.array(datanc.variables['model_hs'][:])
        model_wnd = np.array(datanc.variables['model_wnd'][:])

        datanc.close()

        day0=0
        day=1
        endday2=endday
        if model[m] == "multi1":
           endday2 = 7

        while day <= endday2:
          f0 = day0*24
          f1 = day*24
          indx=np.where(( fhrs <= f1 ) & ( fhrs > f0 ))

          #seasons[k],satelites,model, (all, above 4hs, above 7hs), stats:bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N
          allstats_hs[k,j,m,0,day0,:] = mvalstats.metrics(model_hs[indx],obs_hs[indx]) 
          allstats_wnd[k,j,m,0,day0,:] = mvalstats.metrics(model_wnd[indx],obs_wnd[indx])
          allstats_hs_cal[k,j,m,0,day0,:] = mvalstats.metrics(model_hs[indx],obs_hs_cal[indx])
          allstats_wnd_cal[k,j,m,0,day0,:] = mvalstats.metrics(model_wnd[indx],obs_wnd_cal[indx]) 

          #above 4m 
          indx=np.where( (obs_hs >= 4) & ( fhrs <= f1 ) & ( fhrs > f0 ) )
          allstats_hs[k,j,m,1,day0,:] = mvalstats.metrics(model_hs[indx],obs_hs[indx])
          allstats_wnd[k,j,m,1,day0,:] = mvalstats.metrics(model_wnd[indx],obs_wnd[indx])

          indx=np.where( (obs_hs_cal >= 4) & ( fhrs <= f1 ) & ( fhrs > f0 ) )
          allstats_hs_cal[k,j,m,1,day0,:] = mvalstats.metrics(model_hs[indx],obs_hs_cal[indx])
          allstats_wnd_cal[k,j,m,1,day0,:] = mvalstats.metrics(model_wnd[indx], obs_wnd_cal[indx])

          #above 7m 
          indx=np.where( (obs_hs >= 4) & ( fhrs <= f1 ) & ( fhrs > f0 ) )
          allstats_hs[k,j,m,1,day0,:] = mvalstats.metrics(model_hs[indx],obs_hs[indx])
          allstats_wnd[k,j,m,1,day0,:] = mvalstats.metrics(model_wnd[indx],obs_wnd[indx])

          indx=np.where( (obs_hs_cal >= 4) & ( fhrs <= f1 ) & ( fhrs > f0 ) )
          allstats_hs_cal[k,j,m,1,day0,:] = mvalstats.metrics(model_hs[indx],obs_hs_cal[indx])
          allstats_wnd_cal[k,j,m,1,day0,:] = mvalstats.metrics(model_wnd[indx], obs_wnd_cal[indx])


          day0 = day
          day = day + 1



  #seasons[k],satelites,model, (all, above 4hs, above 7hs),day, stats:bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N
  ncout = nc.Dataset('myfile.nc','w','NETCDF4'); # using netCDF3 for output format
  ncout.createDimension('season',3);
  ncout.createDimension('sat',4);
  ncout.createDimension('model',6);
  ncout.createDimension('maxhsval',3);
  ncout.createDimension('days',16);
  ncout.createDimension('stats',9);
  myvar1 = ncout.createVariable('allstats_hs','float32',('season','sat','model','maxhsval','days','stats'))
  myvar1[:] = allstats_hs;
  myvar2 = ncout.createVariable('allstats_wnd','float32',('season','sat','model','maxhsval','days','stats'))
  myvar2[:] = allstats_wnd;
  myvar3 = ncout.createVariable('allstats_hs_cal','float32',('season','sat','model','maxhsval','days','stats'))
  myvar3[:] = allstats_hs_cal;
  myvar4 = ncout.createVariable('allstats_wnd_cal','float32',('season','sat','model','maxhsval','days','stats'))
  myvar4[:] = allstats_wnd_cal;
  ncout.close();



if __name__ == '__main__':
    main()
