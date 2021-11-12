# Time series plots of WAVEWATCHIII results.
# Mandatory Inputs: fileName(netcdf format) and StationName(or ID, starting with 0)
# if you want to change the resolution (for publications), edit savefig
# examples:
# python3 ww3pointimeseries.py ww3gefs.20160928_tab.nc 41002
# nohup python3 ww3pointimeseries.py ww3gefs.20160928_tab.nc 41002 >> nohup_ww3pointimeseries_20160928.txt 2>&1 &

import matplotlib
matplotlib.use('Agg')
import xarray as xr
import numpy as np
from pylab import *
from matplotlib.mlab import *
import matplotlib.pyplot as plt
import sys
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

# Input arguments
if len(sys.argv) < 3 :
	sys.exit(' Two inputs must be provided: fileName and StationName')
elif len(sys.argv) == 3:
	fname=np.str(sys.argv[1]); stname=np.str(sys.argv[2])
elif len(sys.argv) > 3:
	sys.exit(' Too many inputs')

# --- READ FILE ---
# NetCDF format
ds = xr.open_dataset(fname)
auxstationname=ds['station_name'].values[:,:]; stationname=[]
for i in range(0,auxstationname.shape[0]):
	stationname=np.append(stationname,"".join(np.array(auxstationname[i,:]).astype('str')))


if str.isnumeric(stname) == True:
	if np.int(stname) < 1000:
		inds=np.int(stname); stname=np.str(stationname[inds])
	else:
		inds=np.where(stationname[:]==stname)
		if size(inds)>0:
			inds=np.int(inds[0][0]); stname=np.str(stationname[inds])
		else:
			sys.exit(' Station '+stname+' not included in the output file, or wrong station ID')


ntime = ds['time'].values[:]

# PLOT ---
wvars=list(ds.keys())[3::]
fig = plt.figure(figsize=(15,7))
for i in range(0,size(wvars)):
	fpos=421+i
	ax = fig.add_subplot(fpos)
	ax.plot_date(ntime,ds[wvars[i]].values[:,inds],'k.')
	ax.set_xlim( ntime[0], ntime[-1] )
	ax.xaxis.set_major_formatter( DateFormatter('%b%d') )
	ax.fmt_xdata = DateFormatter('%b%d')
	plt.ylabel(wvars[i]+' ('+ds[wvars[i]].units+')', fontsize=9)
	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(7)

	plt.grid(linewidth=0.5, color='grey', alpha=0.5, linestyle='--')

fig.tight_layout()

plt.savefig('wtimeseries_'+stname+'_'+np.str(pd.to_datetime(ntime[0]).strftime('%Y%m%d%H'))+'.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', 
	papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)

plt.close(); ds.close()
del ax, fig, ds

