#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
get_ndbc_stdmet.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022

PURPOSE:
 Download the historical NDBC buoy data in stdmet format
  https://www.ndbc.noaa.gov/rsa.shtml
  https://www.ndbc.noaa.gov/measdes.shtml
  https://www.ndbc.noaa.gov/stndesc.shtml
 The url used for download refers to quality-controlled observations,
  and the code is meant for historical data (archive). Recent data (last
  week or month) is not included.
 For near-real time data, a different url source and script must be used.
 NOAA National Data Buoy Center:
  https://www.ndbc.noaa.gov/

USAGE:
 Four arguments: firstYear, lastYear, station list, output dir
 Examples (from linux/terminal command line):
  python3 get_ndbc_stdmet.py 2010 2020 allbstations.dat /home/name/Dowloads
  nohup python3 get_ndbc_stdmet.py 2010 2020 allbstations.dat /home/name/Dowloads >> nohup_get_ndbc_stdmet.out 2>&1 &

OUTPUT:
 Text files with NDBC metocean data for the buoys listed in 
  the station list (third argument).
 One file per buoy.

DEPENDENCIES:
 See dependencies.py and the imports below.

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import os
import sys
import numpy as np
import pandas as pd
import urllib.request
from pathlib2 import Path
from pylab import *

iyear=np.str(sys.argv[1]) # first year
fyear=np.str(sys.argv[2]) # final year
lstations=np.str(sys.argv[3]) # path+name of file list of stations
odir=np.str(sys.argv[4]) # outputdir

# Read list of buoys
#lndbc = open('allbstations.dat'); content = lndbc.readlines()
dfabs = pd.read_csv(lstations, comment='$'); bidabs=[]
for i in range(0,dfabs.values.shape[0]-1):
	aux=np.str(dfabs.values[i]).split()
	bidabs=np.append(bidabs,np.str(aux[3][1::]).split("'")[0])

# Dowload each buoy/station for each year
for st in range(0,size(bidabs)):

	namest = str(bidabs[st]).split()[0].lower()
	tb=0
	for year in range(np.int(iyear),np.int(fyear)+1):

		url = 'http://www.ndbc.noaa.gov/view_text_file.php?filename='+namest+'h'+repr(year)+'.txt.gz&dir=data/historical/stdmet/'

		try:
			response = urllib.request.urlopen(url) 
		except :
			print(url+"   does not exist")
		else:

			data = response.read() 
			text = data.decode('utf-8')
			del data, response, url

			if len(text.splitlines()) >10:
				text_file = open(odir+"/NDBC_Historical_stdmet_"+namest.upper()+"_"+repr(year)+".dat", "a")
				text_file.write(text)
				text_file.close()
				tb=tb+1

			del text

	if tb>0:
		os.system("cat "+odir+"/NDBC_Historical_stdmet_"+namest.upper()+"_*.dat >> "+odir+"/NDBC_historical_stdmet_"+namest.upper()+".txt")
		os.system("rm -f "+odir+"/NDBC_Historical_stdmet_"+namest.upper()+"_*.dat")

		path = Path(odir+"/NDBC_historical_stdmet_"+namest.upper()+".txt")
		text = path.read_text()
		text = text.replace("YYYY","#YYYY")
		path.write_text(text)
		del path, text

	else:
		print(" "); print(" No data for "+namest+"  from "+iyear+" to "+fyear); print(" ")

	del tb, namest

