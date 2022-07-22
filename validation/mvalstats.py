#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mvalstats.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022
 v1.1  05/22/2022

PURPOSE:
 Group of python functions for statistical analyses and validation. 
 Users can import as a standard python function, and use it accordingly:
 For example:
  import mvalstats
  mvalstats.metrics(model,obs,0.01,20,15)
 Users can help() each function to obtain information about inputs/outputs
  help(mvalstats.metrics)

USAGE:
 functions
   smrstat
   metrics
 Explanation of each function is contained in the headers, including
  examples.

OUTPUT:
 summary statistics values; and error metrics 

DEPENDENCIES:
 See dependencies.py and the imports below.

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.
 05/22/2022: Ricardo M. Campos, removed L-moments from smrstat.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

from pylab import *
from matplotlib.mlab import *
import numpy as np
import scipy.stats
# import lmoments3 as lm
# pip install git+https://github.com/OpenHydrology/lmoments3.git
import sys
import warnings; warnings.filterwarnings("ignore")

# Table summary statistics
def smrstat(*args):
	'''
	Summary Statistics
	Input: one array of interest
	Output: mean, variance, skewness, kurtosis, min, max, percentile80, percentile90, percentile95, percentile99, percentile99.9
	In addition to the data array, users can enter minimum and maximum values, for a simple quality-control.
	Example:
	  import mvalstats
	  sresult = mvalstats.smrstat(hs)
	  sresult = mvalstats.smrstat(hs,0,20)
	'''
	
	vmin=-np.inf; vmax=np.inf
	if len(args) == 1:
		x=np.array(args[0])
	elif len(args) == 2:
		x=np.array(args[0]); vmin=np.float(args[1])
	elif len(args) == 3:
		x=np.array(args[0]); vmin=np.float(args[1]); vmax=np.float(args[2])
	elif len(args) > 3:
		sys.exit(' Too many inputs')

	ind=np.where((np.isnan(x)==False) & (x>vmin) & (x<vmax))
	if np.any(ind):
		x=np.copy(x[ind[0]])
	else:
		sys.exit(' Array without valid numbers.')

	# ferr=np.zeros((14),'f')*np.nan
	ferr=np.zeros((11),'f')*np.nan
	ferr[0] = np.mean(x)
	ferr[1] = np.var(x)
	ferr[2] = scipy.stats.skew(x)
	ferr[3] = scipy.stats.kurtosis(x)
	ferr[4] = np.min(x)
	ferr[5] = np.max(x)
	ferr[6] = np.percentile(x,80)
	ferr[7] = np.percentile(x,90)
	ferr[8] = np.percentile(x,95)
	ferr[9] = np.percentile(x,99)
	ferr[10] = np.percentile(x,99.9)
	# Hosking & Wallis L-moment ratios
	# pip install git+https://github.com/OpenHydrology/lmoments3.git
	# hwlm = lm.lmom_ratios(x, nmom=5)
	# ferr[11] = hwlm[1]/hwlm[0]
	# ferr[12] = hwlm[2]
	# ferr[13] = hwlm[3]

	return ferr


# Error metrics =============================

def metrics(*args):
	'''
	Error Metrics. Equations from Mentaschi et al. (2013)
	 https://doi.org/10.1016/j.ocemod.2013.08.003
	Input: two arrays of model and observation, respectively.
		They must have the same size
	In addition to the data arrays, users can enter minimum and maximum values, 
		and maximum_difference, for a simple quality-control.
	Output: numpy array with shape equal to 8:
		bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC
	Example:
	  import mvalstats
	  mresult = mvalstats.metrics(model_hs,obs_hs)
	  mresult = mvalstats.metrics(model_hs,obs_hs,0,20)
	'''

	vmin=-np.inf; vmax=np.inf; maxdiff=np.inf
	if len(args) < 2:
		sys.exit(' Need two arrays with model and observations.')
	elif len(args) == 2:
		model=np.array(args[0]); obs=np.array(args[1])
	elif len(args) == 3:
		model=np.array(args[0]); obs=np.array(args[1])
		vmin=np.float(args[2])
	elif len(args) == 4:
		model=np.array(args[0]); obs=np.array(args[1])
		vmin=np.float(args[2]); vmax=np.float(args[3])
	elif len(args) == 5:
		model=np.array(args[0]); obs=np.array(args[1])
		vmin=np.float(args[2]); vmax=np.float(args[3]); maxdiff=np.float(args[4])
	elif len(args) > 5:
		sys.exit(' Too many inputs')

	model=np.atleast_1d(model); obs=np.atleast_1d(obs)
	if model.shape != obs.shape:
		sys.exit(' Model and Observations with different size.')
	if vmax<=vmin:
		sys.exit(' vmin cannot be higher than vmax.')

	ind=np.where((np.isnan(model)==False) & (np.isnan(obs)==False) & (model>vmin) & (model<vmax) & (obs>vmin) & (obs<vmax) & (np.abs(model-obs)<=maxdiff) )
	if np.any(ind) or model.shape[0]==1:
		model=np.copy(model[ind[0]]); obs=np.copy(obs[ind[0]])
	else:
		sys.exit(' Array without valid numbers.')

	ferr=np.zeros((8),'f')*np.nan
	ferr[0] = model.mean()-obs.mean() # Bias
	ferr[1] = (((model-obs)**2).mean())**0.5 # RMSE
	if obs.mean()!=0.:
		ferr[2] = ferr[0] / np.abs(obs.mean()) # Normalized Bias
	ferr[3] = ( ((model-obs)**2).sum() / (obs**2).sum() )**0.5  # Normalized RMSE
	# ferr[4] = ((((model-model.mean())-(obs-obs.mean()))**2).mean())**0.5   # Scatter Component of RMSE
	if ( (ferr[1]**2) - (ferr[0]**2) ) >= 0.:
		ferr[4] = ( (ferr[1]**2) - (ferr[0]**2) )**0.5
	ferr[5] = ( (((model-model.mean())-(obs-obs.mean()))**2).sum() / (obs**2).sum() )**0.5  # Scatter Index
	ferr[6] = ( ((model - obs)**2).sum() / (model * obs).sum() )**0.5  # HH
	ferr[7]=np.corrcoef(model,obs)[0,1]  #  Correlation Coefficient 

	return ferr


