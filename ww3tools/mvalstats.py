#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
mvalstats.py

VERSION AND LAST UPDATE:
 v1.0  04/04/2022
 v1.1  05/22/2022
 v1.2  07/13/2023

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
 See setup.py and the imports below.

AUTHOR and DATE:
 04/04/2022: Ricardo M. Campos, first version.
 05/22/2022: Ricardo M. Campos, removed L-moments from smrstat.
 07/13/2023: Ricardo M. Campos, argument reading has been modified and
   L-moments are optional in smrstat.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

from pylab import *
from matplotlib.mlab import *
import numpy as np
import scipy.stats
import sys
import warnings; warnings.filterwarnings("ignore")

# Table summary statistics
def smrstat(x, vmin=-np.inf, vmax=np.inf, lmoments='no'):
    '''
    Summary Statistics
    Input: one array of interest
    Output: mean, variance, skewness, kurtosis, min, max, percentile80, percentile90, percentile95, percentile99, percentile99.9
    In addition to the data array, users can enter minimum and maximum values, for a simple quality-control.
    If you want to obtain the L-moments (Hosking & Wallis) as well, just include lmoments='yes'
    Examples:
      import mvalstats
      sresult = mvalstats.smrstat(hs)
      sresult = mvalstats.smrstat(hs,0,20)
      sresult = mvalstats.smrstat(hs,vmax=15)
      sresult = mvalstats.smrstat(hs,vmin=0,vmax=20)
      sresult = mvalstats.smrstat(hs,vmin=0,vmax=20,lmoments='yes')
    '''
    x=np.array(x).astype('float')
    ind=np.where((np.isnan(x)==False) & (x>vmin) & (x<vmax))
    if np.any(ind):
        x=np.copy(x[ind[0]])
    else:
        raise ValueError(' Array without valid numbers.')

    if lmoments=='no':
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
    else:
        ferr=np.zeros((14),'f')*np.nan
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
        import lmoments3 as lm
        hwlm = lm.lmom_ratios(x, nmom=5)
        ferr[11] = hwlm[1]/hwlm[0]
        ferr[12] = hwlm[2]
        ferr[13] = hwlm[3]

    return ferr


# Error metrics =============================

def metrics(model,obs,vmin=-np.inf,vmax=np.inf,maxdiff=np.inf, pctlerr='no'):
    '''
    Error Metrics. Equations from Mentaschi et al. (2013)
     https://doi.org/10.1016/j.ocemod.2013.08.003
    Input: two arrays of model and observation, respectively.
        They must have the same size
    In addition to the data arrays, users can enter minimum and maximum values, 
        and maximum_difference, for a simple quality-control. 
    Additional output of 95% can be obtained if pctlerr!='no'
    Output: numpy array with shape equal to 9:
        bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N 
    Examples:
      import mvalstats
      mresult = mvalstats.metrics(model_hs,obs_hs)
      mresult = mvalstats.metrics(model_hs,obs_hs,0,20)
      mresult = mvalstats.metrics(model_hs,obs_hs,vmax=15)
      mresult = mvalstats.metrics(model,obs,vmin=0.1,vmax=20)
      mresult = mvalstats.metrics(model,obs,vmin=0.1,maxdiff=15)
      mresult = mvalstats.metrics(model,obs,vmin=0.1,vmax=20,maxdiff=15)
    '''

    model=np.atleast_1d(np.array(model)); obs=np.atleast_1d(np.array(obs))
    if model.shape != obs.shape:
        raise ValueError(' Model and Observations with different size.')
    if vmax<=vmin:
        raise ValueError(' vmin cannot be higher than vmax.')

    ind=np.where((np.isnan(model)==False) & (np.isnan(obs)==False) & (model>vmin) & (model<vmax) & (obs>vmin) & (obs<vmax) & (np.abs(model-obs)<=maxdiff) )
    if np.size(ind)>1:
        model=np.copy(model[ind[0]]); obs=np.copy(obs[ind[0]])
    else:
        raise ValueError(' Array without valid numbers.')

    del ind
    ferr=np.zeros((9),'f')*np.nan
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
    ferr[8]=len(obs)  #number of match-ups between model/obs used in stat calculations 

    # Bias and RMSE for severe conditions above the 95 percentile.
    if pctlerr!='no':
        ind=np.where((np.isnan(model)==False) & (np.isnan(obs)==False) & (model>vmin) & (model<vmax) & (obs>vmin) & (obs<vmax) & (np.abs(model-obs)<=maxdiff) & (obs>np.nanpercentile(obs,95)) )
        if np.size(ind)>1:
            ferr=np.append(ferr, model[ind[0]].mean()-obs[ind[0]].mean()) # Bias
            ferr=np.append(ferr, ((((model[ind[0]]-obs[ind[0]])**2).mean())**0.5)) # RMSE
            ferr=np.append(ferr, len(model[ind[0]])) #number of obs in 95 percentile 
        else:
            ferr=np.append(ferr, np.nan); ferr=np.append(ferr, np.nan); ferr=np.append(ferr, np.nan)

    return ferr


