import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from matplotlib.mlab import *
import numpy as np
import copy
import scipy.stats
import lmoments3 as lm
# pip install git+https://github.com/OpenHydrology/lmoments3.git
import sys
import warnings; warnings.filterwarnings("ignore")

# Plot Style ---
sl=14
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})


# Table summary statistics
def smrstat(*args):
	'''
	Summary Statistics
	Input: one array of interest
	Output: mean, variance, skewness, kurtosis, min, max, percentile80, percentile90, percentile95, percentile99, percentile99.9
	'''
	vmin=-np.inf; vmax=np.inf
	if len(args) == 1:
		x=copy.copy(args[0])
	elif len(args) == 2:
		x=copy.copy(args[0]); vmin=copy.copy(args[1])
	elif len(args) == 3:
		x=copy.copy(args[0]); vmin=copy.copy(args[1]); vmax=copy.copy(args[2])
	elif len(args) > 3:
		sys.exit(' Too many inputs')

	ind=np.where((np.isnan(x)==False) & (x>vmin) & (x<vmax))
	if np.any(ind):
		x=np.copy(x[ind[0]])
	else:
		sys.exit(' Array without valid numbers.')

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
	hwlm = lm.lmom_ratios(x, nmom=5)
	ferr[11] = hwlm[1]/hwlm[0]
	ferr[12] = hwlm[2]
	ferr[13] = hwlm[3]

	return ferr


# Error metrics =============================

def metrics(*args):
	'''
	Error Metrics. Mentaschi et al. (2013)
	Input: two arrays of model and observation, respectively.
		They must have the same size
	Output: ferr array with shape equal to 8
		bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC
	'''

	vmin=-np.inf; vmax=np.inf; maxdiff=np.inf
	if len(args) < 2:
		sys.exit(' Need two arrays with model and observations.')
	elif len(args) == 2:
		model=copy.copy(args[0]); obs=copy.copy(args[1])
	elif len(args) == 3:
		model=copy.copy(args[0]); obs=copy.copy(args[1])
		vmin=copy.copy(args[2])
	elif len(args) == 4:
		model=copy.copy(args[0]); obs=copy.copy(args[1])
		vmin=copy.copy(args[2]); vmax=copy.copy(args[3])
	elif len(args) == 5:
		model=copy.copy(args[0]); obs=copy.copy(args[1])
		vmin=copy.copy(args[2]); vmax=copy.copy(args[3]); maxdiff=copy.copy(args[4]);
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


