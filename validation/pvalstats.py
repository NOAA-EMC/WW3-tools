#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pvalstats.py

VERSION AND LAST UPDATE:
 v1.0  05/24/2022
 v1.1  06/30/2022
 v1.2  11/23/2022

PURPOSE:
 Group of python functions for visualization of model/observation data
  and validation and assessment results.
 Users can import as a standard python function, and use it accordingly:
 For example:
  import pvalstats
  pvalstats.qqplot(model,obs)
 Users can help() each function to obtain information about inputs/outputs
  help(pvalstats.scatterplot)

USAGE:
 functions
   timeseries
   qqplot
   scatterplot
   taylordiagram
   combinerrors
   pdf

 Explanation of each function is contained in the headers, including
  examples.

OUTPUT:
 png figures saved in the local directory or in the path given through tag/prefix.

DEPENDENCIES:
 See dependencies.py and the imports below.

AUTHOR and DATE:
 05/24/2022: Ricardo M. Campos, first version.
 06/30/2022: Ricardo M. Campos, TaylorDiagram and PDF plots corrected
   for arrays containing NaN.
 11/23/2022: Ricardo M. Campos, new functions errXftime and interp_nan,
   plus some format improvements.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
import math
import pandas as pd
from matplotlib.mlab import *
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import gaussian_kde, linregress
from scipy import interpolate
import matplotlib.ticker
from matplotlib.dates import DateFormatter
from mpl_toolkits.basemap import cm
colormap = cm.GMT_polar
palette = plt.cm.jet
palette.set_bad('aqua', 10.0)
import sys
import warnings; warnings.filterwarnings("ignore")
# Font size and style
sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

def interp_nan(*args):
	'''
	Fill NaN values with linear interpolation.
	User can enter one or two inputs:
	  1) time-series containing NaN values to fill
	  2) maximum number of consecutive NaN values to interpolate (to
	     avoid interpolating long segments)
	'''

	data=np.array(args[0])
	if len(args)==1:
		lmt=10**50
	elif len(args)>1:
		lmt=np.int(args[1])

	if data.ndim>1:
		sys.exit(' Input array with too many dimensions. Only time-series (1 dimension) allowed.')
	else:
		# using pandas
		A=pd.Series(data)
		# B=A.interpolate(method="polynomial",order=2,limit=lmt)
		B=A.interpolate(method="linear",limit=lmt)
		B=np.array(B.values)

	return B
	del lmt,A,data


# Time-Series plot
def timeseries(*args):
	'''
	Simple time-series plots of model (lines) and observations (discrete points).
	The x-axis is time and y-axis is the variable and unit of the array provided.
	Inputs:
	  Mandatory: model(s); observations; time array (datetime).
	   the model array can include one or more model results, through the number of columns, while
	   the observation array must be one-dimensional, with the same number of lines as the model.
	  Optional: y-axis name (string with name and/or unit of variable plotted)
	            tag (string with a name or path/name to save the figure)
	            labels (string array with one or mode model labels that allow recognizing the model in the legend)
	            color array
	            linestyle stype array

	Output: png figure saved in the local directory where python is running or in the path given through ftag.
	Example:
	  import pvalstats
	  pvalstats.timeseries(model,buoydata,ftime)
	  pvalstats.timeseries(np.c_[model1,model2].T,buoydata,ftime,'Hs (m)','/home/rmc/test/WW3ST4Bmax_',['WW3B133','WW3B137'])
	'''

	yaxname=[]; ftag=''; mlabels=[]
	mmark=np.array(np.atleast_1d(['-','--','-.','--','-.','--','-.','--','-.','--','-.','--','-.','--','-.'])).astype('str')
	ccol=np.array(np.atleast_1d(['navy','firebrick','darkgreen','fuchsia','gold',
		'blue','salmon','lime','darkviolet','yellow',
		'cornflowerblue','red','green','violet','orange'])).astype('str')

	if len(args) < 3:
		sys.exit(' At least three arrays must be given: model, observation, and time.')
	if len(args) >= 3:
		model=np.array(np.atleast_2d(args[0])).astype('float')
		obs=np.array(np.atleast_2d(args[1])).astype('float')
		ftime=np.array(args[2])
	if len(args) >= 4:
		yaxname=np.array(args[3]).astype('str')
	if len(args) >= 5:
		ftag=np.str(args[4])
	if len(args) >= 6:
		mlabels=np.array(np.atleast_1d(args[5])).astype('str')
	if len(args) >= 7:
		ccol=np.array(np.atleast_1d(args[6])).astype('str')
	if len(args) >= 8:
		mmark=np.array(np.atleast_1d(args[7])).astype('str')	
	if len(args) > 8:
		sys.exit(' Too many inputs')

	if (size(model.shape)>2) | (size(obs.shape)>2):
		sys.exit(' Too many dimensions in the model or observation arrays.')

	if model.shape[0]>model.shape[1]:
		model=np.array(model).T
	if obs.shape[0]>obs.shape[1]:
		obs=np.array(obs).T
	if (model.shape[1] != obs.shape[1]):
		sys.exit(' Model and observation arrays must have the same index/time size.')

	# plot
	# plt.close('all')
	fig1 = plt.figure(1,figsize=(9,4)); ax = fig1.add_subplot(111)
	if size(mlabels)>0 and size(np.where(obs[0,:]>-999))>1:
		ax.plot(ftime,obs[0,:],color='dimgray',marker='.',linestyle='',linewidth=2.,label='Obs',zorder=2)
	else:
		ax.plot(ftime,obs[0,:],color='dimgray',marker='.',linestyle='',linewidth=2.,zorder=2)

	if size(np.where(obs[0,:]>-999))>1:
		ax.fill_between(ftime, 0., obs[0,:], color='silver',alpha=0.5,zorder=1)

	for i in range(0,model.shape[0]):
		if size(mlabels)>0:
			if size(np.where(model[i,:]>-999))>1:
				ax.plot(ftime,interp_nan(model[i,:],10), color=ccol[i],linestyle=mmark[i],linewidth=2.,label=mlabels[i],alpha=0.8,zorder=3)
		else:
			if size(np.where(model[i,:]>-999))>1:
				ax.plot(ftime,interp_nan(model[i,:],10), color=ccol[i],linestyle=mmark[i],linewidth=2.,alpha=0.8,zorder=3)

		if size(np.where(model[i,:]>-999))>1:
			ax.plot(ftime,model[i,:],marker='.',markersize=5,markerfacecolor=ccol[i],markeredgewidth=0.5,markeredgecolor='black', linestyle='',linewidth=0.,zorder=3)

	ax.set_xlim(ftime[0],ftime[-1])
	ax.xaxis.set_major_formatter( DateFormatter('%b%d') )
	ax.fmt_xdata = DateFormatter('%b%d')
	# plt.ylabel(wvars[i]+' ('+ds[wvars[i]].units+')', fontsize=9)
	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(sl-2)

	ax.set_xlabel("Time")
	if size(yaxname)>0:
		ax.set_ylabel(yaxname);
	else:
		ax.set_ylabel("Model and Observations")

	if size(mlabels)>0:
		if size(mlabels)<7:
			plt.legend(fontsize=sl-3)
		else:
			plt.legend(fontsize=sl-5)

	plt.tight_layout(); # plt.axis('tight') 
	plt.grid(c='grey', ls='--', alpha=0.3,zorder=1)
	plt.savefig(ftag+'TimeSeries.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig1); del fig1, ax

# QQ-plot
def qqplot(*args):
	'''
	Quantile-Quantile plot.
	Inputs:
	  Mandatory: model(s); observations (or any reference).
	   the model array can include one or more model results, through the number of columns, while
	   the observation array must be one-dimensional, with the same number of lines as the model.
	  Optional: tag (string with a name or path/name to save the figure)
	            labels (string array with one or mode model labels that allow recognizing the model in the legend)
	            color array
	            marker stype array

	Output: png figure saved in the local directory where python is running or in the path given through ftag.
	Example:
	  import pvalstats
	  pvalstats.qqplot(model,satdata)
	  pvalstats.qqplot(np.c_[model1,model2].T,satdata,'/home/rmc/test/WW3PR_',['WW3PR3','WW3PR1'])
	'''

	ftag=''; mlabels=[]
	mmark=np.array(np.atleast_1d(['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',])).astype('str')
	ccol=np.array(np.atleast_1d(['navy','firebrick','darkgreen','fuchsia','gold',
		'blue','salmon','lime','darkviolet','yellow',
		'cornflowerblue','red','green','violet','orange'])).astype('str')

	if len(args) < 2:
		sys.exit(' At least two arrays must be given: model and observation.')
	if len(args) >= 2:
		model=np.array(np.atleast_2d(args[0])).astype('float'); obs=np.array(np.atleast_2d(args[1])).astype('float')
	if len(args) >= 3:
		ftag=np.str(args[2])
	if len(args) >= 4:
		mlabels=np.array(np.atleast_1d(args[3])).astype('str')
	if len(args) >= 5:
		ccol=np.array(np.atleast_1d(args[4])).astype('str')
	if len(args) >= 6:
		mmark=np.array(np.atleast_1d(args[5])).astype('str')
	if len(args) > 6:
		sys.exit(' Too many inputs')

	if (size(model.shape)>2) | (size(obs.shape)>2):
		sys.exit(' Too many dimensions in the model or observation arrays.')

	if model.shape[0]>model.shape[1]:
		model=np.array(model).T
	if obs.shape[0]>obs.shape[1]:
		obs=np.array(obs).T
	if (model.shape[1] != obs.shape[1]):
		sys.exit(' Model and observation arrays must have the same index/time size.')

	# percentiles
	p = np.arange(1,96,1); p=np.append(p,[95.5,96.,96.5,97.,97.3,97.6,97.9,98.2,98.5,98.7,98.9,99.1,99.3,99.5,99.7])
	qm = np.zeros((model.shape[0],p.shape[0]),'f')*np.nan
	qobs = np.zeros((model.shape[0],p.shape[0]),'f')*np.nan
	for i in range(0,p.shape[0]):
		for j in range(0,model.shape[0]):
			qobs[j,i] = np.nanpercentile(obs[:],p[i])
			qm[j,i] = np.nanpercentile(model[j,:],p[i])

	a=math.floor(np.nanmin([qobs,qm])*100.)/100. ; b=math.ceil(np.nanmax([qobs,qm])*100.)/100.
	aux=np.linspace(a-0.2*a,b+0.2*a,p.shape[0])
	# plot
	# plt.close('all')
	fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)
	ax.plot(aux,aux,'k', linewidth=2.,alpha=0.4,zorder=2)  # main diagonal
	for i in range(0,model.shape[0]):
		if size(mlabels)>0:
			ax.plot(qobs[i,:],qm[i,:], color=ccol[i], marker=mmark[i], linestyle='--',label=mlabels[i],linewidth=1.,alpha=0.8,zorder=3)
		else:
			ax.plot(qobs[i,:],qm[i,:], color=ccol[i], marker=mmark[i], linestyle='--',linewidth=1.,alpha=0.8,zorder=3)

	plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
	for i in np.array([50,80,90,95,99]):
		plt.axvline(x= np.nanpercentile(obs,np.int(i)),ls='--',color='grey',linewidth=1.,alpha=0.9,zorder=1)
		plt.text(np.nanpercentile(obs,np.int(i)),(aux.max()-aux.min())/15 + aux.min(),np.str(np.int(i))+'th',color='dimgrey',fontsize=sl-7)

	plt.ylim(ymax = aux.max(), ymin = aux.min())
	plt.xlim(xmax = aux.max(), xmin = aux.min())
	plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7) 
	ax.set_xlabel("Observations"); ax.set_ylabel("Model")
	if size(mlabels)>0:
		plt.legend(loc="upper left",fontsize=sl-2)

	plt.tight_layout(); #plt.axis('tight') 
	plt.savefig(ftag+'QQplot.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig1); del fig1, ax

# Scatter Plot
def scatterplot(*args):
	'''
	Scatter plot.
	Inputs:
	  Mandatory: model(s); observations (or any reference).
	   the model array can include one or more model results, through the number of columns, while
	   the observation array must be one-dimensional, with the same number of lines as the model.
	  Optional: tag (string with a name or path/name to save the figure)
	            labels (string array with one or mode model labels that allow recognizing the model in the legend)
	            color array

	Output: png figure saved in the local directory where python is running or in the path given through ftag.
	Example:
	  import pvalstats
	  pvalstats.scatterplot(model,satdata)
	  pvalstats.scatterplot(np.c_[model1,model2].T,satdata,'/home/rmc/test/WW3PR_',['WW3PR3','WW3PR1'])
	'''

	ftag=''; mlabels=[]
	mmark=np.array(np.atleast_1d(['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',])).astype('str')
	ccol=np.array(np.atleast_1d(['navy','firebrick','darkgreen','fuchsia','gold',
		'blue','salmon','lime','darkviolet','yellow',
		'cornflowerblue','red','green','violet','orange'])).astype('str')

	if len(args) < 2:
		sys.exit(' At least two arrays must be given: model and observation.')
	if len(args) >= 2:
		model=np.array(np.atleast_2d(args[0])).astype('float'); obs=np.array(np.atleast_2d(args[1])).astype('float')
	if len(args) >= 3:
		ftag=np.str(args[2])
	if len(args) >= 4:
		mlabels=np.array(np.atleast_1d(args[3])).astype('str')
	if len(args) >= 5:
		ccol=np.array(np.atleast_1d(args[4])).astype('str')
	if len(args) > 5:
		sys.exit(' Too many inputs')

	if (size(model.shape)>2) | (size(obs.shape)>2):
		sys.exit(' Too many dimensions in the model or observation arrays.')

	if model.shape[0]>model.shape[1]:
		model=np.array(model).T
	if obs.shape[0]>obs.shape[1]:
		obs=np.array(obs).T
	if obs.shape[0]>1:
		sys.exit(' You can enter multiple model results for comparison but one array of observations must be used, to be consistent.')

	if (model.shape[1] != obs.shape[1]):
		sys.exit(' Model and observation arrays must have the same index/time size.')

	if (model.shape[0]>1) & (size(mlabels)==0):
		print(" Warning. You should use labels to differentiate model results.")

	# plt.close('all')
	fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)
	a=math.floor(np.nanmin(np.append(obs,model))*100.)/100. ; b=math.ceil(np.nanmax(np.append(obs,model))*100.)/100.
	famin=a-0.1*a; famax=b+0.1*a
	faux=np.linspace(famin,famax,100); del a,b
	for i in range(0,model.shape[0]):
		b=np.array(obs[0,:]); a=np.array(model[i,:])
		ind=np.where((a*b)>-999.)[0]; a=np.copy(a[ind]); b=np.copy(b[ind]); del ind
		if a.shape[0]>50000:
			sk=np.int(np.round(np.float(a.shape[0])/30000.,0))
			a=np.copy(a[::sk]); b=np.copy(b[::sk])

		#amax=math.ceil(np.nanmax(np.array([a,b]))*100.)/100.
		#amin=math.floor(np.nanmin(np.array([a,b]))*100.)/100.
		#aux=np.linspace(amin,amax,100)	

		if (a.shape[0]<30) | (model.shape[0]>1):
			if size(mlabels)>0:
				ax.scatter(b,a,color=ccol[i],marker=mmark[i],label=mlabels[i],zorder=2) #, linestyle='--', linewidth=1.
			else:
				ax.scatter(b,a,color=ccol[i],marker=mmark[i],zorder=2)
		else:
			xy = np.vstack([a,b]); z = gaussian_kde(xy)(xy)
			if size(mlabels)>0:
				ax.scatter(b,a, c=z, s=5,cmap=plt.cm.jet,label=mlabels[i],zorder=2)
			else:
				ax.scatter(b,a, c=z, s=5,cmap=plt.cm.jet,zorder=2)

		# famax=np.append(famax,amax); famin=np.append(famin,amin)
		del a,b
	
	ax.plot(faux,faux,'k--', linewidth=1.,alpha=0.9,zorder=3)  # main diagonal
	plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7)
	ax.set_xlabel("Observations"); ax.set_ylabel("Model")
	plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
	for i in np.array([50,80,90,95,99]):
		plt.axvline(x= np.nanpercentile(obs,np.int(i)),ls='--',color='grey',linewidth=1.,alpha=0.7,zorder=1)
		plt.text(np.nanpercentile(obs,np.int(i)),((famax-famin)/15)+famin,np.str(np.int(i))+'th',color='grey',fontsize=sl-7,zorder=4)

	plt.gca().set_xlim(left=famin, right=famax); plt.gca().set_ylim(ymin=famin,ymax=famax)

	if size(mlabels)>0:
		plt.legend(loc="upper left",fontsize=sl-2)

	# ax.set_axisbelow(True)
	plt.tight_layout(); #plt.axis('tight')
	plt.savefig(ftag+'ScatterPlot.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig1); del fig1, ax

# Taylor Diagram
def taylordiagram(*args):
	'''
	Taylor Diagram using skill_metrics.
	https://github.com/PeterRochford/SkillMetrics
	Inputs:
	  Mandatory: model(s); observations (or any reference).
	   the model array can include one or more model results, through the number of columns, while
	   the observation array must be one-dimensional, with the same number of lines as the model.
	  Optional: tag (string with a name or path/name to save the figure)
	            labels (string array with one or mode model labels that allow recognizing the model in the legend)
	            color array

	Output: png figure saved in the local directory where python is running or in the path given through ftag.
	Example:
	  import pvalstats
	  pvalstats.taylordiagram(model,satdata)
	  pvalstats.taylordiagram(np.c_[model1,model2].T,satdata,'/home/rmc/test/WW3PR_',['WW3PR3','WW3PR1'])
	'''

	import skill_metrics as sm
	ftag=''; mlabels=[]
	ccol=np.array(np.atleast_1d(['navy','firebrick','darkgreen','fuchsia','gold',
		'blue','salmon','lime','darkviolet','yellow',
		'cornflowerblue','red','green','violet','orange'])).astype('str')

	if len(args) < 2:
		sys.exit(' At least two arrays must be given: model and observation.')
	if len(args) >= 2:
		model=np.array(np.atleast_2d(args[0])).astype('float'); obs=np.array(np.atleast_2d(args[1])).astype('float')
	if len(args) >= 3:
		ftag=np.str(args[2])
	if len(args) >= 4:
		mlabels=np.array(np.atleast_1d(args[3])).astype('str')
	if len(args) >= 5:
		ccol=np.array(np.atleast_1d(args[4])).astype('str')
	if len(args) > 5:
		sys.exit(' Too many inputs')

	if (size(model.shape)>2) | (size(obs.shape)>2):
		sys.exit(' Too many dimensions in the model or observation arrays.')

	if model.shape[0]>model.shape[1]:
		model=np.array(model).T
	if obs.shape[0]>obs.shape[1]:
		obs=np.array(obs).T
	if obs.shape[0]>1:
		sys.exit(' You can enter multiple model results for comparison but one array of observations must be used, to be consistent.')

	if (model.shape[1] != obs.shape[1]):
		sys.exit(' Model and observation arrays must have the same index/time size.')

	if (model.shape[0]>1) & (size(mlabels)==0):
		print(" Warning. You should use labels to differentiate model results.")

	# Organizing labels for the Taylor Diagram
	if size(mlabels)>0:
		label = np.append(['Obs'],mlabels)
	else:
		if model.shape[0] > 1:
			label = np.array(['Obs'])
			for j in range(0,model.shape[0]):
				label = np.append(label,['Model'+repr(j+1)])
		else:
			label = np.array(['Obs','Model'])

	# remove NaN
	ind=np.where( (np.mean(model,axis=0)>-999.) & (np.mean(obs,axis=0)>-999.) )
	if size(ind)>0:
		model=np.copy(model[:,ind[0]]); obs=np.copy(obs[:,ind[0]]); del ind
	else:
		sys.exit(' No quality data available.')

	# plot
	# plt.close('all')
	fig1 = plt.figure(1,figsize=(7,6)); ax = fig1.add_subplot(111)
	# initial statistics
	ts = sm.taylor_statistics(model[0,:],obs[0,:])
	sdev=np.array(ts['sdev'][0]); crmsd=np.array(ts['crmsd'][0]); ccoef=np.array(ts['ccoef'][0])
	for i in range(0,model.shape[0]):
		ts = sm.taylor_statistics(model[i,:],obs[0,:])
		sdev=np.append(sdev,ts['sdev'][1])
		crmsd=np.append(crmsd,ts['crmsd'][1])
		ccoef=np.append(ccoef,ts['ccoef'][1])

	# formatting
	if crmsd.max()>0.8: npr=np.int(2)
	else: npr=np.int(3)
	tRMS = np.array(np.linspace(0,np.round(crmsd.max()+0.1*crmsd.max(),npr),5)).round(npr)

	if sdev.max()>0.8: npr=np.int(2)
	else: npr=np.int(3)
	axmax = np.round(sdev.max()+0.5*sdev.max(),npr)

	if (sdev.max()+0.4*sdev.max())>0.8: npr=np.int(2)
	else: npr=np.int(3)
	tSTD = np.array(np.linspace(sdev.min()-0.3*sdev.max(),np.round(sdev.max()+0.4*sdev.max(),npr),4)).round(npr)
	del npr

	fp1=1.07; fp2=1.03; fp3=0.055; fp4=13
	if size(label)>8:
		fp2=1.07; fp3=0.040; fp4=12
	elif size(label)>5:
		fp2=1.06; fp3=0.045

	slmax=0
	for i in range(0,size(label)):
		if len(label[i])>slmax:
			slmax=np.int(len(label[i]))

	if slmax>11:
		fp1=1.04; fp4=12
	elif slmax>8:
		fp1=1.04

	# loop through models
	for i in range(0,model.shape[0]):
		if i==0:
			sm.taylor_diagram(sdev.take([0,i+1]), crmsd.take([0,i+1]), ccoef.take([0,i+1]), locationColorBar = 'EastOutside',  markerColor = ccol[i],
				# markerLegend = 'on', markerLabel = list(label.take([0,i+1])),
				styleOBS = '-', colOBS = 'k', markerobs = 'o', titleOBS = 'Obs', 
				markerSize = 17, tickRMS = tRMS, axismax = axmax, alpha=0.7, tickSTD = tSTD,
				tickRMSangle = 87, showlabelsRMS = 'on',
				colRMS='dimgrey', styleRMS=':', widthRMS=1.0, titleRMS='off',
				colSTD='k', styleSTD='-.', widthSTD=1.0, titleSTD ='on',
				colCOR='k', styleCOR='--', widthCOR=1.0, titleCOR='on')

		else:
			sm.taylor_diagram(sdev.take([0,i+1]), crmsd.take([0,i+1]), ccoef.take([0,i+1]), overlay = 'on', locationColorBar = 'EastOutside',  markerColor = ccol[i],
				# markerLegend = 'on', markerLabel = list(label.take([0,i+1])), 
				styleOBS = '-', colOBS = 'k', markerobs = 'o', titleOBS = 'Obs', 
				markerSize = 17, tickRMS = tRMS, axismax = axmax, alpha=0.5, tickSTD = tSTD,
				tickRMSangle = 87, showlabelsRMS = 'on',
				colRMS='dimgrey', styleRMS=':', widthRMS=1.0, titleRMS='off',
				colSTD='k', styleSTD='-.', widthSTD=1.0, titleSTD ='on',
				colCOR='k', styleCOR='--', widthCOR=1.0, titleCOR='on')

		if size(mlabels)>0: 
			ax.text(fp1, fp2-(fp3*i), label[i+1], verticalalignment='bottom', horizontalalignment='right',
				transform=ax.transAxes,color=ccol[i], fontsize=fp4)

	plt.tight_layout()
	plt.savefig(ftag+'TaylorDiagram.png', dpi=250, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig1); del fig1, ax

# Combined Errors (Normalized bias, Scatter Index, Correlation Coefficient)
def combinerrors(*args):
	'''
	Combined Errors, Scatter Index and Normalized Bias
	Inputs:
	  Mandatory: model(s); observations (or any reference).
	   the model array can include one or more model results, through the number of columns, while
	   the observation array must be one-dimensional, with the same number of lines as the model.
	  Optional: tag (string with a name or path/name to save the figure)
	            labels (string array with one or mode model labels that allow recognizing the model in the legend)
	            color array

	Output: png figure saved in the local directory where python is running or in the path given through ftag.
	Example:
	  import pvalstats
	  pvalstats.combinerrors(model,satdata)
	  pvalstats.combinerrors(np.c_[model1,model2].T,satdata,'/home/rmc/test/WW3PR_',['WW3PR3','WW3PR1'])
	'''

	import mvalstats
	ftag=''; mlabels=[]
	ccol=np.array(np.atleast_1d(['navy','firebrick','darkgreen','fuchsia','gold',
		'blue','salmon','lime','darkviolet','yellow',
		'cornflowerblue','red','green','violet','orange'])).astype('str')

	if len(args) < 2:
		sys.exit(' At least two arrays must be given: model and observation.')
	if len(args) >= 2:
		model=np.array(np.atleast_2d(args[0])).astype('float'); obs=np.array(np.atleast_2d(args[1])).astype('float')
	if len(args) >= 3:
		ftag=np.str(args[2])
	if len(args) >= 4:
		mlabels=np.array(np.atleast_1d(args[3])).astype('str')
	if len(args) >= 5:
		ccol=np.array(np.atleast_1d(args[4])).astype('str')
	if len(args) > 5:
		sys.exit(' Too many inputs')

	if (size(model.shape)>2) | (size(obs.shape)>2):
		sys.exit(' Too many dimensions in the model or observation arrays.')

	if model.shape[0]>model.shape[1]:
		model=np.array(model).T
	if obs.shape[0]>obs.shape[1]:
		obs=np.array(obs).T
	if obs.shape[0]>1:
		sys.exit(' You can enter multiple model results for comparison but one array of observations must be used, to be consistent.')

	if (model.shape[1] != obs.shape[1]):
		sys.exit(' Model and observation arrays must have the same index/time size.')

	if (model.shape[0]>1) & (size(mlabels)==0):
		print(" Warning. You should use labels to differentiate model results.")

	# plot
	# plt.close('all')
	fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)
	fauxnb=[];fauxsi=[];
	for i in range(0,model.shape[0]):
		aux = mvalstats.metrics(model[i,:],obs[0,:])
		if size(mlabels)>0:
			ax.scatter(aux[2],aux[5], s=40+150*aux[7]**5,color=ccol[i],edgecolors='k',linewidth=1,alpha=0.7,label=mlabels[i],zorder=3)
		else:
			ax.scatter(aux[2],aux[5], s=40+150*aux[7]**5,color=ccol[i],edgecolors='k',linewidth=1,alpha=0.7,zorder=3)

		fauxnb=np.append(fauxnb,np.abs(aux[2]))
		fauxsi=np.append(fauxsi,np.abs(aux[5]))
		del aux

	plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
	plt.axvline(x=0,ls='--',color='dimgrey',linewidth=1.,alpha=0.9,zorder=2)
	plt.locator_params(axis='y', nbins=5) ; plt.locator_params(axis='x', nbins=5) 
	ax.set_xlabel("Normalized Bias"); ax.set_ylabel("Scatter Index")
	if (np.nanmax(fauxnb)-np.nanmin(fauxnb))<0.001:
		plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
	if (np.nanmax(fauxsi)-np.nanmin(fauxsi))<0.001:
		plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
	if size(mlabels)>0:
		plt.legend(fontsize=sl-2)

	plt.tight_layout(); #plt.axis('tight') 
	plt.savefig(ftag+'NBiasSI.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig1); del fig1, ax

# Probability Density Function plots
def pdf(*args):
	'''
	Probability Density Plot. Observations (shaded grey) and model(s)
	Inputs:
	  Mandatory: model(s); observations (or any reference).
	   the model array can include one or more model results, through the number of columns, while
	   the observation array must be one-dimensional, with the same number of lines as the model.
	  Optional: x-axis name (string with name and/or unit of variable plotted)
	            tag (string with a name or path/name to save the figure)
	            labels (string array with one or mode model labels that allow recognizing the model in the legend)
	            color array

	Output: two png figure saved in the local directory where python is running or in the path given through ftag.
	Example:
	  import pvalstats
	  pvalstats.pdf(model,satdata)
	  pvalstats.pdf(np.c_[model1,model2].T,satdata,'Hs (m)','/home/rmc/test/WW3PR_',['WW3PR3','WW3PR1'])
	'''

	from scipy.stats import gaussian_kde
	xaxname=[]; ftag=''; mlabels=[]
	mmark=np.array(np.atleast_1d(['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',])).astype('str')
	ccol=np.array(np.atleast_1d(['navy','firebrick','darkgreen','fuchsia','gold',
		'blue','salmon','lime','darkviolet','yellow',
		'cornflowerblue','red','green','violet','orange'])).astype('str')

	if len(args) < 2:
		sys.exit(' At least two arrays must be given: model and observation.')
	if len(args) >= 2:
		model=np.array(np.atleast_2d(args[0])).astype('float'); obs=np.array(np.atleast_2d(args[1])).astype('float')
	if len(args) >= 3:
		xaxname=np.array(args[2]).astype('str')
	if len(args) >= 4:
		ftag=np.str(args[3])
	if len(args) >= 5:
		mlabels=np.array(np.atleast_1d(args[4])).astype('str')
	if len(args) >= 6:
		ccol=np.array(np.atleast_1d(args[5])).astype('str')
	if len(args) > 7:
		sys.exit(' Too many inputs')

	if (size(model.shape)>2) | (size(obs.shape)>2):
		sys.exit(' Too many dimensions in the model or observation arrays.')

	if model.shape[0]>model.shape[1]:
		model=np.array(model).T
	if obs.shape[0]>obs.shape[1]:
		obs=np.array(obs).T
	if obs.shape[0]>1:
		sys.exit(' You can enter multiple model results for comparison but one array of observations must be used, to be consistent.')

	if (model.shape[1] != obs.shape[1]):
		sys.exit(' Model and observation arrays must have the same index/time size.')

	if (model.shape[0]>1) & (size(mlabels)==0):
		print(" Warning. You should use labels to differentiate model results.")

	# remove NaN
	ind=np.where( (np.mean(model,axis=0)>-999.) & (np.mean(obs,axis=0)>-999.) )
	if size(ind)>0:
		model=np.copy(model[:,ind[0]]); obs=np.copy(obs[:,ind[0]]); del ind
	else:
		sys.exit(' No quality data available.')

	# plot
	px=np.linspace(np.nanmin([np.nanmin(obs),np.nanmin(model)]),np.nanmax([np.nanmax(obs),np.nanmax(model)])*0.99,100)

	# plt.close('all')
	fig1 = plt.figure(1,figsize=(5,5)); ax = fig1.add_subplot(111)
	dx = gaussian_kde(obs[0,:])
	ax.fill_between(px, 0., gaussian_filter(dx(px), 1.), color='grey', alpha=0.7,zorder=1)
	ax.plot(px,gaussian_filter(dx(px), 1.),color='dimgrey',linestyle='-',linewidth=0.5,alpha=0.7,zorder=2)
	for i in range(0,model.shape[0]):
		de = gaussian_kde(model[i,:])
		if size(mlabels)>0:
			ax.plot(px,gaussian_filter(de(px), 1.),color=ccol[i],marker=mmark[i],linestyle='-',linewidth=1.,alpha=0.7,label=mlabels[i],zorder=3)
		else:
			ax.plot(px,gaussian_filter(de(px), 1.),color=ccol[i],marker=mmark[i],linestyle='-',linewidth=1.,alpha=0.7,zorder=3)

		del de

	plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
	plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7) 
	ax.set_ylabel("Probability Density")

	if size(xaxname)>0:
		ax.set_xlabel(xaxname);
	else:
		ax.set_xlabel("Model and Observations")

	if size(mlabels)>0:
		plt.legend(fontsize=sl-2)

	plt.tight_layout(); plt.axis('tight'); plt.ylim(ymin = -0.005) 
	plt.savefig(ftag+'PDF.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	ax.set_yscale( "log" ); plt.axis('tight');plt.ylim(ymin = 0.001);plt.ylim(ymax = 1.)
	plt.savefig(ftag+'PDFlogscale.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig1); del fig1, ax


# Plot error versus forecast lead time
def errXftime(*args):
	'''
	Plot of previously calculated errors (y-axis) versus forecast lead time (x-axis).
	Inputs:
	  Mandatory: model(s) and forecast time array.
	   the model array can include one or more model results, through the number of columns, while
	   the forecast lead time array must be one-dimensional (any dimension, hours, days, weeks).
	  Optional: y-axis name (string with name and/or unit of error metric plotted)
	            x-axis name (string with time unit name, hours, days, weeks etc)
	            tag (string with a name or path/name to save the figure)
	            labels (string array with one or mode model labels that allow recognizing the model in the legend)
	            color array
	            linestyle stype array

	Output: png figure saved in the local directory where python is running or in the path given through ftag.
	Example:
	  import pvalstats
	  pvalstats.errXftime(modelerr,ftime)
	  pvalstats.errXftime(np.c_[modelerr1,modelerr2].T,ftime,'RMSE (m)','Lead Time (hours)','/home/rmc/test/WW3ST4Bmax_',['WW3B133','WW3B137'])
	'''

	yaxname=[]; xaxname=[]; ftag=''; mlabels=[]
	mmark=np.array(np.atleast_1d(['-','--','-.','--','-.','--','-.','--','-.','--','-.','--','-.','--','-.'])).astype('str')
	ccol=np.array(np.atleast_1d(['navy','firebrick','darkgreen','fuchsia','gold',
		'blue','salmon','lime','darkviolet','yellow',
		'cornflowerblue','red','green','violet','orange'])).astype('str')

	if len(args) < 2:
		sys.exit(' At least three arrays must be given: model, observation, and time.')
	if len(args) >= 2:
		model=np.array(np.atleast_2d(args[0])).astype('float')
		ftime=np.array(args[1])
	if len(args) >= 3:
		yaxname=np.array(args[2]).astype('str')
	if len(args) >= 4:
		xaxname=np.array(args[3]).astype('str')
	if len(args) >= 5:
		ftag=np.str(args[4])
	if len(args) >= 6:
		mlabels=np.array(np.atleast_1d(args[5])).astype('str')
	if len(args) >= 7:
		ccol=np.array(np.atleast_1d(args[6])).astype('str')
	if len(args) >= 8:
		mmark=np.array(np.atleast_1d(args[7])).astype('str')	
	if len(args) > 8:
		sys.exit(' Too many inputs')

	if (size(model.shape)>2) | (size(ftime.shape)>1):
		sys.exit(' Too many dimensions in the model or ftime arrays.')

	#if model.shape[0]>model.shape[1]:
	#	model=np.array(model).T
	#if (model.shape[0] != ftime.shape[0]):
	#	sys.exit(' Model and ftime arrays must have the same time size.')

	# plot
	# plt.close('all')
	fig1 = plt.figure(1,figsize=(8,5)); ax = fig1.add_subplot(111)
	for i in range(0,model.shape[0]):
		if size(mlabels)>0:
			# ax.plot(ftime,gaussian_filter(model[i,:],1), color=ccol[i],linestyle=mmark[i],linewidth=2.,label=mlabels[i],alpha=0.8,zorder=3)
			ax.plot(ftime,model[i,:], color=ccol[i],linestyle=mmark[i],linewidth=2.,label=mlabels[i],alpha=0.8,zorder=3)
		else:
			# ax.plot(ftime,gaussian_filter(model[i,:],1), color=ccol[i],linestyle=mmark[i],linewidth=2.,alpha=0.8,zorder=3)
			ax.plot(ftime,model[i,:], color=ccol[i],linestyle=mmark[i],linewidth=2.,alpha=0.8,zorder=3)

	ax.set_xlim(ftime[0],ftime[-1])
	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(sl-2)

	if size(xaxname)>0:
		ax.set_xlabel(xaxname)
	else:
		ax.set_xlabel("Forecast Lead Time")

	if size(yaxname)>0:
		ax.set_ylabel(yaxname)
	else:
		ax.set_ylabel("Model Error")

	if size(mlabels)>0:
		plt.legend(fontsize=sl-3)

	plt.tight_layout(); # plt.axis('tight') 
	plt.grid(c='grey', ls='--', alpha=0.3,zorder=1)
	plt.savefig(ftag+'ErrXFTime.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
	plt.close(fig1); del fig1, ax


# Plot Error X Percentiles

# Contourf of Error X Forecast Time X Quantiles

# Chiclet Chart variable and error (Contourf of Error X forecast time X time)

# Error X Latitude (satellite)

