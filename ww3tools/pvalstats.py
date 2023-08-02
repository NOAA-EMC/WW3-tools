#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pvalstats.py

VERSION AND LAST UPDATE:
 v1.0  05/24/2022
 v1.1  06/30/2022
 v1.2  11/23/2022
 v1.2  07/17/2023

PURPOSE:
 Group of functionalities for visualization of model/observation data
  and validation.
 Users can import as a standard python module, and use it accordingly:
 For example:

  from pvalstats import ModelObsPlot
  mop=ModelObsPlot(model,obs)
  mop.qqplot()
  mop.scatterplot()
  mop.taylordiagram()
  mop.combinerrors()
  mop.pdf()

  mop=ModelObsPlot(model,obs,dtime)
  mop.timeseries()

USAGE:
 Class ModelObsPlot
 Functions:
   timeseries
   qqplot
   scatterplot
   taylordiagram
   combinerrors
   pdf

 The explanation for each function is contained in the headers, including
  examples,
 help(pvalstats)

OUTPUT:
 png figures saved in the local directory or in the path given through the tag/prefix.

DEPENDENCIES:
 See setup.py and the imports below.

AUTHOR and DATE:
 05/24/2022: Ricardo M. Campos, first version.
 06/30/2022: Ricardo M. Campos, TaylorDiagram and PDF plots corrected
   for arrays containing NaN.
 11/23/2022: Ricardo M. Campos, new functions errXftime and interp_nan,
   plus some format improvements.
 07/17/2023: Ricardo M. Campos, modification to use python Class and OOP.

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pylab import *
from matplotlib.mlab import *
import math
import pandas as pd
from datetime import datetime, timedelta
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import gaussian_kde, linregress
import matplotlib.ticker
from matplotlib.dates import DateFormatter
from mpl_toolkits.basemap import cm
colormap = cm.GMT_polar
palette = plt.cm.jet
palette.set_bad('aqua', 10.0)
import warnings; warnings.filterwarnings("ignore")
# Font size and style
sl=13
matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl) 
matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

def interp_nan(data,lmt=10**50):
    '''
    Fill NaN values with linear interpolation.
    User can enter one or two inputs:
      1) time-series containing NaN values to fill in
      2) maximum number of consecutive NaN values to interpolate (to
         avoid interpolating long segments)
    '''
    data=np.array(data)
    lmt=int(lmt)

    if data.ndim>1:
        raise ValueError(' Input array with too many dimensions. Only time-series (1 dimension) allowed.')
    else:
        # using pandas
        A=pd.Series(data)
        # B=A.interpolate(method="polynomial",order=2,limit=lmt)
        B=A.interpolate(method="linear",limit=lmt)
        B=np.array(B.values)

    return B
    del lmt,A,data

class ModelObsPlot:

    def __init__(self,model=None,obs=None,time=None,axisnames=None,mlabels=[],linreg=None,vaxisname=None,linestyle=None,marker=None,color=None,ftag=''):
        '''
          Initialization of the main object used for the plot functions below.
          model and obs are mandatory.
          time is mandatory for the timeseries plot only.
          Optional: 
          - axisnames adds specific model and observation axis names to qqplot and scatterplot,
              for ex axisnames=["WW3","Buoy"]
          - vaxisname adds specific axis names when only one axis is need, for timeseries (y-axis) and pdf (x-axis),
              for ex vaxisname="Hs(m), WW3 and Buoy"
          - mlabels adds model labels for the legend, to differenciate simulations,
              for ex mlabels=["WW3T1","WW3T2"]
          - ftag is used to include prefix (can include the full path or not),
              for example ftag="/home/ricardo/test/NewRunWW3T2_"
        '''

        if model is None or obs is None:
            raise ValueError("Missing required input arguments.")
            return
        else:
            self.model = np.array(np.atleast_2d(model)).astype('float')
            self.obs = np.array(np.atleast_2d(obs)).astype('float')

        self.ftag = ftag
        self.mlabels = mlabels
        self.linreg = linreg
        self.vaxisname = vaxisname if vaxisname is not None else "Model and Observations"
        self.axisnames = axisnames if axisnames is not None else ["Model","Observations"]
        self.linestyle = (linestyle if linestyle is not None else
            np.array(['-', '--', '-.', '--', '-.', '--', '-.', '--', '-.', '--', '-.', '--', '-.', '--', '-.']))

        self.marker = (marker if marker is not None else
            np.array(np.atleast_1d(['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',])).astype('str'))

        self.color = (color if color is not None else
            np.array(['navy', 'firebrick', 'darkgreen', 'fuchsia', 'gold', 'blue', 'salmon', 'lime', 'darkviolet', 'yellow',
                'cornflowerblue', 'red', 'green', 'violet', 'orange']))

        if self.model.ndim > 2 or self.obs.ndim > 2:
            raise ValueError('Too many dimensions in the model or observation arrays.')

        if self.model.shape[0] > self.model.shape[1]:
            self.model = self.model.T
        if self.obs.shape[0] > self.obs.shape[1]:
            self.obs = self.obs.T

        if self.model.shape[1] != self.obs.shape[1]:
            raise ValueError('Model and observation arrays must have the same index/time size.')

        if (self.model.shape[0]>1) & (np.size(self.mlabels)<=1):
            print(" Warning. You should use labels to differentiate model results.")

        if np.size(time)>1:
            self.time = time
            if self.time.shape[0] != self.model.shape[1]:
                raise ValueError('Time array does not match the model size.')
            else:
                ind=np.where(self.time>-np.inf)
                if np.size(ind)>1:
                    self.model=self.model[:,ind[0]]
                    self.obs=self.obs[:,ind[0]]
                    self.time=self.time[ind[0]]
                else:
                    raise ValueError('Time array full of NaN.')

                if isinstance(self.time[0], datetime):
                    self.dtime=self.time
                elif isinstance(self.time[0], float):
                    epoch_start = datetime(1970, 1, 1)
                    self.dtime = [epoch_start + timedelta(seconds=int(ts)) for ts in self.time]
                else:
                    raise ValueError('Time array format not recognized.')

        if self.linreg:
            from scipy.stats import linregress

    def timeseries(self):
        '''
        Time-series plots of model (lines) and observations (discrete points).
        The x-axis is time and y-axis is the variable and unit of the array provided.
        Inputs:
          Mandatory: model(s); observations; time array (datetime).
           the model array can include one or more model results, through the number of columns, while
           the observation array must be one-dimensional, with the same number of lines as the model and
           the time array.
          Optional: see object construction above.
        Output: png figure saved in the local directory where python is running or in the path given through ftag.
        Example:
          from pvalstats import ModelObsPlot

          mop=ModelObsPlot(ww3gfs,obs,dtime)
          mop.timeseries()

          mop=ModelObsPlot(model=np.c_[model1[:],model2[:]],obs=buoydata[:],time=t[:],
              vaxisname="Hs (m)",mlabels=["WW3T1","WW3T2"],ftag="NewRun_")
          mop.timeseries()

        '''

        if self.dtime is None:
            raise ValueError("Missing time array input argument.")
            return

        fig1, ax = plt.subplots(figsize=(9, 4))
        if np.size(self.mlabels) > 0 and np.size(np.where(self.obs[0, :] > -999)) > 1:
            ax.plot(self.dtime, self.obs[0, :], color='dimgray', marker='.', linestyle='', linewidth=2., label='Obs', zorder=2)
        else:
            ax.plot(self.dtime, self.obs[0, :], color='dimgray', marker='.', linestyle='', linewidth=2., zorder=2)

        if np.size(np.where(self.obs[0, :] > -999)) > 1:
            ax.fill_between(self.dtime, 0., self.obs[0, :], color='silver', alpha=0.5, zorder=1)

        for i in range(self.model.shape[0]):
            if np.size(self.mlabels) > 0:
                if np.size(np.where(self.model[i, :] > -999)) > 1:
                    ax.plot(self.dtime, interp_nan(self.model[i, :], 10), color=self.color[i], linestyle=self.linestyle[i],
                            linewidth=2., label=self.mlabels[i], alpha=0.8, zorder=3)
            else:
                if np.size(np.where(self.model[i, :] > -999)) > 1:
                    ax.plot(self.dtime, interp_nan(self.model[i, :], 10), color=self.color[i], linestyle=self.linestyle[i],
                            linewidth=2., alpha=0.8, zorder=3)

            if np.size(np.where(self.model[i, :] > -999)) > 1:
                ax.plot(self.dtime, self.model[i, :], marker='.', markersize=5, markerfacecolor=self.color[i],
                        markeredgewidth=0.5, markeredgecolor='black', linestyle='', linewidth=0., zorder=3)

        ax.set_xlim(self.dtime[0], self.dtime[-1])
        ax.xaxis.set_major_formatter(DateFormatter('%b%d'))
        ax.fmt_xdata = DateFormatter('%b%d')

        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(sl - 2)

        ax.set_xlabel("Time")
        if np.size(self.vaxisname) > 0:
            ax.set_ylabel(self.vaxisname)
        else:
            ax.set_ylabel("Model and Observations")

        if np.size(self.mlabels) > 0:
            if np.size(self.mlabels) < 7:
                ax.legend(fontsize=sl - 3)
            else:
                ax.legend(fontsize=sl - 5)

        # plt.tight_layout()
        plt.grid(c='grey', ls='--', alpha=0.3, zorder=1)
        plt.savefig(self.ftag + 'TimeSeries.png', dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
                    format='png', bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1)

    def qqplot(self):
        '''
        Quantile-Quantile plot.
        Inputs:
          Mandatory: model(s); observations (or any reference).
           the model array can include one or more model results, through the number of columns, while
           the observation array must be one-dimensional, with the same number of lines as the model.
          Optional: see object construction above.
        Output: png figure saved in the local directory where python is running or in the path given through ftag.
        Example:
          from pvalstats import ModelObsPlot
          mop=ModelObsPlot(ww3gfs,obs)
          mop.qqplot()

          mop=ModelObsPlot(model=np.c_[model1[:],model2[:]],obs=buoydata[:],axisnames=["WW3","Buoy"],
              mlabels=["WW3T1","WW3T2"],ftag="/home/ricardo/testWW3/NewRun_")
          mop.qqplot()
        '''

        # percentiles
        p = np.arange(1,96,1); p=np.append(p,[95.5,96.,96.5,97.,97.3,97.6,97.9,98.2,98.5,98.7,98.9,99.1,99.3,99.5,99.7])
        qm = np.zeros((self.model.shape[0],p.shape[0]),'f')*np.nan
        qobs = np.zeros((self.model.shape[0],p.shape[0]),'f')*np.nan
        for i in range(0,p.shape[0]):
            for j in range(0,self.model.shape[0]):
                qobs[j,i] = np.nanpercentile(self.obs[:],p[i])
                qm[j,i] = np.nanpercentile(self.model[j,:],p[i])

        a=math.floor(np.nanmin([qobs,qm])*100.)/100. ; b=math.ceil(np.nanmax([qobs,qm])*100.)/100.
        aux=np.linspace(a-0.2*a,b+0.2*a,p.shape[0])
	# plot
        fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)
        ax.plot(aux,aux,'k', linewidth=2.,alpha=0.4,zorder=2); del a,b
        for i in range(0,self.model.shape[0]):
            if np.size(self.mlabels)>0:
                ax.plot(gaussian_filter(qobs[i,:],1),gaussian_filter(qm[i,:],1), color=self.color[i], marker=self.marker[i], linestyle='--',label=self.mlabels[i],linewidth=1.,alpha=0.8,zorder=3)
            else:
                ax.plot(gaussian_filter(qobs[i,:],1),gaussian_filter(qm[i,:],1), color=self.color[i], marker=self.marker[i], linestyle='--',linewidth=1.,alpha=0.8,zorder=3)

            if self.linreg:
                a=np.array(qm[i,:]); b=np.array(qobs[i,:]); ind=np.where((a*b)>-999.)
                if np.size(ind)>2:
                    a=np.copy(a[ind]); b=np.copy(b[ind]);
                    r = linregress(b,a)
                    aregr = np.array(r.slope*aux + r.intercept )
                    ax.plot(aux,aregr,color=self.color[i],ls='-',linewidth=1.,alpha=0.8,zorder=4)
                    ax.plot(aux,aregr,color='k',ls=':',linewidth=0.7,alpha=0.7,zorder=4)
                    if np.size(self.mlabels)>0:
                        print(self.ftag+"QQplot "+self.mlabels[i]+": Slope "+np.str(np.round(float(r.slope),5))+", Intercept "+np.str(np.round(float(r.intercept),5)))
                    else:
                        print(self.ftag+"QQplot: Slope "+np.str(np.round(float(r.slope),5))+", Intercept "+np.str(np.round(float(r.intercept),5)))                      

                    del a,b,r,aregr

                del ind

        plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
        for i in np.array([50,80,90,95,99]):
            plt.axvline(x= np.nanpercentile(self.obs,int(i)),ls='--',color='grey',linewidth=1.,alpha=0.9,zorder=1)
            plt.text(np.nanpercentile(self.obs,int(i)),(aux.max()-aux.min())/15 + aux.min(),str(int(i))+'th',color='dimgrey',fontsize=sl-7,zorder=4)
            plt.text(np.nanpercentile(self.obs,int(i)),(aux.max()-aux.min())/1.05 + aux.min(),str(int(i))+'th',color='dimgrey',fontsize=sl-7)

        plt.ylim(ymax = aux.max(), ymin = aux.min())
        plt.xlim(xmax = aux.max(), xmin = aux.min())
        plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7) 

        ax.set_xlabel(self.axisnames[1]); ax.set_ylabel(self.axisnames[0])

        if np.size(self.mlabels)>0:
            plt.legend(loc="upper left",fontsize=sl-2)

        plt.tight_layout()
        plt.savefig(self.ftag+'QQplot.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax

    def scatterplot(self):
        '''
        Scatter plot.
        Inputs:
          Mandatory: model(s); observations (or any reference).
           the model array can include one or more model results, through the number of columns, while
           the observation array must be one-dimensional, with the same number of lines as the model.
          Optional: see object construction above.
        Output: png figure saved in the local directory where python is running or in the path given through ftag.
        Example:
          from pvalstats import ModelObsPlot
          mop=ModelObsPlot(ww3gfs,obs)
          mop.scatterplot()

          mop=ModelObsPlot(model=np.c_[model1[:],model2[:]],obs=buoydata[:],axisnames=["WW3","Buoy"],
              mlabels=["WW3T1","WW3T2"],ftag="/home/ricardo/testWW3/NewRun_")
          mop.scatterplot()
        '''

        if self.obs[0,:].shape[0]>50000:
            sk=int(np.round(np.float(self.obs[0,:].shape[0])/30000.,0))
        else:
            sk=1

        a=math.floor(np.nanmin(np.append(self.obs[:,::sk],self.model[:,::sk]))*100.)/100. ; b=math.ceil(np.nanmax(np.append(self.obs[:,::sk],self.model[:,::sk]))*100.)/100.
        famin=a-0.1*a; famax=b+0.1*a
        aux=np.linspace(famin,famax,100); del a,b

        # plot
        fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)

        for i in range(0,self.model.shape[0]):
            b=np.array(self.obs[0,::sk]); a=np.array(self.model[i,::sk])
            ind=np.where((a*b)>-999.)[0]; a=np.copy(a[ind]); b=np.copy(b[ind]); del ind

            if (a.shape[0]<30) | (self.model.shape[0]>1):
                if np.size(self.mlabels)>0:
                    if self.mlabels[0] != '':
                        ax.scatter(b,a,color=self.color[i],marker=self.marker[i],label=self.mlabels[i],zorder=2)
                else:
                    ax.scatter(b,a,color=self.color[i],marker=self.marker[i],zorder=2)

            elif (np.size(self.color)==1) & (self.model.shape[0]==1):
                ax.scatter(b,a,color=self.color[i],marker=self.marker[i],zorder=2)

            else:
                xy = np.vstack([a,b]); z = gaussian_kde(xy)(xy)
                if np.size(self.mlabels)>0:
                    if self.mlabels[0] != '':
                        ax.scatter(b,a, c=z, s=5,cmap=plt.cm.jet,label=self.mlabels[i],zorder=2)
                else:
                    ax.scatter(b,a, c=z, s=5,cmap=plt.cm.jet,zorder=2)

            if self.linreg:
                r = linregress(b,a)
                aregr = np.array(r.slope*aux + r.intercept )
                ax.plot(aux,aregr,color=self.color[i],ls='-',linewidth=1.,alpha=0.8,zorder=4)
                ax.plot(aux,aregr,color='k',ls=':',linewidth=0.7,alpha=0.7,zorder=4)
                if np.size(self.mlabels)>0:
                    print(self.ftag+"ScatterPlot "+self.mlabels[i]+": Slope "+np.str(np.round(float(r.slope),5))+", Intercept "+np.str(np.round(float(r.intercept),5)))
                else:
                    print(self.ftag+"ScatterPlot: Slope "+np.str(np.round(float(r.slope),5))+", Intercept "+np.str(np.round(float(r.intercept),5)))                      

                del r,aregr

            del a,b

        ax.plot(aux,aux,'k--', linewidth=1.,alpha=0.9,zorder=3)  # main diagonal
        plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7)

        ax.set_xlabel(self.axisnames[1]); ax.set_ylabel(self.axisnames[0])

        plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
        for i in np.array([50,80,95,99]):
            plt.axvline(x= np.nanpercentile(self.obs,int(i)),ls='--',color='grey',linewidth=1.,alpha=0.7,zorder=1)
            plt.text(np.nanpercentile(self.obs,int(i)),((famax-famin)/15)+famin,str(int(i))+'th',color='dimgrey',fontsize=sl-7,zorder=4)
            plt.text(np.nanpercentile(self.obs,int(i)),((famax-famin)/1.05)+famin,str(int(i))+'th',color='dimgrey',fontsize=sl-7,zorder=4)

        plt.gca().set_xlim(left=famin, right=famax); plt.gca().set_ylim(ymin=famin,ymax=famax)

        if np.size(self.mlabels)>0:
            if self.mlabels[0] != '':
                plt.legend(loc="upper left",fontsize=sl-2)

        plt.tight_layout()
        plt.savefig(self.ftag+'ScatterPlot.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax

    def taylordiagram(self):
        '''
        Taylor Diagram using skill_metrics.
        https://github.com/PeterRochford/SkillMetrics
        Inputs:
          Mandatory: model(s); observations (or any reference).
           the model array can include one or more model results, through the number of columns, while
           the observation array must be one-dimensional, with the same number of lines as the model.
          Optional: see object construction above.
        Output: png figure saved in the local directory where python is running or in the path given through ftag.
        Example:
        Example:
          from pvalstats import ModelObsPlot
          mop=ModelObsPlot(ww3gfs,obs)
          mop.taylordiagram()

          mop=ModelObsPlot(model=np.c_[model1[:],model2[:]],obs=buoydata[:],
              mlabels=["WW3T1","WW3T2"],ftag="/home/ricardo/testWW3/NewRun_")
          mop.taylordiagram()
        '''

        import skill_metrics as sm

        # Organizing labels for the Taylor Diagram
        if np.size(self.mlabels)>0:
            label = np.append(['Obs'],self.mlabels)
        else:
            if self.model.shape[0] > 1:
                label = np.array(['Obs'])
                for j in range(0,self.model.shape[0]):
                    label = np.append(label,['Model'+repr(j+1)])
            else:
                label = np.array(['Obs','Model'])

        # remove NaN
        ind=np.where( (np.mean(self.model,axis=0)>-999.) & (np.mean(self.obs,axis=0)>-999.) )
        if np.size(ind)>0:
            model=np.copy(self.model[:,ind[0]]); obs=np.copy(self.obs[:,ind[0]]); del ind
        else:
            raise ValueError(' No quality data available.')

        # plot
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
        if crmsd.max()>0.8: npr=int(2)
        else: npr=int(3)
        tRMS = np.array(np.linspace(0,np.round(crmsd.max()+0.1*crmsd.max(),npr),5)).round(npr)

        if sdev.max()>0.8: npr=int(2)
        else: npr=int(3)
        axmax = np.round(sdev.max()+0.5*sdev.max(),npr)

        if (sdev.max()+0.4*sdev.max())>0.8: npr=int(2)
        else: npr=int(3)
        tSTD = np.array(np.linspace(sdev.min()-0.3*sdev.max(),np.round(sdev.max()+0.4*sdev.max(),npr),4)).round(npr)
        del npr

        fp1=1.07; fp2=1.03; fp3=0.055; fp4=13
        if np.size(label)>8:
            fp2=1.07; fp3=0.040; fp4=12
        elif np.size(label)>5:
            fp2=1.06; fp3=0.045

        slmax=0
        for i in range(0,np.size(label)):
            if len(label[i])>slmax:
                slmax=int(len(label[i]))

        if slmax>11:
            fp1=1.04; fp4=12
        elif slmax>8:
            fp1=1.04

        # loop through models
        for i in range(0,model.shape[0]):
            if i==0:
                sm.taylor_diagram(sdev.take([0,i+1]), crmsd.take([0,i+1]), ccoef.take([0,i+1]), locationColorBar = 'EastOutside',  markerColor = self.color[i],
                    # markerLegend = 'on', markerLabel = list(label.take([0,i+1])),
                    styleOBS = '-', colOBS = 'k', markerobs = 'o', titleOBS = 'Obs', 
                    markerSize = 17, tickRMS = tRMS, axismax = axmax, alpha=0.7, tickSTD = tSTD,
                    tickRMSangle = 87, showlabelsRMS = 'on',
                    colRMS='dimgrey', styleRMS=':', widthRMS=1.0, titleRMS='off',
                    colSTD='k', styleSTD='-.', widthSTD=1.0, titleSTD ='on',
                    colCOR='k', styleCOR='--', widthCOR=1.0, titleCOR='on')

            else:
                sm.taylor_diagram(sdev.take([0,i+1]), crmsd.take([0,i+1]), ccoef.take([0,i+1]), overlay = 'on', locationColorBar = 'EastOutside',  markerColor = self.color[i],
                    # markerLegend = 'on', markerLabel = list(label.take([0,i+1])), 
                    styleOBS = '-', colOBS = 'k', markerobs = 'o', titleOBS = 'Obs', 
                    markerSize = 17, tickRMS = tRMS, axismax = axmax, alpha=0.5, tickSTD = tSTD,
                    tickRMSangle = 87, showlabelsRMS = 'on',
                    colRMS='dimgrey', styleRMS=':', widthRMS=1.0, titleRMS='off',
                    colSTD='k', styleSTD='-.', widthSTD=1.0, titleSTD ='on',
                    colCOR='k', styleCOR='--', widthCOR=1.0, titleCOR='on')

            if np.size(self.mlabels)>0: 
                ax.text(fp1, fp2-(fp3*i), label[i+1], verticalalignment='bottom', horizontalalignment='right',
                    transform=ax.transAxes,color=self.color[i], fontsize=fp4)

        plt.tight_layout()
        plt.savefig(self.ftag+'TaylorDiagram.png', dpi=250, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax

    def combinerrors(self):
        '''
        Combined Errors, Scatter Index and Normalized Bias
        Inputs:
          Mandatory: model(s); observations (or any reference).
           the model array can include one or more model results, through the number of columns, while
           the observation array must be one-dimensional, with the same number of lines as the model.
          Optional: see object construction above.
        Output: png figure saved in the local directory where python is running or in the path given through ftag.
        Example:
          from pvalstats import ModelObsPlot
          mop=ModelObsPlot(ww3gfs,obs)
          mop.combinerrors()

          mop=ModelObsPlot(model=np.c_[model1[:],model2[:]],obs=buoydata[:],axisnames=["WW3","Buoy"],
              mlabels=["WW3T1","WW3T2"],ftag="/home/ricardo/testWW3/NewRun_")
          mop.combinerrors()
        '''

        import mvalstats

        fig1 = plt.figure(1,figsize=(5,4.5)); ax = fig1.add_subplot(111)
        fauxnb=[];fauxsi=[];
        for i in range(0,self.model.shape[0]):
            aux = mvalstats.metrics(self.model[i,:],self.obs[0,:])
            if np.size(self.mlabels)>0:
                ax.scatter(aux[2],aux[5], s=40+150*aux[7]**5,color=self.color[i],edgecolors='k',linewidth=1,alpha=0.7,label=self.mlabels[i],zorder=3)
            else:
                ax.scatter(aux[2],aux[5], s=40+150*aux[7]**5,color=self.color[i],edgecolors='k',linewidth=1,alpha=0.7,zorder=3)

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
        if np.size(self.mlabels)>0:
            plt.legend(fontsize=sl-2)

        plt.tight_layout()
        plt.savefig(self.ftag+'NBiasSI.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax

    def pdf(self):

        '''
        Probability Density Plot. Observations (shaded grey) and model(s)
        Inputs:
          Mandatory: model(s); observations (or any reference).
           the model array can include one or more model results, through the number of columns, while
           the observation array must be one-dimensional, with the same number of lines as the model.
          Optional: see object construction above.

        Output: two png figure saved in the local directory where python is running or in the path given through ftag.
        Example:
          from pvalstats import ModelObsPlot
          mop=ModelObsPlot(ww3gfs,obs)
          mop.pdf()

          mop=ModelObsPlot(model=np.c_[model1[:],model2[:]],obs=buoydata[:],vaxisname=["Hs(m), WW3 and Buoy"],
              mlabels=["WW3T1","WW3T2"],ftag="/home/ricardo/testWW3/NewRun_")
          mop.pdf()
        '''

        from scipy.stats import gaussian_kde

        # remove NaN
        ind=np.where( (np.mean(self.model,axis=0)>-999.) & (np.mean(self.obs,axis=0)>-999.) )
        if np.size(ind)>0:
            model=np.copy(self.model[:,ind[0]]); obs=np.copy(self.obs[:,ind[0]]); del ind
        else:
            raise ValueError(' No quality data available.')

        # plot
        px=np.linspace(np.nanmin([np.nanmin(obs),np.nanmin(model)]),np.nanmax([np.nanmax(obs),np.nanmax(model)])*0.99,100)

        # plt.close('all')
        fig1 = plt.figure(1,figsize=(5,5)); ax = fig1.add_subplot(111)
        dx = gaussian_kde(obs[0,:])
        ax.fill_between(px, 0., gaussian_filter(dx(px), 1.), color='grey', alpha=0.7,zorder=1)
        ax.plot(px,gaussian_filter(dx(px), 1.),color='dimgrey',linestyle='-',linewidth=0.5,alpha=0.7,zorder=2)
        mmax=np.array([1.])
        for i in range(0,model.shape[0]):
            de = gaussian_kde(model[i,:])
            if np.size(self.mlabels)>0:
                ax.plot(px,gaussian_filter(de(px), 1.),color=self.color[i],marker=self.marker[i],linestyle='-',linewidth=1.,alpha=0.7,label=self.mlabels[i],zorder=3)
            else:
                ax.plot(px,gaussian_filter(de(px), 1.),color=self.color[i],marker=self.marker[i],linestyle='-',linewidth=1.,alpha=0.7,zorder=3)

            mmax=np.append(mmax,np.nanmax(gaussian_filter(de(px), 1.)))
            del de

        plt.grid(c='grey', ls=':', alpha=0.5,zorder=1)
        plt.locator_params(axis='y', nbins=7) ; plt.locator_params(axis='x', nbins=7) 
        ax.set_ylabel("Probability Density")
        ax.set_xlabel(self.vaxisname)

        if np.size(self.mlabels)>0:
            plt.legend(fontsize=sl-2)

        plt.tight_layout(); plt.axis('tight'); plt.ylim(ymin = -0.005) 
        plt.savefig(self.ftag+'PDF.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        ax.set_yscale( "log" ); plt.axis('tight');plt.ylim(ymin = 0.001);plt.ylim(ymax = np.nanmax(mmax))
        plt.savefig(self.ftag+'PDFlogscale.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax, mmax

