#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pvalstats.py

VERSION AND LAST UPDATE:
 v1.0  05/24/2022
 v1.1  06/30/2022
 v1.2  11/23/2022
 v1.2  07/17/2023
 v1.3  08/03/2023
 v1.4  08/04/2023

PURPOSE:
 Group of functionalities for visualization of model/observation data
  and validation.
 Users can import as a standard python module, and use it accordingly:
- For example (ModelObsPlot):

  from pvalstats import ModelObsPlot
  mop=ModelObsPlot(model,obs)
  mop.qqplot()
  mop.scatterplot()
  mop.taylordiagram()
  mop.combinerrors()
  mop.pdf()

  mop=ModelObsPlot(model,obs,dtime)
  mop.timeseries()

  fctintervals=np.arange(0,1+np.ceil(np.max(forecastLeadTime)/24)*24,24).astype('int')
  fctxt=np.array(frintervalsd/24).astype('int')[0:-1]+1 # x-axis plot array
  fctxticks=np.array([1,2,3,4,5,7,10,15,20,25,30,35]) # x-axis figure ticks
  outpath="/home/ricardo/validation/"
  mop=ModelObsPlot(np.c_[modelControl,modelEnsMean].T,obs=obs,
    fctime=forecastLeadTime,fctintervals=fctintervals,fctxt=fctxt,
    fctxticks=fctxticks,fctunits="days",
    ftag=outpath+"Results_Hs_",vaxisname="Hs")

  mop.errxfctime()

- For example (GlobalMapPlot):
  from pvalstats import GlobalMapPlot, gsmooth

  gdata = np.random.rand(10, 10)  # Example data
  lonm = np.linspace(0, 360, 10)  # Example longitude values
  latm = np.linspace(-90, 90, 10)  # Example latitude values
  gmp = GlobalMapPlot(gdata=gsmooth(gdata,1),lat=latm,lon=lonm)
  gmp.plot()

  outpath="/home/ricardo/validation/"
  figname=outpath+"GlobalMap_RMSE_Hs.png"
  palette = plt.cm.jet; pextend="max"
  levels = np.array(np.linspace(0.1,2.3,31)).round(2)
  gmp = GlobalMapPlot(gdata=gsmooth(RmseGefsControl,1),latm=lat,lonm=lon,
    latmin=-67.,latmax=67.,palette = palette, levels=levels,ftag=figname,pextend=pextend)
  gmp.plot()   

USAGE:
 Class ModelObsPlot
 Functions:
   timeseries
   qqplot
   scatterplot
   taylordiagram
   combinerrors
   pdf
   monthlystats
   errxfctime
   rankhist
   spreadxfctime
   crpsxfctime
   brierxfctime
   rocaucxfctime

 Class GlobalMapPlot
 Functions:
   gsmooth
   plot

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
 08/03/2023: Ricardo M. Campos, new functions included rankhist,monthlystats,
   errxfctime,spreadxfctime,crpsxfctime
 08/04/2023: Ricardo M. Campos, new class GlobalMapPlot and function gsmooth

PERSON OF CONTACT:
 Ricardo M Campos: ricardo.campos@noaa.gov

"""

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
from datetime import datetime, timedelta
from scipy.ndimage.filters import gaussian_filter
from scipy.stats import gaussian_kde, linregress
import matplotlib.ticker
from matplotlib.dates import DateFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
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

    def __init__(self,model=None,obs=None,time=None,
        fctime=None,fctintervals=None,fctxt=None,fctxticks=None,fctunits=None,month=None,linreg=None,axisnames=None,
        mlabels=[],vaxisname=None,linestyle=None,linestyle_html=None,marker=None,color=None,ftag='',figformat=None):

        """
        Initialization of the main object used for the plot functions below.
        model and obs are mandatory.
        time is mandatory for the timeseries plot only.
        Optional: 
        :param axisnames: adds specific model and observation axis names to qqplot and scatterplot,
            for ex axisnames=["WW3","Buoy"]
        :param vaxisname: adds specific axis names when only one axis is need, for timeseries (y-axis), month (y-axis), and pdf (x-axis),
            for ex vaxisname="Hs(m)" or vaxisname="Hs(m), WW3 and Buoy"
        :param mlabels: adds model labels for the legend, to differenciate simulations,
            for ex mlabels=["WW3T1","WW3T2"]
        :param ftag: is used to include prefix (can include the full path or not),
            for example ftag="/home/ricardo/test/NewRunWW3T2_"

        """

        if model is None or obs is None:
            raise ValueError("Missing required input arguments.")
            return
        else:
            self.model = np.array(np.atleast_2d(model)).astype('float') 
            self.obs = np.array(np.atleast_2d(obs)).astype('float')

        self.figformat = figformat
        self.ftag = ftag
        self.mlabels = mlabels
        self.linreg = linreg
        self.vaxisname = vaxisname if vaxisname is not None else "Model and Observations"
        self.axisnames = axisnames if axisnames is not None else ["Model","Observations"]
        self.linestyle = (linestyle if linestyle is not None else
            np.array(['-', '--', '-.', '--', '-.', '--', '-.', '--', '-.', '--', '-.', '--', '-.', '--', '-.']))

        self.linestyle_html = (linestyle_html if linestyle_html is not None else
            np.array(['solid','dashed','dashdot','dashed','dashdot','dashed','dashdot','dashed','dashdot','dashed','dashdot','dashed','dashdot','dashed','dashdot']))

        self.marker = (marker if marker is not None else
            np.array(np.atleast_1d(['.','.','.','.','.','.','.','.','.','.','.','.','.','.','.',])).astype('str'))

        self.color = (color if color is not None else
            np.array(['darkblue', 'darkred', 'darkgreen', 'darkorange', 'purple', 'deeppink', 'brown', 'salmon', 'lime', 'darkviolet', 'yellow',
                'cornflowerblue', 'red', 'green', 'violet', 'orange']))

        if self.model.ndim > 2 or self.obs.ndim > 2:
            raise ValueError('Too many dimensions in the model or observation arrays.')

        if self.model.shape[0] > self.model.shape[1]:
            self.model = self.model.T
        if self.obs.shape[0] > self.obs.shape[1]:
            self.obs = self.obs.T

        if self.model.shape[1] != self.obs.shape[1]:
            raise ValueError('Model and observation arrays must have the same index/time size.')

        if np.size(fctime)>1:
            self.fctime = fctime
            if self.fctime.shape[0] != self.model.shape[1]:
                raise ValueError('Forecast time array does not match the model/obs size.')

        if np.size(month)>1:
            self.fmonth = month
            if self.fmonth.shape[0] != self.model.shape[1]:
                raise ValueError('Month array does not match the model/obs size.')

        if np.size(time)>1:
            self.time = time
            if self.time.shape[0] != self.model.shape[1]:
                raise ValueError('Time array does not match the model/obs size.')
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

        self.fctintervals = fctintervals
        self.fctxt = fctxt
        self.fctunits = fctunits
        self.fctxticks = fctxticks
        # name of the 9 error metrics
        self.nerrm=np.array(['bias','RMSE','NBias','NRMSE','SCrmse','SI','HH','CC','N'])

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
                    ax.plot(self.dtime, interp_nan(self.model[i, :], 3), color=self.color[i], linestyle=self.linestyle[i],
                            linewidth=2., label=self.mlabels[i], alpha=0.8, zorder=3)
            else:
                if np.size(np.where(self.model[i, :] > -999)) > 1:
                    ax.plot(self.dtime, interp_nan(self.model[i, :], 3), color=self.color[i], linestyle=self.linestyle[i],
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

        # Html plot with bokeh
        if self.figformat == 'html':
            from bokeh.plotting import figure, output_file, show
            plt.close('all')
            p = figure(x_axis_type="datetime", plot_width=1400, plot_height=600)
            p.title.text = 'Click on legend entries to hide the corresponding lines'
            p.line(self.dtime,self.obs[0, :],line_width=3, color='black', alpha=0.8, legend_label='Obs')
            for i in range(self.model.shape[0]):
                if np.size(self.mlabels) > 0:
                    if np.size(np.where(self.model[i, :] > -999)) > 1:
                        p.line(self.dtime,interp_nan(self.model[i, :], 3), line_width=2, color=self.color[i], line_dash=self.linestyle_html[i], alpha=0.8, legend_label=self.mlabels[i])

                else:
                    if np.size(np.where(self.model[i, :] > -999)) > 1:
                        p.line(self.dtime,interp_nan(self.model[i, :], 3), line_width=2, color=self.color[i], line_dash=self.linestyle_html[i], alpha=0.8)

            if np.size(self.mlabels) > 0:
                p.legend.location = "top_left" ; p.legend.click_policy = "hide"

            p.xaxis.axis_label = "Time"

            if np.size(self.vaxisname) > 0:
                p.yaxis.axis_label = self.vaxisname
            else:
                p.yaxis.axis_label = "Model and Observations"
         
            output_file(self.ftag + 'TimeSeries.html'); show(p); del p

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


    def scatterplot(self, dwscl='no'):
        '''
        Scatter plot.
        Inputs:
          Mandatory: model(s); observations (or any reference).
           the model array can include one or more model results, through the number of columns, while
           the observation array must be one-dimensional, with the same number of lines as the model.
          Optional: see object construction above.
          - mlabels: List containing labels for each model. Default is an empty list.
          - ftag: Path to save the figure. Default is the current directory.
          - dwscl: Whether to apply downsampling. Default is 'no'. If set to 'no', downsampling is not applied.
                If set to any other value, downsampling is applied.
        Output: png figure saved in the local directory where python is running or in the path given through ftag.
        Example:
          from pvalstats import ModelObsPlot
          mop=ModelObsPlot(ww3gfs,obs)
          mop.scatterplot()

          mop=ModelObsPlot(model=np.c_[model1[:],model2[:]],obs=buoydata[:],axisnames=["WW3","Buoy"],
              mlabels=["WW3T1","WW3T2"],ftag="/home/ricardo/testWW3/NewRun_")
          mop.scatterplot(dwscl='yes') yes=downsampling; no: No downsampling
          mop.scatterplot()
        '''

        if self.obs[0,:].shape[0] > 50000 and dwscl != 'no':
            sk = int(np.round(float(self.obs[0,:].shape[0]) / 30000., 0))
        else:
            sk = 1

        a = math.floor(np.nanmin(np.append(self.obs[:,::sk], self.model[:,::sk])) * 100.) / 100.
        b = math.ceil(np.nanmax(np.append(self.obs[:,::sk], self.model[:,::sk])) * 100.) / 100.
        famin = a - 0.1 * a
        famax = b + 0.1 * a
        aux = np.linspace(famin, famax, 100)


        # plot
        fig1 = plt.figure(1, figsize=(5, 4.5))
        ax = fig1.add_subplot(111)

        for i in range(0, self.model.shape[0]):
            b = np.array(self.obs[0,::sk])
            a = np.array(self.model[i,::sk])
            ind = np.where((a*b) > -999.)[0]
            a = np.copy(a[ind])
            b = np.copy(b[ind])

            if (a.shape[0] < 30) or (self.model.shape[0] > 1):
                if np.size(self.mlabels) > 0:
                    if self.mlabels[0] != '':
                        ax.scatter(b, a, color=self.color[i], marker=self.marker[i], label=self.mlabels[i], zorder=2)
                    else:
                        ax.scatter(b, a, color=self.color[i], marker=self.marker[i], zorder=2)
                else:
                    ax.scatter(b, a, color=self.color[i], marker=self.marker[i], zorder=2)
            elif (np.size(self.color) == 1) and (self.model.shape[0] == 1):
                ax.scatter(b, a, color=self.color[i], marker=self.marker[i], zorder=2)
            else:
                xy = np.vstack([a, b])
                z = gaussian_kde(xy)(xy)
                if np.size(self.mlabels) > 0:
                    if self.mlabels[0] != '':
                        ax.scatter(b, a, c=z, s=5, cmap=plt.cm.jet, label=self.mlabels[i], zorder=2)
                    else:
                        ax.scatter(b, a, c=z, s=5, cmap=plt.cm.jet, zorder=2)
                else:
                    ax.scatter(b, a, c=z, s=5, cmap=plt.cm.jet, zorder=2)

            if self.linreg:
                r = linregress(b, a)
                aregr = np.array(r.slope * aux + r.intercept)
                ax.plot(aux, aregr, color=self.color[i], ls='-', linewidth=1., alpha=0.8, zorder=4)
                ax.plot(aux, aregr, color='k', ls=':', linewidth=0.7, alpha=0.7, zorder=4)
                if np.size(self.mlabels) > 0:
                    print(self.ftag + "ScatterPlot " + self.mlabels[i] + ": Slope " + np.str(np.round(float(r.slope), 5))
                          + ", Intercept " + np.str(np.round(float(r.intercept), 5)))
                else:
                    print(self.ftag + "ScatterPlot: Slope " + np.str(np.round(float(r.slope), 5)) + ", Intercept "
                          + np.str(np.round(float(r.intercept), 5)))
                del r, aregr

            del a, b

        ax.plot(aux, aux, 'k--', linewidth=1., alpha=0.9, zorder=3)  # main diagonal
        plt.locator_params(axis='y', nbins=7)
        plt.locator_params(axis='x', nbins=7)

        ax.set_xlabel(self.axisnames[1])
        ax.set_ylabel(self.axisnames[0])

        plt.grid(c='grey', ls=':', alpha=0.5, zorder=1)
        for i in np.array([50, 80, 95, 99]):
            plt.axvline(x=np.nanpercentile(self.obs, int(i)), ls='--', color='grey', linewidth=1., alpha=0.7, zorder=1)
            plt.text(np.nanpercentile(self.obs, int(i)), ((famax - famin) / 15) + famin, str(int(i)) + 'th', color='dimgrey',
                     fontsize=sl - 7, zorder=4)
            plt.text(np.nanpercentile(self.obs, int(i)), ((famax - famin) / 1.05) + famin, str(int(i)) + 'th',
                     color='dimgrey', fontsize=sl - 7, zorder=4)

        plt.gca().set_xlim(left=famin, right=famax)
        plt.gca().set_ylim(ymin=famin, ymax=famax)

        if np.size(self.mlabels) > 0:
            if self.mlabels[0] != '':
                plt.legend(loc="upper left", fontsize=sl - 2)

        plt.tight_layout()
        plt.savefig(self.ftag + 'ScatterPlot.png', dpi=200, facecolor='w', edgecolor='w', orientation='portrait',
                    format='png', transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1)
        del fig1, ax

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

    def pdf(self,result="no"):
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


    def monthlystats(self,result="no"):

        if (np.size(self.fmonth)!=self.obs.shape[1]):
            raise ValueError(' Array with months should have the same length as the array of observations and model(s).')

        self.fmonth=np.array(self.fmonth).astype('int')
        import mvalstats
        months=np.arange(1,13,1)
        monthlyval=np.zeros((months.shape[0],self.model.shape[0],np.size(self.nerrm)),'f')*np.nan
        for i in range(0,months.shape[0]):
            indm=np.where((self.fmonth==months[i]) & (np.mean(self.model[:,:],axis=0)>-999) & (np.mean(self.obs[:,:],axis=0)>-999))
            if np.size(indm)>1:
                indm=indm[0]
                for j in range(0,self.model.shape[0]):
                    monthlyval[i,j,:] = mvalstats.metrics(self.model[j,indm],self.obs[0,indm])

        for i in range(0,len(self.nerrm)):

            fig1 = plt.figure(1,figsize=(9,4)); ax = fig1.add_subplot(111)
            for j in range(0,self.model.shape[0]):
                if np.size(self.mlabels)>0:
                    if self.mlabels[0] != '':
                        ax.plot(months,monthlyval[:,j,i],self.color[j],marker=self.marker[j],label=self.mlabels[j], linewidth=2.,zorder=2)
                    else:
                        ax.plot(months,monthlyval[:,j,i],self.color[j],marker=self.marker[j],linewidth=2.,zorder=2)

                else:
                    ax.plot(months,monthlyval[:,j,i],self.color[j],marker=self.marker[j],linewidth=2.,zorder=2)

            ax.set_xticks(months)
            ax.set_xlabel('Month',size=sl)

            if np.size(self.mlabels)>0:
                plt.legend(fontsize=sl-2)

            if self.vaxisname != "Model and Observations":
                ax.set_ylabel(self.nerrm[i]+" "+self.vaxisname)
            else:
                ax.set_ylabel(self.nerrm[i])

            plt.tight_layout();plt.axis('tight')
            ax.grid(c='dimgrey', ls='--', alpha=0.3)
            plt.xlim(xmin = 0.9, xmax = 12.1)
            plt.savefig(self.ftag+"MonthlyError_"+self.nerrm[i]+".png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
            plt.close(fig1); del fig1, ax

        if result!="no":
            print(" Returning two arrays with error_metrics and metrics_names, from ModelObsPlot.monthlystats")
            print(" The error_metrics array has dimensions [month,models,metrics_names]")
            return monthlyval,self.nerrm

    def errxfctime(self,result="no"):

        if np.size(self.fctime)<=1:
            raise ValueError('Forecast time array must be informed.')

        import mvalstats

        if np.size(self.fctintervals)<=1:
            self.fctintervals = np.unique(self.fctime)

        if np.size(self.fctxt)<=1:
            self.fctxt=[]
            for i in range(len(self.fctintervals) - 1):
                self.fctxt.append(f"{self.fctintervals[i]}-{self.fctintervals[i+1]}")

            self.fctxt=np.array(self.fctxt)

        if self.fctxt.shape[0] == (self.fctintervals.shape[0]):
            valarr=0
        elif self.fctxt.shape[0] == (self.fctintervals.shape[0]-1):
            valarr=1
        else:
            raise ValueError('The size of fctxt must be equal to np.size(fctintervals) or np.size(fctintervals-1) ')

        errm=np.zeros((self.fctxt.shape[0],self.model.shape[0],np.size(self.nerrm)),'f')*np.nan
        for i in range(0,self.fctxt.shape[0]):
            if valarr==1:
                ind=np.where( (self.fctime>=self.fctintervals[i]) & (self.fctime<=self.fctintervals[i+1]) )
            else:
                ind=np.where(self.fctime==self.fctintervals[i])

            if np.size(ind)>1:
                ind=ind[0]
                obs=np.array(self.obs[0,ind])
                for j in range(0,self.model.shape[0]):
                    errm[i,j,:] = mvalstats.metrics(self.model[j,ind],obs)

                del obs

            del ind

        # plot
        for i in range(0,len(self.nerrm)):

            fig1 = plt.figure(1,figsize=(9,4)); ax = fig1.add_subplot(111)
            for j in range(0,self.model.shape[0]):
                if np.size(self.mlabels)>0:
                    if self.mlabels[0] != '':
                        ax.plot(self.fctxt,errm[:,j,i],self.color[j],marker=self.marker[j],label=self.mlabels[j], linewidth=2.,zorder=2)
                    else:
                        ax.plot(self.fctxt,errm[:,j,i],self.color[j],marker=self.marker[j],linewidth=2.,zorder=2)

                else:
                    ax.plot(self.fctxt,errm[:,j,i],self.color[j],marker=self.marker[j],linewidth=2.,zorder=2)

            ax.set_xticks(self.fctxt)

            if self.fctunits:
                ax.set_xlabel('Forecast Lead Time ('+self.fctunits+')',size=sl)
            else:
                ax.set_xlabel('Forecast Lead Time ',size=sl)

            if self.vaxisname != "Model and Observations":
                ax.set_ylabel(self.nerrm[i]+" "+self.vaxisname)
            else:
                ax.set_ylabel(self.nerrm[i])

            if np.size(self.mlabels)>0:
                plt.legend(loc='best',fontsize=sl-2)

            ax.grid(c='dimgrey', ls='--', alpha=0.3)
            plt.tight_layout();plt.axis('tight')

            if np.size(self.fctxticks)>1:
                ax.set_xticks(self.fctxticks)
                plt.xlim(xmin = np.min(self.fctxticks)*0.9, xmax = np.max(self.fctxticks)*1.01)

            plt.savefig(self.ftag+"ErrXForecastTime_"+self.nerrm[i]+".png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
            plt.close(fig1); del fig1, ax

        if result!="no":
            print(" Returning two arrays with error_metrics and metrics_names, from ModelObsPlot.errxfctime")
            print(" The error_metrics array has dimensions [forecastLeadTime,models,metrics_names]")
            return errm,self.nerrm

    # --- Ensemble Validation ---

    # Rank Histogram (Talagrand diagram) 
    # https://journals.ametsoc.org/view/journals/mwre/129/3/1520-0493_2001_129_0550_iorhfv_2.0.co_2.xml
    def rankhist(self,result="no"):

        from scipy.stats import rankdata
        # remove NaN
        ind=np.where( (np.mean(self.model,axis=0)>-999.) & (np.mean(self.obs,axis=0)>-999.) )
        if np.size(ind)>0:
            model=np.copy(self.model[:,ind[0]]); obs=np.copy(self.obs[0,ind[0]]); del ind
        else:
            raise ValueError(' No quality data available.')

        combined=np.vstack((obs[np.newaxis],model))
        ranks=np.apply_along_axis(lambda x: rankdata(x,method='min'),0,combined)
        ties=np.sum(ranks[0]==ranks[1:], axis=0)
        ranks=ranks[0]
        tie=np.unique(ties)
        for i in range(1,len(tie)):
            index=ranks[ties==tie[i]]
            ranks[ties==tie[i]]=[np.random.randint(index[j],index[j]+tie[i]+1,tie[i])[0] for j in range(len(index))]

        result = np.histogram(ranks, bins=np.linspace(0.5, combined.shape[0]+0.5, combined.shape[0]+1))
        # plot
        fig1 = plt.figure(1,figsize=(6,5)); ax1 = fig1.add_subplot(111)
        fresult=result[0]/np.sum(result[0])
        rxaxis=range(1,len(result[0])+1)
        plt.bar(rxaxis,fresult,color='grey', edgecolor='k')
        # plt.ylim(ymax = 1.01, ymin = -0.01)
        ax1.set_xticks(np.arange(1,13,1))
        ax1.set_xlabel('Rank',size=sl); ax1.set_ylabel('Probability',size=sl)
        plt.tight_layout(); # plt.axis('tight') 
        plt.grid(c='grey', ls='--', alpha=0.3)
        plt.savefig(self.ftag+'RankHistogram.png', dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax1, result

        if result!="no":
            print(" Returning two arrays with probabilities (y-axis) and ensemble_members+1 (x-axis), from ModelObsPlot.rankhist")
            return fresult,rxaxis


    def spreadxfctime(self,result="no"):

        if np.size(self.fctime)<=1:
            raise ValueError('Forecast time array must be informed.')

        if self.model.shape[0]<=1:
            raise ValueError('Forecast must have at least 2 members.')

        if np.size(self.fctintervals)<=1:
            self.fctintervals = np.unique(self.fctime)

        if np.size(self.fctxt)<=1:
            self.fctxt=[]
            for i in range(len(self.fctintervals) - 1):
                self.fctxt.append(f"{self.fctintervals[i]}-{self.fctintervals[i+1]}")

            self.fctxt=np.array(self.fctxt)

        if self.fctxt.shape[0] == (self.fctintervals.shape[0]):
            valarr=0
        elif self.fctxt.shape[0] == (self.fctintervals.shape[0]-1):
            valarr=1
        else:
            raise ValueError('The size of fctxt must be equal to np.size(fctintervals) or np.size(fctintervals-1) ')

        mspread=np.zeros((self.fctxt.shape[0]),'f')*np.nan
        for i in range(0,self.fctxt.shape[0]):
            if valarr==1:
                ind=np.where( (self.fctime>=self.fctintervals[i]) & (self.fctime<=self.fctintervals[i+1]) )
            else:
                ind=np.where(self.fctime==self.fctintervals[i])

            if np.size(ind)>1:
                ind=ind[0]
                # unbiased std using Bessel's correction (ddof=1)
                mspread[i]=np.nanmean(np.nanstd(self.model[:,ind],axis=0,ddof=1))

        # plot Spread
        fig1 = plt.figure(1,figsize=(9,4)); ax = fig1.add_subplot(111)
        ax.plot(self.fctxt,mspread,'k',marker='.',linewidth=2.,zorder=2)
        ax.set_xticks(self.fctxt)
        if self.fctunits:
            ax.set_xlabel('Forecast Lead Time ('+self.fctunits+')',size=sl)
        else:
            ax.set_xlabel('Forecast Lead Time ',size=sl)

        if self.vaxisname != "Model and Observations":
            ax.set_ylabel("Spread "+self.vaxisname)
        else:
            ax.set_ylabel("Spread")

        ax.grid(c='dimgrey', ls='--', alpha=0.3)
        plt.tight_layout();plt.axis('tight')

        if np.size(self.fctxticks)>1:
            ax.set_xticks(self.fctxticks)
            plt.xlim(xmin = np.min(self.fctxticks)*0.9, xmax = np.max(self.fctxticks)*1.01)

        plt.savefig(self.ftag+"SpreadXForecastTime.png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax

        if result!="no":
            print(" Returning SPREAD array, from ModelObsPlot.spreadxfctime")
            return mspread


    def crpsxfctime(self,result="no"):

        import properscoring as ps

        if np.size(self.fctime)<=1:
            raise ValueError('Forecast time array must be informed.')

        if self.model.shape[0]<=1:
            raise ValueError('Forecast must have at least 2 members.')

        if np.size(self.fctintervals)<=1:
            self.fctintervals = np.unique(self.fctime)

        if np.size(self.fctxt)<=1:
            self.fctxt=[]
            for i in range(len(self.fctintervals) - 1):
                self.fctxt.append(f"{self.fctintervals[i]}-{self.fctintervals[i+1]}")

            self.fctxt=np.array(self.fctxt)

        if self.fctxt.shape[0] == (self.fctintervals.shape[0]):
            valarr=0
        elif self.fctxt.shape[0] == (self.fctintervals.shape[0]-1):
            valarr=1
        else:
            raise ValueError('The size of fctxt must be equal to np.size(fctintervals) or np.size(fctintervals-1) ')

        crps=np.zeros((self.fctxt.shape[0]),'f')*np.nan
        for i in range(0,self.fctxt.shape[0]):
            if valarr==1:
                ind=np.where( (self.fctime>=self.fctintervals[i]) & (self.fctime<=self.fctintervals[i+1]) )
            else:
                ind=np.where(self.fctime==self.fctintervals[i])

            if np.size(ind)>1:
                ind=ind[0]
                crps[i]=np.nanmean(ps.crps_ensemble( np.array(self.obs[0,ind]),self.model[:,ind].T))

        # plot CRPS
        fig1 = plt.figure(1,figsize=(9,4)); ax = fig1.add_subplot(111)
        ax.plot(self.fctxt,crps,'k',marker='.',linewidth=2.,zorder=2)
        ax.set_xticks(self.fctxt)
        if self.fctunits:
            ax.set_xlabel('Forecast Lead Time ('+self.fctunits+')',size=sl)
        else:
            ax.set_xlabel('Forecast Lead Time ',size=sl)

        if self.vaxisname != "Model and Observations":
            ax.set_ylabel("CRPS "+self.vaxisname)
        else:
            ax.set_ylabel("CRPS")

        ax.grid(c='dimgrey', ls='--', alpha=0.3)
        plt.tight_layout();plt.axis('tight')

        if np.size(self.fctxticks)>1:
            ax.set_xticks(self.fctxticks)
            plt.xlim(xmin = np.min(self.fctxticks)*0.9, xmax = np.max(self.fctxticks)*1.01)

        plt.savefig(self.ftag+"CRPSXForecastTime.png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax

        if result!="no":
            print(" Returning CRPS array, from ModelObsPlot.crpsxfctime")
            return crps

    # https://www.psl.noaa.gov/people/tom.hamill/ROC_problems.pdf
    def brierxfctime(self,threshold=None,result="no"):

        from sklearn.metrics import brier_score_loss

        if np.size(self.fctime)<=1:
            raise ValueError('Forecast time array must be informed.')

        if self.model.shape[0]<=1:
            raise ValueError('Forecast must have at least 2 members.')

        if np.size(self.fctintervals)<=1:
            self.fctintervals = np.unique(self.fctime)

        if np.size(self.fctxt)<=1:
            self.fctxt=[]
            for i in range(len(self.fctintervals) - 1):
                self.fctxt.append(f"{self.fctintervals[i]}-{self.fctintervals[i+1]}")

            self.fctxt=np.array(self.fctxt)

        if self.fctxt.shape[0] == (self.fctintervals.shape[0]):
            valarr=0
        elif self.fctxt.shape[0] == (self.fctintervals.shape[0]-1):
            valarr=1
        else:
            raise ValueError('The size of fctxt must be equal to np.size(fctintervals) or np.size(fctintervals-1) ')

        if threshold==None:
            threshold=np.nanmean(self.obs[0,:])

        num_ensemble_members = int(self.model.shape[0])
        prob_forecasts = [1/num_ensemble_members] * num_ensemble_members

        brier=np.zeros((self.fctxt.shape[0]),'f')*np.nan
        for i in range(0,self.fctxt.shape[0]):
            if valarr==1:
                ind=np.where( (self.fctime>=self.fctintervals[i]) & (self.fctime<=self.fctintervals[i+1]) )
            else:
                ind=np.where(self.fctime==self.fctintervals[i])

            if np.size(ind)>1:
                ind=ind[0]
                if self.obs[0,ind].shape[0]>50000:
                    sk=int(np.round(np.float(self.obs[0,ind].shape[0])/30000.,0))
                else:
                    sk=1

                if i==0 and sk>1:
                    print(" Calculating Brier Score for a large array, this may take some time ...")

                model = self.model[:,ind][:,::sk]
                obs = self.obs[0,ind][::sk]

                true_binary=[]; binary_forecasts=[]; fprob_forecasts=[]
                for j in range(0,obs.shape[0]):
                    atrue_binary = 1 if obs[j] >= threshold else 0
                    abinary_forecasts = [1 if auxv >= threshold else 0 for auxv in model[:,j]]
                    atrue_binary = [atrue_binary]*len(abinary_forecasts)        
                    true_binary=np.append(true_binary,atrue_binary)
                    binary_forecasts=np.append(binary_forecasts,abinary_forecasts)
                    fprob_forecasts=np.append(fprob_forecasts,prob_forecasts)
                    del atrue_binary,abinary_forecasts

                brier[i] = brier_score_loss(true_binary, binary_forecasts, sample_weight=fprob_forecasts)
                del model,obs,true_binary,binary_forecasts,fprob_forecasts

                if sk>1:
                    print("   - Brier Score, Ok: "+repr(self.fctxt[i])+"  of "+repr(self.fctxt.shape[0]))

            del ind

        # plot brier_score
        fig1 = plt.figure(1,figsize=(9,4)); ax = fig1.add_subplot(111)
        ax.plot(self.fctxt,brier,'k',marker='.',linewidth=2.,zorder=2)
        ax.set_xticks(self.fctxt)
        if self.fctunits:
            ax.set_xlabel('Forecast Lead Time ('+self.fctunits+')',size=sl)
        else:
            ax.set_xlabel('Forecast Lead Time ',size=sl)

        if self.vaxisname != "Model and Observations":
            ax.set_ylabel("Brier Score "+self.vaxisname)
        else:
            ax.set_ylabel("Brier Score")

        ax.grid(c='dimgrey', ls='--', alpha=0.3)
        plt.tight_layout();plt.axis('tight')

        if np.size(self.fctxticks)>1:
            ax.set_xticks(self.fctxticks)
            plt.xlim(xmin = np.min(self.fctxticks)*0.9, xmax = np.max(self.fctxticks)*1.01)

        plt.savefig(self.ftag+"BrierScoreXForecastTime.png", dpi=200, facecolor='w', edgecolor='w',orientation='portrait', format='png',transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close(fig1); del fig1, ax

        if result!="no":
            print(" Returning brier_score array, from ModelObsPlot.brierxfctime")
            return brier


#    def rocaucxfctime(self,result="no"):


# --------------------------------
# ------ Global Map --------------

def gsmooth(data,gflev=1):
    """
    Gaussian filter for spatial smoothness.
    data is a 2-dimensional field. data.shape must be equal to 2.
    gflev is the level of smoothing.    
    """
    from scipy.ndimage.filters import gaussian_filter
    data=np.array(data)
    gflev=float(gflev)

    fdata1=gaussian_filter(data,gflev); fdata2=gaussian_filter(data,gflev/2)
    fdata=np.nanmean(np.append(np.append(np.array([data]),np.array([fdata1]),axis=0),np.array([fdata2]),axis=0),axis=0)

    return fdata
    del fdata1,fdata

class GlobalMapPlot:
    def __init__(self,lat=None,lon=None,gdata=None,latmin=-85.,latmax=85.,
        ftag=None,levels=None,palette=None,pextend="max", contcolor="lightgrey",sl=13):

        self.lat = lat
        self.lon = lon
        self.latmin = latmin
        self.latmax = latmax
        self.gdata = gdata
        self.levels = levels
        self.pextend = pextend
        self.contcolor = contcolor
        self.sl = sl 

        self.ftag = ftag if ftag is not None else str(os.getcwd()) + '/GlobalMap.png'
        self.palette = palette if palette is not None else plt.cm.jet

        if self.gdata is None or self.lon is None or self.lat is None:
            raise ValueError("Missing required input arguments.")
            return

        matplotlib.rcParams.update({'font.size': sl}); plt.rc('font', size=sl)
        matplotlib.rc('xtick', labelsize=sl); matplotlib.rc('ytick', labelsize=sl); matplotlib.rcParams.update({'font.size': sl})

    def plot(self):

        if self.levels is None:
            if self.pextend == "both":
                aux=np.nanpercentile(np.abs(self.gdata), 99.)
                self.levels = np.linspace(-aux,aux, 101); del aux
            else:
                self.levels = np.linspace(0, np.nanpercentile(self.gdata, 99.), 101)       	

        plt.figure(figsize=(7, 4))
        # ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=-90))
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
        ax.set_extent([0., 360., self.latmin, self.latmax], crs=ccrs.PlateCarree())
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='grey', alpha=0.5, linestyle='--')
        gl.xlabel_style = {'size': self.sl-4, 'color': 'k', 'rotation': 0}
        gl.ylabel_style = {'size': self.sl-4, 'color': 'k', 'rotation': 0}
        cs = ax.contourf(self.lon, self.lat, self.gdata, self.levels, cmap=self.palette, zorder=1, extend=self.pextend, transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.OCEAN, facecolor=("white"))
        ax.add_feature(cfeature.LAND, facecolor=(self.contcolor), edgecolor='grey', linewidth=0.2, zorder=2)
        ax.add_feature(cfeature.BORDERS, edgecolor='dimgrey', linestyle='-', linewidth=0.2, alpha=1, zorder=3)
        ax.coastlines(resolution='110m', color='dimgrey', linewidth=0.2, linestyle='-', alpha=1, zorder=4)
        plt.tight_layout()
        ax = plt.gca()
        pos = ax.get_position()
        l, b, w, h = pos.bounds
        cax = plt.axes([l + 0.07, b - 0.075, w - 0.12, 0.025])
        cbar = plt.colorbar(cs, cax=cax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=self.sl-3)
        plt.axes(ax)
        plt.savefig(self.ftag, dpi=200, facecolor='w', edgecolor='w',
                    orientation='portrait', format='png', transparent=False, bbox_inches='tight', pad_inches=0.1)
        plt.close('all')
        del ax


