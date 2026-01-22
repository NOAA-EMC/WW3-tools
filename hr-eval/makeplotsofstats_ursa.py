import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import argparse

import os
import datetime as dt
from dateutil.relativedelta import relativedelta

from emcpy.plots.plots import Scatter, LinePlot
from emcpy.plots.create_plots import CreatePlot, CreateFigure




'''
Read stats and plot 
'''

def main():


  datapath = "/scratch4/NCEPDEV/marine/Ming.Chen/wave_eval/HR_eval/hr-eval/myfile.nc"

  datanc  = nc.Dataset(datapath)

  allstats_hs = np.array(datanc.variables['allstats_hs'][:])
  allstats_wnd = np.array(datanc.variables['allstats_wnd'][:])
  allstats_hs_cal = np.array(datanc.variables['allstats_hs_cal'][:])
  allstats_wnd_cal = np.array(datanc.variables['allstats_wnd_cal'][:])

  #allstats_wnd  _cal
  datanc.close()

    #seasons[k],satelites,model, (all, above 4hs, above 7hs),day, stats:bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N

  season=['winter']
  satelites=['JASON3', 'CRYOSAT2', 'SARAL', 'SENTINEL3A'] #JASON3,JASON2,CRYOSAT2,JASON1,HY2,SARAL,SENTINEL3A,ENVISAT,ERS1,ERS2,GEOSAT,GFO,TOPEX,SENTINEL3B,CFOSAT
  colorsformodel=['darkred','blue','darkgreen','darkorange','deeppink','purple']
  stats=['bias', 'RMSE','NBias', 'NRMSE', 'SCrmse', 'SI', 'HH','CC', 'N']
  for k in range(len(season)):
    if season[k] == "winter":
       startdate = dt.datetime(2024,12,1,0)
       enddate = dt.datetime(2024,12,17,18)
       datestride = 6
       endday = 16
       #model=['Retrotests', 'GFSv16']
       model=['GFSv16','Retrotests']
       orderofplots=[0,1,2,3,4]
    elif season[k] == "summer":
       startdate = dt.datetime(2020,6,1)
       enddate = dt.datetime(2020,8,30)
       datestride = 3
       endday = 16
       model=['multi1','HR1', 'HR2', 'HR3a', 'HR3b', 'GFSv16']
       orderofplots=[0,5,1,2,3,4]
    elif season[k] == "hurricane":
       startdate = dt.datetime(2020,7,20)
       enddate = dt.datetime(2020,11,20)
       datestride = 1
       endday = 7
       model=['multi1', 'HR1', 'HR2', 'HR3a', 'HR3b', 'GFSv16']
       orderofplots=[0,5,1,2,3,4]
    xday=np.arange(1, endday+1, 1)*24
    yday=np.arange(0, endday, 1)
    xticksday=np.arange(1, 16+1, 1)*24
    for j in range(len(satelites)):
     for s in range(len(stats)):
      plot1 = CreatePlot()  # Create  Plot
      plt_list = []  # initialize emtpy plot list

      # Bottom plot
      plot2 = CreatePlot()  # Create Plot
      plt_list2 = []  # initialize empty plot list


      for mm in range(len(model)):
        m=orderofplots[mm]
        print(season[k])
        print(model[m])
        #seasons[k],satelites,model, (all, above 4hs, above 7hs),day, stats:bias, RMSE, NBias, NRMSE, SCrmse, SI, HH, CC, N
        print(satelites[j]) 
        print(allstats_hs[k,j,m,0,:,1]) 
  

        # Top (HS) plot
        lp = LinePlot(xday, allstats_hs[k,j,m,0,yday,s])  # Create line plot object
        lp.color = colorsformodel[m]   # line color
        lp.linestyle = "-"  # line style
        lp.linewidth = 1.5  # line width
        lp.marker = "o"  # marker type
        lp.markersize = 4  # markersize
        lp.alpha = None  # transparency
        lp.label = f"{model[m]} (all)"  # give it a label
        plt_list.append(lp)  # Add line plot object to list

        lp = LinePlot(xday, allstats_hs[k,j,m,1,yday,s])  # Create line plot object
        lp.color = colorsformodel[m]   # line color
        lp.linestyle = "--"  # line style
        lp.linewidth = 1.5  # line width
        lp.marker = "o"  # marker type
        lp.markersize = 4  # markersize
        lp.alpha = None  # transparency
        lp.label = f"{model[m]} (Hs>4m)"  # give it a label
        plt_list.append(lp)  # Add line plot object to list

        #lp = LinePlot(xday, allstats_hs[k,j,m,2,yday,s])  # Create line plot object
        #lp.color = colorsformodel[m]   # line color
        #lp.linestyle = "-."  # line style
        #lp.linewidth = 1.5  # line width
        ##lp.marker = "o"  # marker type
        ##lp.markersize = 4  # markersize
        #lp.alpha = None  # transparency
        #lp.label = "line2"  # give it a label
        #plt_list.append(lp)  # Add line plot object to list


        # Bottom plot
        lp = LinePlot(xday, allstats_wnd_cal[k,j,m,0,yday,s])  # Create line plot object
        lp.color = colorsformodel[m]  # line color
        lp.linestyle = "-"  # line style
        lp.linewidth = 1.5  # line width
        lp.marker = "o"  # marker type
        lp.markersize = 4  # markersize
        lp.alpha = None  # transparency
        lp.label = f"{model[m]}"  # give it a label
        plt_list2.append(lp)  # Add line plot object to list

        #lp = LinePlot(xday, allstats_wnd[k,j,m,1,yday,s])  # Create line plot object
        #lp.color = colorsformodel[m]  # line color
        #lp.linestyle = "--"  # line style
        #lp.linewidth = 1.5  # line width
        ##lp.marker = "o"  # marker type
        ##lp.markersize = 4  # markersize
        #lp.alpha = None  # transparency
        #lp.label = model[m]  # give it a label
        #plt_list2.append(lp)  # Add line plot object to list


      plot1.plot_layers = plt_list  # draw plot1 (the top plot)
      plot2.plot_layers = plt_list2  # draw plot2 (the bottom plot)

      # Add plot features
      plot1.add_title(label="Significant Wave Height")
      plot1.add_xlabel(xlabel="Forecast Hour")
      plot1.add_ylabel(ylabel=stats[s])
      plot1.add_grid()
      plot1.set_xticks(xticksday)
      plot1.set_xlim(0,384) #endday*24)
      #stats=['bias', 'RMSE','NBias', 'NRMSE', 'SCrmse', 'SI', 'HH','CC', 'N']
      #         0        1     2        3        4        5     6    7    8
      if s == 0:
        plot1.set_ylim(-1.75,1.75)
      elif s == 1: 
        plot1.set_ylim(0,3)
      elif s==7:
        plot1.set_ylim(0,1)
      #plot1.set_xticklabels([str(item) for item in x1], rotation=0)
      #yticks = np.arange(np.min(y2), np.max(y2) + 1, 1)
      #plot1.set_yticks(yticks)
      #plot1.set_yticklabels([str(item) for item in yticks], rotation=0)
      plot1.add_legend(loc="lower left", fancybox=True, framealpha=0.80)
      ##plot1.add_legend(loc="upper left", fancybox=True, framealpha=1)
      # Add plot features
      plot2.add_title(label="Wind Speed")
      plot2.add_xlabel(xlabel="Forecast Hour")
      plot2.add_ylabel(ylabel=stats[s])
      plot2.add_grid()
      plot2.set_xticks(xticksday)#xday)
      plot2.set_xlim(0,384)
      if s == 0:
        plot2.set_ylim(-1,0)
      elif s == 1:
        plot2.set_ylim(0,5)
      elif s==7:
        plot2.set_ylim(0,1)
      #plot2.set_xticklabels([str(item) for item in x2], rotation=0)
      #yticks = np.arange(np.min(y2), np.max(y2) + 1, 1)
      #plot2.set_yticks(yticks)
      #plot2.set_yticklabels([str(item) for item in yticks], rotation=0)
      plot2.add_legend(loc="lower left", fancybox=True, framealpha=0.80)
      #plot2.add_legend(loc="upper left", fancybox=True, framealpha=1)


      # Return matplotlib figure
      fig = CreateFigure(nrows=2, ncols=1, figsize=(8, 6))
      fig.plot_list = [plot1, plot2]
      fig.create_figure()
      supertitle=f"fig_{stats[s]}_{satelites[j]}_{season[k]}"
      fig.add_suptitle(supertitle)
      fig.tight_layout()
      #plotpath = outpath+'/scatter'+'%03d.png' % (ihr)
      #fig.save_figure(plotpath)
      plotpath=f"fig_{stats[s]}_{satelites[j]}_{season[k]}.png"
      fig.save_figure(plotpath)

      fig.close_figure()




if __name__ == '__main__':
    main()

