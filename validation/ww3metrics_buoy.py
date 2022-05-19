#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3metrics_buoy.py - calculation of validation metrics for WW3 point output and buoy obs.
"""
import sys
import netCDF4 as nc
import numpy as np
import merrevalstats as mer


def ww3_metrics_buoy(*args):
        '''
        PURPOSE: calculate validation metrics from buoy collocation input.
        USAGE:  ww3_metrics_buoy <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]
        '''
        vmin=0; vmax=20; parm='hs'
        if   len(args) == 0:
                sys.exit('ww3_metrics_buoy <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')
        elif len(args) == 1:
                dat=np.str(args[0])
        elif len(args) == 2:
                dat=np.str(args[0]); parm=np.str(args[1])
        elif len(args) == 3:
                dat=np.str(args[0]); parm=np.str(args[1]); vmin=np.int(args[2])
        elif len(args) == 4:
                dat=np.str(args[0]); parm=np.str(args[1]); vmin=np.int(args[2]); vmax=np.int(args[3])
        elif len(args)  > 4:
                sys.exit('ww3_metrics_buoy <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')

        print()
        print('ww3_metrics_buoy')
        print('----------------')
        print('  Buoy collocation file: ' + dat)
        print('  WW3 parameter: ' + parm)
        print('  Parameter Min: ' + str(vmin))
        print('  Parameter Max: ' + str(vmax))
        print('  Calculating...')
        mod_parm  = 'model_' + parm
        buoy_parm = 'obs_'   + parm

        f       =nc.Dataset(dat)
        m_vals  =f.variables[mod_parm]
        b_vals  =f.variables[buoy_parm]
        M_VALS  =np.reshape(m_vals,(m_vals.shape[0]*m_vals.shape[1]))
        B_VALS  =np.reshape(b_vals,(b_vals.shape[0]*b_vals.shape[1]))
        ind     =np.where((M_VALS>0 )&(B_VALS>0))
        mets    =np.array([mer.metrics(M_VALS[ind],B_VALS[ind],vmin,vmax)])

        met_out ='buoy_metrics_' + parm + '.txt'
        np.savetxt(met_out, mets, fmt='%10.4f', delimiter=' ', newline='\n', footer='', comments='')
        print('  Metric output:  ' + met_out + '\n')
