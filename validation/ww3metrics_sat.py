#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3metrics_sat.py - calculation of validation metrics for WW3 point output and satellite obs.
"""
import sys
import netCDF4 as nc
import numpy as np
import merrevalstats as mer


def ww3_metrics_sat(*args):
        '''
        PURPOSE: calculate validation metrics from satellite collocation input.
        USAGE:  ww3_metrics_sat <mod_sat_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]
        '''
        vmin=0; vmax=20; parm='hs'
        if   len(args) == 0:
                sys.exit('ww3_metrics_sat <mod_sat_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')
        elif len(args) == 1:
                dat=str(args[0])
        elif len(args) == 2:
                dat=str(args[0]); parm=str(args[1])
        elif len(args) == 3:
                dat=str(args[0]); parm=str(args[1]); vmin=np.int(args[2])
        elif len(args) == 4:
                dat=str(args[0]); parm=str(args[1]); vmin=np.int(args[2]); vmax=np.int(args[3])
        elif len(args)  > 4:
                sys.exit('ww3_metrics_sat <mod_sat_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')

        print()
        print('ww3_metrics_sat')
        print('---------------')
        print('  Sat collocation file: ' + dat)
        print('  WW3 parameter: ' + parm)
        print('  Parameter Min: ' + str(vmin))
        print('  Parameter Max: ' + str(vmax))
        print('  Calculating...')
        mod_parm  = 'model_' + parm
        sat_parm  = 'obs_'   + parm

        f       =nc.Dataset(dat)
        m_vals  =f.variables[mod_parm]
        s_vals  =f.variables[sat_parm]
        M_VALS  =np.array(m_vals[:])
        S_VALS  =np.array(s_vals[:])
        ind     =np.where((M_VALS>0 )&(S_VALS>0))
        mets    =np.array([mer.metrics(M_VALS[ind],S_VALS[ind],vmin,vmax)])

        met_out ='sat_metrics_' + parm + '.txt'
        np.savetxt(met_out, mets, fmt='%10.4f', delimiter=' ', newline='\n', footer='', comments='')
        print('  Metric output:  ' + met_out + '\n')
