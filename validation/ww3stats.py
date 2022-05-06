#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3stats.py - calculation of summary statistics for WW3 point output.
"""
import sys
import copy
import netCDF4 as nc
import numpy as np
import merrevalstats as mer


def ww3_pnt_stats(*args):
        '''
        PURPOSE: calculate summary statistics from buoy collocation input.
        USAGE:   ww3_pnt_stats <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]
        '''
        vmin=-np.inf; vmax=np.inf; parm='hs'
        if   len(args) == 0:
                sys.exit('ww3_pnt_stats <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')
        elif len(args) == 1:
                dat=copy.copy(args[0])
        elif len(args) == 2:
                dat=copy.copy(args[0]); parm=copy.copy(args[1])
        elif len(args) == 3:
                dat=copy.copy(args[0]); parm=copy.copy(args[1]); vmin=copy.copy(args[2])
        elif len(args) == 4:
                dat=copy.copy(args[0]); parm=copy.copy(args[1]); vmin=copy.copy(args[2]); vmax=copy.copy(args[3])
        elif len(args)  > 4:
                sys.exit('ww3_pnt_stats <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')

        print()
        print('ww3_pnt_stats')
        print('-------------')
        print('  Buoy collocation file: ' + dat)
        print('  WW3 parameter: ' + parm)
        print('  Parameter Min: ' + str(vmin))
        print('  Parameter Max: ' + str(vmax))
        print('  Calculating...')
        mod_parm  = 'model_' + parm
        buoy_parm = 'buoy_'  + parm

        f      =nc.Dataset(dat)
        m_vals =f.variables[mod_parm]
        b_vals =f.variables[buoy_parm]
        M_VALS =np.reshape(m_vals,(m_vals.shape[0]*m_vals.shape[1]))
        B_VALS =np.reshape(b_vals,(b_vals.shape[0]*b_vals.shape[1]))
        ind    =np.where((M_VALS>0 )&(B_VALS>0))

        m_stats=np.array([mer.smrstat(M_VALS[ind],vmin,vmax)])
        b_stats=np.array([mer.smrstat(B_VALS[ind],vmin,vmax)])

        m_out  ='stats_model_' + parm + '.txt'
        b_out  ='stats_buoy_' + parm + '.txt'
        np.savetxt(m_out, m_stats, fmt='%10.4f', delimiter=' ', newline='\n', footer='', comments='')
        np.savetxt(b_out, b_stats, fmt='%10.4f', delimiter=' ', newline='\n', footer='', comments='')
        print('  Model output: ' + m_out)
        print('  Buoy output:  ' + b_out + '\n')
