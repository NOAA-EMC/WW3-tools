#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3stats_buoy.py - calculation of summary statistics for WW3 point output and buoy obs.
"""
import sys
import netCDF4 as nc
import numpy as np
import merrevalstats as mer


def ww3_stats_buoy(*args):
        '''
        PURPOSE: calculate summary statistics from buoy collocation input.
        USAGE:   ww3_stats_buoy <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]
        '''
        vmin=0; vmax=20; parm='hs'
        if   len(args) == 0:
                sys.exit('ww3_stats_buoy <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')
        elif len(args) == 1:
                dat=np.str(args[0])
        elif len(args) == 2:
                dat=np.str(args[0]); parm=np.str(args[1])
        elif len(args) == 3:
                dat=np.str(args[0]); parm=np.str(args[1]); vmin=np.int(args[2])
        elif len(args) == 4:
                dat=np.str(args[0]); parm=np.str(args[1]); vmin=np.int(args[2]); vmax=np.int(args[3])
        elif len(args)  > 4:
                sys.exit('ww3_stats_buoy <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')

        print()
        print('ww3_stats_buoy')
        print('--------------')
        print('  Buoy collocation file: ' + dat)
        print('  WW3 parameter: ' + parm)
        print('  Parameter Min: ' + str(vmin))
        print('  Parameter Max: ' + str(vmax))
        print('  Calculating...')
        mod_parm  = 'model_' + parm
        buoy_parm = 'obs_'   + parm

        f      =nc.Dataset(dat)
        m_vals =f.variables[mod_parm]
        b_vals =f.variables[buoy_parm]
        M_VALS =np.reshape(m_vals,(m_vals.shape[0]*m_vals.shape[1]))
        B_VALS =np.reshape(b_vals,(b_vals.shape[0]*b_vals.shape[1]))
        ind    =np.where((M_VALS>0 )&(B_VALS>0))

        m_stats=np.array([mer.smrstat(M_VALS[ind],vmin,vmax)])
        b_stats=np.array([mer.smrstat(B_VALS[ind],vmin,vmax)])

        m_out  ='buoy_stats_model_' + parm + '.txt'
        b_out  ='buoy_stats_obs_'   + parm + '.txt'
        np.savetxt(m_out, m_stats, fmt='%10.4f', delimiter=' ', newline='\n', footer='', comments='')
        np.savetxt(b_out, b_stats, fmt='%10.4f', delimiter=' ', newline='\n', footer='', comments='')
        print('  Model output: ' + m_out)
        print('  Buoy output:  ' + b_out + '\n')
