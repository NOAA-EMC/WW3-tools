#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3stats_sat.py - calculation of summary statistics for WW3 point output and satellite obs.
"""
import sys
import netCDF4 as nc
import numpy as np
import merrevalstats as mer


def ww3_stats_sat(*args):
        '''
        PURPOSE: calculate summary statistics from satellite collocation input.
        USAGE:   ww3_stats_sat <mod_sat_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]
        '''
        vmin=0; vmax=20; parm='hs'
        if   len(args) == 0:
                sys.exit('ww3_stats_sat <mod_sat_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')
        elif len(args) == 1:
                dat=str(args[0])
        elif len(args) == 2:
                dat=str(args[0]); parm=str(args[1])
        elif len(args) == 3:
                dat=str(args[0]); parm=str(args[1]); vmin=np.int(args[2])
        elif len(args) == 4:
                dat=str(args[0]); parm=str(args[1]); vmin=np.int(args[2]); vmax=np.int(args[3])
        elif len(args)  > 4:
                sys.exit('ww3_stats_sat <mod_sat_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')

        print()
        print('ww3_stats_sat')
        print('-------------')
        print('  Sat collocation file: ' + dat)
        print('  WW3 parameter: ' + parm)
        print('  Parameter Min: ' + str(vmin))
        print('  Parameter Max: ' + str(vmax))
        print('  Calculating...')
        mod_parm  = 'model_' + parm
        sat_parm  = 'obs_'   + parm

        f      =nc.Dataset(dat)
        m_vals =f.variables[mod_parm]
        s_vals =f.variables[sat_parm]
        M_VALS =np.array(m_vals[:])
        S_VALS =np.array(s_vals[:])
        ind    =np.where((M_VALS>0.)&(S_VALS>0.))

        m_stats=np.array([mer.smrstat(M_VALS[ind],vmin,vmax)])
        s_stats=np.array([mer.smrstat(S_VALS[ind],vmin,vmax)])

        m_out  ='sat_stats_model_' + parm + '.txt'
        s_out  ='sat_stats_obs_'   + parm + '.txt'
        np.savetxt(m_out, m_stats, fmt='%10.4f', delimiter=' ', newline='\n', footer='', comments='')
        np.savetxt(s_out, s_stats, fmt='%10.4f', delimiter=' ', newline='\n', footer='', comments='')
        print('  Model output: ' + m_out)
        print('  Sat output:   ' + s_out + '\n')
