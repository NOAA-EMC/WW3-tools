#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_sat.py - driver for 1. sat collocation, 2. statistics, and 3. validation metric functions.
"""
import sys
import numpy as np
from ww3stats_sat   import ww3_stats_sat
from ww3metrics_sat import ww3_metrics_sat
from ww3colloc_sat  import ww3_colloc_sat

# Satellite altimeters
sat_alt1='AltimeterGridded_CRYOSAT2.nc'
sat_alt2='AltimeterGridded_JASON3.nc'
sat_alt3='AltimeterGridded_SARAL.nc'
sat_alt4='AltimeterGridded_SENTINEL3A.nc'

def main(argv):
    '''
    PURPOSE: Main for calling sat collocation, statistics, and metrics functions.
    USAGE:   python3 run_sat.py <ww3-grd.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]
    '''
    args=sys.argv[1:]
    vmin=0; vmax=20; parm='hs'
    if   len(args) == 0:
        sys.exit('python3 run_sat.py <ww3-grd.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')
    elif len(args) == 1:
        ww3_grd=str(args[0])
    elif len(args) == 2:
        ww3_grd=str(args[0]); parm=str(args[1])
    elif len(args) == 3:
        ww3_grd=str(args[0]); parm=str(args[1]); vmin=np.int(args[2])
    elif len(args) == 4:
        ww3_grd=str(args[0]); parm=str(args[1]); vmin=np.int(args[2]); vmax=np.int(args[3])
    elif len(args)  > 4:
        sys.exit('python3 run_sat.py <ww3-grd.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')


    print('\n*  RUN SAT  *\n')

    print('\nCalling ww3_colloc_sat...')
    mod_sat_col = ww3_colloc_sat(ww3_grd, sat_alt1, sat_alt2, sat_alt3, sat_alt4)

    print('\nCalling ww3_stats_sat...')
    ww3_stats_sat(mod_sat_col, parm, vmin, vmax)

    print('\nCalling ww3_metrics_sat...')
    ww3_metrics_sat(mod_sat_col, parm, vmin, vmax)

    print('\nRUN SAT Complete.')



if __name__ == "__main__":
    main(sys.argv[1:])
