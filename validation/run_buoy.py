#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
run_buoy.py - driver for 1. buoy collocation, 2. statistics, and 3. validation metric functions.
"""

import sys
import numpy as np
from ww3stats_buoy   import ww3_stats_buoy
from ww3metrics_buoy import ww3_metrics_buoy
from ww3colloc_buoy  import ww3_colloc_buoy


def main(argv):
    '''
    PURPOSE: Main for calling buoy collocation, statistics, and metrics functions.
    USAGE:   python3 run_buoy.py <ww3-point-tab.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]
    '''
    args=sys.argv[1:]
    vmin=0; vmax=20; parm='hs'
    if   len(args) == 0:
        sys.exit('python3 run_buoy.py <ww3-point-tab.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')
    elif len(args) == 1:
        ww3_tab=np.str(args[0])
    elif len(args) == 2:
        ww3_tab=np.str(args[0]); parm=np.str(args[1])
    elif len(args) == 3:
        ww3_tab=np.str(args[0]); parm=np.str(args[1]); vmin=np.int(args[2])
    elif len(args) == 4:
        ww3_tab=np.str(args[0]); parm=np.str(args[1]); vmin=np.int(args[2]); vmax=np.int(args[3])
    elif len(args)  > 4:
        sys.exit('python3 run_buoy.py <ww3-point-tab.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')


    print('\n*  RUN BUOY  *\n')

    print('\nCalling ww3_colloc_buoy...')
    mod_buoy_col = ww3_colloc_buoy(ww3_tab)

    print('\nCalling ww3_stats_buoy...')
    ww3_stats_buoy(mod_buoy_col, parm, vmin, vmax)

    print('\nCalling ww3_metrics_buoy...')
    ww3_metrics_buoy(mod_buoy_col, parm, vmin, vmax)

    print('\nRUN BUOY Complete.')



if __name__ == "__main__":
    main(sys.argv[1:])
