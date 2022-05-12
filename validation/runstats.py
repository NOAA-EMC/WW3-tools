#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
runstats.py - sample driver for statistics and validation metric functions.
"""

import sys
import copy
from ww3stats   import ww3_pnt_stats
from ww3metrics import ww3_pnt_metrics


def main(argv):
    '''
    PURPOSE: Main for calling statistics and metrics functions.
    USAGE:   python3 runstats.py <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]
    '''
    args=sys.argv[1:]
    vmin=0; vmax=20; parm='hs'
    if   len(args) == 0:
        sys.exit('python3 runstats.py <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')
    elif len(args) == 1:
        dat=copy.copy(args[0])
    elif len(args) == 2:
        dat=copy.copy(args[0]); parm=copy.copy(args[1])
    elif len(args) == 3:
        dat=copy.copy(args[0]); parm=copy.copy(args[1]); vmin=copy.copy(args[2])
    elif len(args) == 4:
        dat=copy.copy(args[0]); parm=copy.copy(args[1]); vmin=copy.copy(args[2]); vmax=copy.copy(args[3])
    elif len(args)  > 4:
        sys.exit('python3 runstats.py <mod_buoy_colloc.nc> [ <ww3_param> [ <param_min> [ <param_max> ]]]')


    print('\n*  RUNSTATS  *\n')

    print('\nCalling ww3_pnt_stats...')
    ww3_pnt_stats(dat, parm, vmin, vmax)

    print('\nCalling ww3_pnt_metrics...')
    ww3_pnt_metrics(dat, parm, vmin, vmax)

    print('\nRUNSTATS Complete.')



if __name__ == "__main__":
    main(sys.argv[1:])
