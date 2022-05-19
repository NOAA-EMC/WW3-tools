#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3colloc_buoy.py - collocation of WW3 model points and buoy stations.
"""
import sys
import glob
import subprocess
import numpy as np

def ww3_colloc_buoy(*args):
        '''
        PURPOSE: determine collocation of WW3 model output points with buoy stations.
        USAGE:   ww3_colloc_buoy <ww3-point-tab>
        RETURN:  model-buoy collocation file name
        '''
        if   len(args) == 0:
                sys.exit('ww3_colloc_buoy <ww3-point-tab>')
        elif len(args) == 1:
                ww3_tab=np.str(args[0])
        elif len(args)  > 1:
                sys.exit('ww3_colloc_buoy <ww3-point-tab>')

        print()
        print('ww3_colloc_buoy')
        print('---------------')
        print('  WW3 point tab file: ' + ww3_tab)

        listfile='./ww3list.txt'
        try:
                with open(listfile, 'w') as f: f.write(ww3_tab + '\n' + ww3_tab + '\n')
        except:
                sys.exit("unable to write list file: {}".format(listfile))
        try:
                subprocess.call("./modelBuoy_collocation_hindcast.py", shell=True)
        except:
                sys.exit("problem calling: modelBuoy_collocation_hindcast.py")
        try:
                colloc_files=glob.glob('./ww3buoy_collocation*.nc')
        except:
                sys.exit("unable to locate model-buoy collocation file.")
        if len(colloc_files) > 1:
                print('more than 1 model-buoy collocation file found.')
                print(colloc_files)
                sys.exit('aborting.')
        colloc_file = colloc_files[0]
        print('  Model-buoy collocation file: ' + colloc_file)
        return colloc_file
