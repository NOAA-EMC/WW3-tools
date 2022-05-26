#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ww3colloc_sat.py - collocation of WW3 model field output and satellite tracks
"""
import sys
import glob
import os
from   os import remove
from   os.path import exists as file_exists
import subprocess


def ww3_colloc_sat(ww3grd, *args):
        '''
        PURPOSE: determine collocation of WW3 model field output with satellites
        USAGE:   ww3_colloc_sat <ww3-grd> <sat_alt1> [<sat_alt2> ... <sat_altN>]
        RETURN:  model-sat collocation file name
        '''
        if   len(args) < 2:
                sys.exit('ww3_colloc_sat <ww3-grd> <sat_alt1> [<sat_alt2> ... <sat_altN>]')
        elif len(args) >= 2:
                ww3_grd=str(ww3grd)


        # input lists: satlist.txt, ww3list.txt
        satlist='satlist.txt'
        ww3list='ww3list.txt'
        if file_exists(satlist):
                try:
                         os.remove(satlist)
                except:
                         sys.exit("unable to remove existing sat list file: {}. aborting.".format(satlist))
        if file_exists(ww3list):
                try:
                         os.remove(ww3list)
                except:
                         sys.exit("unable to remove existing ww3 list file: {}. aborting.".format(ww3list))
                         

        print()
        print('ww3_colloc_sat')
        print('--------------')
        print('  WW3 field output file:  ' + ww3_grd)
        for arg in args:
                print('  Sat altimeter file:     ' + arg)
                try:
                         with open(satlist, 'a') as f: f.write(arg + '\n')
                except:
                         sys.exit("unable to write to altimeter file, {}, to sat list file: {}".format(arg,satlist))
                         
        try:
                with open(ww3list, 'w') as f: f.write(ww3_grd + '\n')
        except:
                sys.exit("unable to write list file: {}".format(ww3list))
        try:
                subprocess.call('./modelSat_collocation.py', shell=True)
        except:
                sys.exit('problem calling: modelSat_collocation.py')
        try:
                colloc_files=glob.glob('WW3.Altimeter_*.nc')
        except:
                sys.exit("unable to locate model-sat collocation file.")
        if len(colloc_files) > 1:
                print('more than 1 model-sat collocation file found.')
                print(colloc_files)
                sys.exit('aborting.')
        colloc_file = colloc_files[0]
        print('  Model-sat collocation file: ' + colloc_file)
        return colloc_file
