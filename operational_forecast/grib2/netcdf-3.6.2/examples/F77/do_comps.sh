#!/bin/sh
# This shell script runs the cmp test on the example programs.
# $Id: do_comps.sh,v 1.3 2006/06/14 20:48:36 ed Exp $

set -e
echo ""
echo "*** Testing that F77 examples produced same files as C examples."
echo "*** checking simple_xy.nc..."
cmp simple_xy.nc ../C/simple_xy.nc

echo "*** checking sfc_pres_temp.nc..."
cmp sfc_pres_temp.nc ../C/sfc_pres_temp.nc

echo "*** checking pres_temp_4D.nc..."
cmp pres_temp_4D.nc ../C/pres_temp_4D.nc

echo "*** All F77 example comparisons worked!"
exit 0
