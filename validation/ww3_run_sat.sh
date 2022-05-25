#!/bin/bash

# ww3_run_sat.sh
# PURPOSE: shell script wrapper around python satellite collocation 
#          and generation of statistics and metrics.
#
# USAGE:   ww3_run_sat.sh <run_id> <ww3.grd.nc> [<var> [<var_min> [<var_max>]]]
#
# OUTPUT:  
#
export PATH=/scratch2/NCEPDEV/marine/Matthew.Masarik/wavpy/miniconda3/bin:${PATH}
export PYTHONPATH=/scratch2/NCEPDEV/marine/Matthew.Masarik/wavpy/miniconda3/pkgs:${PYTHONPATH}

usage_str='ww3_run_sat.sh <run_id> <ww3.grd.nc> [<var> [<var_min> [<var_max>]]]'
if [[ $# -lt 2 ]] || [[ $# -gt 5 ]]; then echo $usage_str; exit; fi

run_id="$1"
ww3_grd="$2"
shift; shift
if [ ! -f $ww3_grd ]; then
    echo "file: $ww3_grd, not found"
    exit
fi
echo -e "\nInput file: $ww3_grd"

# existing sat-model collocation file, ww3list.txt, satlist.txt
num_sat_colloc=$(ls -1 WW3.Altimeter_*.nc 2> /dev/null | wc -l)
if [ $num_sat_colloc -ne 0 ]; then
    echo "Found ${num_sat_colloc} model-sat collocation file."
    echo "Remove all files before re-running."
    echo "Exiting."
    exit
fi
if [ -f ww3list.txt ]; then rm -f ww3list.txt; fi
if [ -f satlist.txt ]; then rm -f satlist.txt; fi

# existing statistics and metrics files
if [ -f sat_stats_obs*.txt   ]; then rm -f sat_stats_obs*.txt;   fi
if [ -f sat_stats_model*.txt ]; then rm -f sat_stats_model*.txt; fi
if [ -f sat_metrics*.txt     ]; then rm -f sat_metrics*.txt;     fi

# the call
python3 run_sat.py $ww3_grd $@

# data mgmt
sat_colloc=$(ls -1 WW3.Altimeter_*.nc 2> /dev/null)
if [ ! -f $sat_colloc ]; then
    echo "Model-sat collocation file not found. Exiting." 
    exit
fi
run_sat_colloc=${run_id}.${sat_colloc}
run_ww3list=${run_id}.ww3list.txt
run_satlist=${run_id}.satlist.txt
mv $sat_colloc   $run_sat_colloc
mv ww3list.txt   $run_ww3list
mv satlist.txt   $run_satlist

if [ ! -f sat_stats_obs*.txt ]; then 
    echo 'sat obs stats not found.'
else
    sat_stats_obs=$(ls -1 sat_stats_obs*.txt 2> /dev/null)
    run_stats_obs=${run_id}.${sat_stats_obs}
    mv $sat_stats_obs   $run_stats_obs
fi

if [ ! -f sat_stats_model*.txt ]; then 
    echo 'sat model stats not found.'
else
    sat_stats_model=$(ls -1 sat_stats_model*.txt 2> /dev/null)
    run_stats_model=${run_id}.${sat_stats_model}
    mv $sat_stats_model  $run_stats_model
fi

if [ ! -f sat_metrics*.txt ]; then 
    echo 'sat metrics not found.'
else
    mod_sat_mets=$(ls -1 sat_metrics*.txt 2> /dev/null)
    run_mets=${run_id}.${mod_sat_mets}
    mv $mod_sat_mets   $run_mets
fi

echo -e "\nInput WW3 grd:      $ww3_grd\n"
run_out_dir=run_out_$run_id
mkdir $run_out_dir
if [ -d $run_out_dir ]; then
    mv -f $run_ww3list       $run_out_dir
    mv -f $run_satlist       $run_out_dir
    mv -f $run_sat_colloc    $run_out_dir
    mv -f $run_stats_obs     $run_out_dir
    mv -f $run_stats_model   $run_out_dir
    mv -f $run_mets          $run_out_dir
    echo -e "Output Directory:           ${run_out_dir}/"
fi
echo -e "\tww3list:            $run_ww3list"
echo -e "\tsatlist:            $run_satlist"
echo -e "\tSat collocation:    $run_sat_colloc"
echo -e "\tSat stats:          $run_stats_obs"
echo -e "\tModel stats:        $run_stats_model"
echo -e "\tValidation metrics: $run_mets"
