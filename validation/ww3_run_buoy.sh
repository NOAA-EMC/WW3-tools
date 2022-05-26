#!/bin/bash

# ww3_run_buoy.sh
# PURPOSE: shell script wrapper around python buoy collocation 
#          and generation of statistics and metrics.
#
# USAGE:   ww3_run_buoy.sh <run_id> <ww3.tab.nc> [<var> [<var_min> [<var_max>]]]
#
# OUTPUT:  
#
export PATH=/scratch2/NCEPDEV/marine/Matthew.Masarik/wavpy/miniconda3/bin:${PATH}
export PYTHONPATH=/scratch2/NCEPDEV/marine/Matthew.Masarik/wavpy/miniconda3/pkgs:${PYTHONPATH}

usage_str='ww3_run_buoy.sh <run_id> <ww3.tab.nc> [<var> [<var_min> [<var_max>]]]'
if [[ $# -lt 2 ]] || [[ $# -gt 5 ]]; then echo $usage_str; exit; fi

run_id="$1"
ww3_tab="$2"
shift; shift
if [ ! -f $ww3_tab ]; then
    echo "file: $ww3_tab, not found"
    exit
fi
echo -e "\nInput file: $ww3_tab"

# existing buoy-model collocation file, ww3list.txt
num_buoy_colloc=$(ls -1 WW3.Buoy_*.nc 2> /dev/null | wc -l)
if [ $num_buoy_colloc -ne 0 ]; then
    echo "Found ${num_buoy_colloc} model-buoy collocation file."
    echo "Remove all files before re-running."
    echo "Exiting."
    exit
fi
if [ -f ww3list.txt ]; then rm -f ww3list.txt; fi

# existing statistics and metrics files
if [ -f buoy_stats_obs*.txt   ]; then rm -f buoy_stats_obs*.txt;   fi
if [ -f buoy_stats_model*.txt ]; then rm -f buoy_stats_model*.txt; fi
if [ -f buoy_metrics*.txt     ]; then rm -f buoy_metrics*.txt;     fi

# the call
python3 run_buoy.py $ww3_tab $@

# data mgmt
buoy_colloc=$(ls -1 WW3.Buoy_*.nc 2> /dev/null)
if [ ! -f $buoy_colloc ]; then
    echo "Model-buoy collocation file not found. Exiting." 
    exit
fi
run_buoy_colloc=${run_id}.${buoy_colloc}
run_ww3list=${run_id}.ww3list.txt
mv $buoy_colloc  $run_buoy_colloc
mv ww3list.txt   $run_ww3list

if [ ! -f buoy_stats_obs*.txt ]; then 
    echo 'buoy obs stats not found.'
else
    buoy_stats_obs=$(ls -1 buoy_stats_obs*.txt 2> /dev/null)
    run_stats_obs=${run_id}.${buoy_stats_obs}
    mv $buoy_stats_obs   $run_stats_obs
fi

if [ ! -f buoy_stats_model*.txt ]; then 
    echo 'buoy model stats not found.'
else
    buoy_stats_model=$(ls -1 buoy_stats_model*.txt 2> /dev/null)
    run_stats_model=${run_id}.${buoy_stats_model}
    mv $buoy_stats_model  $run_stats_model
fi

if [ ! -f buoy_metrics*.txt ]; then 
    echo 'buoy metrics not found.'
else
    mod_buoy_mets=$(ls -1 buoy_metrics*.txt 2> /dev/null)
    run_mets=${run_id}.${mod_buoy_mets}
    mv $mod_buoy_mets   $run_mets
fi

echo -e "\nInput WW3 tab:      $ww3_tab\n"
run_out_dir=run_out_$run_id
mkdir $run_out_dir
if [ -d $run_out_dir ]; then
    mv -f $run_ww3list       $run_out_dir
    mv -f $run_buoy_colloc   $run_out_dir
    mv -f $run_stats_obs     $run_out_dir
    mv -f $run_stats_model   $run_out_dir
    mv -f $run_mets          $run_out_dir
    echo -e "Output Directory:           ${run_out_dir}/"
fi
echo -e "\tww3list:            $run_ww3list"
echo -e "\tBuoy collocation:   $run_buoy_colloc"
echo -e "\tBuoy stats:         $run_stats_obs"
echo -e "\tModel stats:        $run_stats_model"
echo -e "\tValidation metrics: $run_mets"
