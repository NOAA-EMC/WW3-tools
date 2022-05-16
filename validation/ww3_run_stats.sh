#!/bin/bash

# ww3_run_stats.sh
# PURPOSE: shell script wrapper around python buoy collocation 
#          and generation of statistics and metrics.
#
# USAGE:   ww3_run_stats.sh <ww3.tab.nc> [<var> [<var_min> [<var_max>]]]
#
export PATH=/scratch2/NCEPDEV/marine/Matthew.Masarik/wavpy/miniconda3/bin:${PATH}
export PYTHONPATH=/scratch2/NCEPDEV/marine/Matthew.Masarik/wavpy/miniconda3/pkgs:${PYTHONPATH}

usage_str='ww3_run_stats.sh <ww3.tab.nc> [<var> [<var_min> [<var_max>]]]'
if [[ $# -lt 1 ]] || [[ $# -gt 4 ]]; then echo $usage_str; exit; fi

ww3_tab="$1"
shift
if [ ! -f $ww3_tab ]; then
    echo "file: $ww3_tab, not found"
    exit
fi
echo -e "\nInput file: $ww3_tab"
run_id=${ww3_tab%%.tab.*}
run_id=${run_id#ww3.}

ww3list=ww3list.txt
if [ -f $ww3list ]; then rm -f $ww3list; fi
touch $ww3list
echo -e "$ww3_tab\n$ww3_tab" >> $ww3list


# Buoy-model collocation
buoy_colloc=$(ls -1 ww3buoy_collocation*.nc 2> /dev/null)
if [ -f $buoy_colloc ]; then rm -f $buoy_colloc; fi
echo 'Performing model-buoy collocation...'
python3 modelBuoy_collocation_hindcast.py
num_colloc=$(ls -1 ww3buoy_collocation*.nc | wc -l 2> /dev/null)
if [ $num_colloc -ne 1 ]; then
    echo "invalid num collocation files: ${num_colloc}."
    exit
fi
buoy_colloc=$(ls -1 ww3buoy_collocation*.nc 2> /dev/null)


# Statistics and metrics
if [ -f stats_buoy*.txt  ]; then rm -f stats_buoy*.txt;  fi
if [ -f stats_model*.txt ]; then rm -f stats_model*.txt; fi
if [ -f metrics*.txt     ]; then rm -f metrics*.txt;     fi
python3 runstats.py $buoy_colloc $@


# Data mgmt
mv $buoy_colloc  ${run_id}.${buoy_colloc}
mv $ww3list      ${run_id}.${ww3list}

if [ ! -f stats_buoy*.txt ]; then 
    echo 'buoy stats not found.'
else
    buoy_stats=$(ls -1 stats_buoy*.txt)
    mv $buoy_stats   ${run_id}.${buoy_stats}
fi

if [ ! -f stats_model*.txt ]; then 
    echo 'model stats not found.'
else
    model_stats=$(ls -1 stats_model*.txt)
    mv $model_stats   ${run_id}.${model_stats}
fi

if [ ! -f metrics*.txt ]; then 
    echo 'metrics not found.'
else
    mod_buoy_mets=$(ls -1 metrics*.txt)
    mv $mod_buoy_mets   ${run_id}.${mod_buoy_mets}
fi

echo -e "\nInput WW3 tab:      $ww3_tab\n"
stats_out_dir=stats_out_$run_id
mkdir $stats_out_dir
if [ -d $stats_out_dir ]; then
    mv -f ${run_id}.${buoy_colloc}   $stats_out_dir
    mv -f ${run_id}.${ww3list}       $stats_out_dir
    mv -f ${run_id}.${buoy_stats}    $stats_out_dir
    mv -f ${run_id}.${model_stats}   $stats_out_dir
    mv -f ${run_id}.${mod_buoy_mets} $stats_out_dir
    echo -e "Output Directory:           ${stats_out_dir}/"
fi
echo -e "\tww3list:            ${run_id}.${ww3list}"
echo -e "\tBuoy collocation:   ${run_id}.${buoy_colloc}"
echo -e "\tBuoy stats:         ${run_id}.${buoy_stats}"
echo -e "\tModel stats:        ${run_id}.${model_stats}"
echo -e "\tValidation metrics: ${run_id}.${mod_buoy_mets}\n"

