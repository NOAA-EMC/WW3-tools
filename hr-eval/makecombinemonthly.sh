#!/bin/bash
#SBATCH --nodes=1
#SBATCH -q batch
#SBATCH -t 08:00:00
#SBATCH -A marine-cpu
#SBATCH -J makecombinemonthly
#SBATCH -o makecombinemonthly.o%j
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# ================= USER'S INPUT =======================
MACHINE="ursa"

satoutdir=/scratch3/NCEPDEV/climate/Jessica.Meixner/WaveEvaluation/processsatdata/out    # processed satellite data need to be combined
basedir=/scratch4/NCEPDEV/marine/Ming.Chen/wave_eval/processsatdata/combineoutmonthly # output combined data

SATS=('SENTINEL6A')

months=("08"       "09"  "10"    "03"   "04"    "05"   "06"   "07"   "08"   "09"   "10"   "11"   "12"   "01"   "02"   "03"   "04"   "05")
nextmonths=("08"   "09"  "10"    "04"   "05"    "06"   "07"   "08"   "09"   "10"   "11"   "12"   "01"   "02"   "03"   "04"   "05"   "06")
years=("2022"      "2022" "2022" "2024" "2024"  "2024" "2024" "2024" "2024" "2024" "2024" "2024" "2024" "2025" "2025" "2025" "2025" "2025")
nextyears=("2022"  "2022" "2022" "2024" "2024"  "2024" "2024" "2024" "2024" "2024" "2024" "2024" "2025" "2025" "2025" "2025" "2025" "2025")
# ======================================================

if [[ "$MACHINE" == "ursa" ]]; then
    module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core
    module use /contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/intel-oneapi-mpi/2021.13-haww6b3/gcc/12.4.0
    module load stack-oneapi/2024.2.1
    module load stack-intel-oneapi-mpi/2021.13
    module load cdo/2.4.4

elif [[ "$MACHINE" == "hercules" ]]; then
    module use /apps/contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.1.0/install/modulefiles/Core
    module use /apps/contrib/spack-stack/spack-stack-1.9.2/envs/ue-oneapi-2024.1.0/install/modulefiles/intel-oneapi-mpi/2021.13-sqiixt7/gcc/13.3.0
    module load stack-oneapi/2024.2.1
    module load stack-intel-oneapi-mpi/2021.13
    module load cdo/2.4.4

else
    echo "ERROR: No module configuration defined for machine $MACHINE"
    exit 1
fi

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# check directories
if [[ ! -d "$satoutdir" ]]; then
    echo "ERROR: Satellite input directory does not exist: $satoutdir"
    exit 1
fi

if [[ ! -d "$basedir" ]]; then
    echo "Output directory does not exist → creating: $basedir"
    mkdir -p "$basedir" || {
        echo "ERROR: Failed to create basedir: $basedir"
        exit 1
    }
fi

for SAT in "${SATS[@]}"; do

  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "Processing satellite: ${SAT}"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  cd ${basedir}
  for i in ${!months[@]}; do
    month=${months[$i]}
    year=${years[$i]}
    nextmonth=${nextmonths[$i]}
    nextyear=${nextyears[$i]}
    date=${year}${month}
    mkdir ${date}${SAT}
    cd ${date}${SAT}
    cp ${satoutdir}/${SAT}/AltimeterAlongTrack_ww3tools_${SAT}_${date}*.nc .
    cp ${satoutdir}/${SAT}/AltimeterAlongTrack_ww3tools_${SAT}_${nextyear}${nextmonth}01*.nc .
    cdo mergetime AltimeterAlongTrack_ww3tools_${SAT}_*.nc Altimeter_${SAT}_${date}.nc
    mv Altimeter_${SAT}_${date}.nc ${basedir}/
    cd ${basedir}
  done
done
echo "All processing completed."
