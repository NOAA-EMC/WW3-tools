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

set -euo pipefail

# ================= USER'S INPUT =======================
MACHINE="ursa"
WORKDIR="/scratch3/NCEPDEV/marine/Ming.Chen/ursa/ww3tools"

satoutdir="/scratch3/NCEPDEV/climate/Jessica.Meixner/WaveEvaluation/processsatdata/out"    # if satoutdir="", it will use default directory as ${WORKDIR}/processsatdata/out; or user can input directory

SATS=("JASON3" "CRYOSAT2")

NTHREADS=8 # Threads for CDO (usually match cpus-per-task)
# ======================================================
basedir="${WORKDIR}/processsatdata/combineoutmonthly"

if [[ -z "${satoutdir}" ]]; then
    satoutdir="${WORKDIR}/processsatdata/out"
fi

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
    mkdir -p "$basedir"
fi

for SAT in "${SATS[@]}"; do

  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "Processing satellite: ${SAT}"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  in_dir="${satoutdir}/${SAT}"
  if [[ ! -d "$in_dir" ]]; then
    echo "WARNING: input directory not found for $SAT: $in_dir  -> skip"
    continue
  fi

  shopt -s nullglob
  all_files=( "${in_dir}/AltimeterAlongTrack_ww3tools_${SAT}_"*.nc )
  shopt -u nullglob

  if (( ${#all_files[@]} == 0 )); then
    echo "WARNING: no input files found for $SAT in $in_dir"
    continue
  fi

  # Extract unique YYYYMM from filenames
  months=$(
    printf "%s\n" "${all_files[@]}" \
      | sed -E 's|.*_([0-9]{6})[0-9]{4}to.*|\1|' \
      | sort -u
  )

  # Merge each month
  while read -r ym; do
    [[ -z "$ym" ]] && continue

    shopt -s nullglob
    month_files=( "${in_dir}/AltimeterAlongTrack_ww3tools_${SAT}_${ym}"*.nc )
    shopt -u nullglob

    if (( ${#month_files[@]} == 0 )); then
      echo "  $ym: no files -> skip"
      continue
    fi

    # Deterministic order
    IFS=$'\n' month_files_sorted=($(printf "%s\n" "${month_files[@]}" | sort))
    unset IFS

    outfile="${basedir}/Altimeter_${SAT}_${ym}.nc"
    echo "  $ym -> $(basename "$outfile") (nfiles=${#month_files_sorted[@]})"

    # Combine along time dimension (correct for along-track data)
    cdo -O -P "${NTHREADS}" mergetime "${month_files_sorted[@]}" "$outfile"

  done <<< "$months"
done

echo ""
echo "All processing completed."
