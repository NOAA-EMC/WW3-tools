#!/bin/bash

set -euo pipefail
set -x


MACHINE="ursa"
WORKDIR="/scratch4/NCEPDEV/marine/Ming.Chen/ww3tools_dev"

SCRIPTDIR="${WORKDIR}/WW3-tools/ush"
PY_SCRIPT="create_ww3list.py"

DATA_DIR="/scratch3/NCEPDEV/climate/Jessica.Meixner/Data/gfsv16"

START_DATE="20250501"
END_DATE="20250510"

CYCLES=("00" "06" "12" "18")

FILENAME_PATTERN="gfswave.t{cycle}z.bull_tar"

OUTDIR="${WORKDIR}/src"
OUTFILE="ww3list.txt"


if [[ "${MACHINE}" == "ursa" ]]; then
    module use /scratch3/NCEPDEV/climate/Jessica.Meixner/general/modulefiles
    module load ww3tools
elif [[ "${MACHINE}" == "orion" || "${MACHINE}" == "hercules" ]]; then
    module use /work2/noaa/marine/jmeixner/general/modulefiles
    module load ww3tools
else
    echo "ERROR: Unsupported MACHINE='${MACHINE}'"
    exit 1
fi

echo "====================================================="
echo "Runtime configuration"
echo "====================================================="
echo "MACHINE            = ${MACHINE}"
echo "WORKDIR            = ${WORKDIR}"
echo "SCRIPTDIR          = ${SCRIPTDIR}"
echo "DATA_DIR           = ${DATA_DIR}"
echo "START_DATE         = ${START_DATE}"
echo "END_DATE           = ${END_DATE}"
echo "CYCLES             = ${CYCLES[*]}"
echo "FILENAME_PATTERN   = ${FILENAME_PATTERN}"
echo "OUTDIR             = ${OUTDIR}"
echo "OUTFILE            = ${OUTFILE}"
echo "====================================================="

CMD=(
    python -u "${SCRIPTDIR}/${PY_SCRIPT}"
    -d "${DATA_DIR}"
    -s "${START_DATE}"
    -e "${END_DATE}"
    -c "${CYCLES[@]}"
    -p "${FILENAME_PATTERN}"
    -o "${OUTDIR}"
    -f "${OUTFILE}"
)

echo "COMMAND:"
printf '%q ' "${CMD[@]}"
echo
echo "====================================================="

"${CMD[@]}"

echo "Done."

