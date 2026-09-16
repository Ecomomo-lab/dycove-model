#!/bin/bash
#
# Portable Linux example for running DYCOVE + D-Flow FM with MPI.
#
# Required:
#   D3D_HOME       Delft3D FM installation root
#
# Optional:
#   NUM_PROCESSES  MPI rank count (default: 2)
#   MPIEXEC        MPI launcher (default: mpiexec)
#   PYTHON         Python executable (default: python)
#
# Example:
#   export D3D_HOME=/path/to/delft3d/lnx64
#   NUM_PROCESSES=4 ./run_parallel_linux.sh

set -euo pipefail

NUM_PROCESSES=${NUM_PROCESSES:-2}
MPIEXEC=${MPIEXEC:-mpiexec}
PYTHON=${PYTHON:-python}

if [[ -z "${D3D_HOME:-}" ]]; then
    echo "ERROR: D3D_HOME is not set."
    echo "Set D3D_HOME to the Delft3D FM installation root."
    exit 1
fi

if [[ ! "$NUM_PROCESSES" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: NUM_PROCESSES must be a positive integer."
    exit 1
fi

WORK_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DFM_DIR="$WORK_DIR/dflowfm"

export PATH="$D3D_HOME/bin:${PATH:-}"
export LD_LIBRARY_PATH="$D3D_HOME/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

echo "============================================================"
echo "DYCOVE + DFM parallel example"
echo "D3D_HOME       = $D3D_HOME"
echo "NUM_PROCESSES  = $NUM_PROCESSES"
echo "MPIEXEC        = $MPIEXEC"
echo "PYTHON         = $PYTHON"
echo "============================================================"

# -------------------------------------------------------------
# 1. Update DIMR process list from NUM_PROCESSES.
# -------------------------------------------------------------
"$PYTHON" - "$WORK_DIR/dimr_config.xml" "$NUM_PROCESSES" <<'PY'
import re
import sys
from pathlib import Path

config = Path(sys.argv[1])
nprocs = int(sys.argv[2])

processes = " ".join(str(i) for i in range(nprocs))

text = config.read_text()

new_text, count = re.subn(
    r"<process>[^<]*</process>",
    f"<process>{processes}</process>",
    text,
    count=1,
)

if count != 1:
    raise RuntimeError(
        "Expected exactly one <process>...</process> entry "
        "in dimr_config.xml."
    )

config.write_text(new_text)

print(f"DIMR process list: {processes}")
PY

# -------------------------------------------------------------
# 2. Remove old partition products.
# -------------------------------------------------------------
rm -f "$DFM_DIR"/FlowFM_[0-9][0-9][0-9][0-9].mdu
rm -f "$DFM_DIR"/*_[0-9][0-9][0-9][0-9]_net.nc

# -------------------------------------------------------------
# 3. Partition the DFM model.
# -------------------------------------------------------------
echo
echo "Partitioning DFM into $NUM_PROCESSES domains..."

(
    cd "$DFM_DIR"

    "$D3D_HOME/bin/dflowfm" \
        --partition:ndomains="$NUM_PROCESSES" \
        FlowFM.mdu
)

# -------------------------------------------------------------
# 4. Launch DYCOVE + DFM with the same number of MPI ranks.
# -------------------------------------------------------------
echo
echo "Starting MPI simulation..."

cd "$WORK_DIR"

"$MPIEXEC" -np "$NUM_PROCESSES" \
    "$PYTHON" run_tide_channel.py

echo
echo "Parallel DYCOVE + DFM simulation finished."
