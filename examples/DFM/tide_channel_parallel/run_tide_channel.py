"""
Simple example of running Delft3D-FM with DYCOVE in MPI parallel mode.

The Delft3D installation root is supplied through the D3D_HOME
environment variable. DFM domain partitioning and MPI launch are
handled by run_parallel_linux.sh.
"""

import os
from pathlib import Path

from dycove import VegetationSpecies, DFM_hydro


work_dir = Path(__file__).parent.resolve()

config_file = work_dir / "dimr_config.xml"
mdu_file = work_dir / "dflowfm" / "FlowFM.mdu"

try:
    DFM_ROOT = Path(os.environ["D3D_HOME"])
except KeyError as exc:
    raise RuntimeError(
        "D3D_HOME must point to the Delft3D FM installation root."
    ) from exc


sim_time = 4
time_unit = "eco-morphodynamic years"

veg_1 = VegetationSpecies(work_dir / "veg1.json")

HydroModel = DFM_hydro.DFM(
    DFM_ROOT,
    config_file,
    mdu_file,
    vegetation=veg_1,
)

HydroModel.run_simulation(
    sim_time,
    sim_time_unit=time_unit,
    ecofac=50,
    n_ets=14,
    veg_interval=43200,
)
