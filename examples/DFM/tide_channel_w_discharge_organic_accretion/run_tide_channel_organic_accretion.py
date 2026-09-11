"""
Simple example script to run Delft3D FM with DYCOVE, morphology,
discharge, and vegetation-driven organic accretion.

Adapted from ANUGA-DYCOVE example.
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from pathlib import Path

from dycove import DFM_hydro, OrganicAccretion, VegetationSpecies

#------------------------------------------------------------------------------
# Define model file locations
#------------------------------------------------------------------------------

work_dir = Path(__file__).parent
config_file = work_dir / 'dimr_config.xml'
mdu_file    = work_dir / 'dflowfm/FlowFM.mdu'

DFM_DLLs = Path('C:/Program Files (x86)/Deltares/Delft3D Flexible Mesh Suite HMWQ (2021.03)/'
                'plugins/DeltaShell.Dimr/kernels/x64')  # path to the local Delft3D software folder

#------------------------------------------------------------------------------
# Run Delft3D-FM with DYCOVE
#------------------------------------------------------------------------------

# define simulation time period
sim_time = 4
time_unit = "eco-morphodynamic years"  # 'hydrodynamic days' or 'eco-morphodynamic years'

# create vegetation species object
veg_1 = VegetationSpecies(work_dir / "veg1.json", mor=1)

# configure vegetation-driven organic accretion
organic = OrganicAccretion(
    bed_type="mixing_layer",
    target_sediment=1,
)

# instantiate DFM model
HydroModel = DFM_hydro.DFM(DFM_DLLs, config_file, mdu_file, vegetation=veg_1, organic=organic)

# do timestepping
HydroModel.run_simulation(sim_time, sim_time_unit=time_unit, n_ets=14, veg_interval=43200)
