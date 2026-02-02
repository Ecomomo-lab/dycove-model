"""
Simple example script to run a Delft3D-FM hydrodynamic simulation with DYCOVE (no morphology).

Adapted from ANUGA-DYCOVE example "simple_tide_ANUGA.py"
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from pathlib import Path

from dycove import VegetationSpecies, DFM_hydro


#------------------------------------------------------------------------------
# Define model file locations
#------------------------------------------------------------------------------

work_dir = Path(__file__).parent
config_file = work_dir / 'dimr_config.xml'
mdu_file    = work_dir / 'dflowfm/FlowFM.mdu'

DFM_DLLs = Path('C:/Program Files (x86)/Deltares/Delft3D Flexible Mesh Suite HM (2021.03)/'
                'plugins/DeltaShell.Dimr/kernels/x64')  # path to the local Delft3D software folder

#------------------------------------------------------------------------------
# Run Delft3D-FM with DYCOVE
#------------------------------------------------------------------------------

# define simulation time period
sim_time = 4
time_unit = "eco-morphodynamic years"  # 'hydrodynamic days' or 'eco-morphodynamic years'

# create vegetation species object
veg_1 = VegetationSpecies("veg1.json")

# instantiate DFM model
HydroModel = DFM_hydro.DFM(DFM_DLLs, config_file, mdu_file, vegetation=veg_1)

# do timestepping
HydroModel.run_simulation(sim_time, sim_time_unit=time_unit)
