"""
Simple example script to run a Delft3D-FM hydrodynamic simulation with DYCOVE (with morphology).

Adapted from ANUGA-DYCOVE example "simple_tide_ANUGA.py"
"""

#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from pathlib import Path

from dycove import VegetationSpecies, DFM_hydro
from dycove import MultipleVegetationSpecies as MultiVeg

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
veg_1 = VegetationSpecies("NelumboLutea.json", 
                          mor=1,
                          rand_seed_frac=0.8,
                          )
veg_2 = VegetationSpecies("ColocasciaEsculenta.json", 
                          mor=1,
                          rand_seed_frac=0.8,
                          )

# instantiate DFM model
HydroModel = DFM_hydro.DFM(DFM_DLLs, config_file, mdu_file, vegetation=MultiVeg([veg_1, veg_2]))

# do timestepping
HydroModel.run_simulation(sim_time, time_unit)
