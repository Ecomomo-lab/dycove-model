"""
Simple example script to run an ANUGA hydrodynamic simulation with DYCOVE.

Adapted from ANUGA example `channel1.py`.

For running in parallel, locally.
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from dycove import VegetationSpecies, ANUGA_hydro
from gen_anuga_domain_pll import RectangSlopeDomainGenerator as RectangDomain

#------------------------------------------------------------------------------
# Create a sloped rectangular ANUGA domain using anuga.rectangular_cross_domain
#------------------------------------------------------------------------------

HydroDomain = RectangDomain("rectang_beach", mesh_spacing=20)

#------------------------------------------------------------------------------
# Run ANUGA with DYCOVE
#------------------------------------------------------------------------------

# define simulation time period
sim_time = 4
time_unit = "eco-morphodynamic years"  # 'hydrodynamic days' or 'eco-morphodynamic years'

# create vegetation species object
veg_1 = VegetationSpecies("veg1.json")

# instantiate ANUGA model
HydroModel = ANUGA_hydro.ANUGA(HydroDomain.domain, vegetation=veg_1)

# do timestepping
HydroModel.run_simulation(sim_time, time_unit)
