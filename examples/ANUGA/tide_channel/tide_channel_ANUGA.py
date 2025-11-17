"""
Simple example script to run an ANUGA hydrodynamic simulation with DYCOVE.

Adapted from ANUGA example `channel1.py`.
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from dycove import VegetationSpecies, ANUGA_hydro
from gen_anuga_domain import RectangSlopeDomainGenerator as RectangDomain

#------------------------------------------------------------------------------
# Create a sloped rectangular ANUGA domain using anuga.rectangular_cross_domain
#------------------------------------------------------------------------------

HydroDomain = RectangDomain("rectang_beach")

#------------------------------------------------------------------------------
# Run ANUGA with DYCOVE
#------------------------------------------------------------------------------

# define simulation time period
sim_time = 1
time_unit = "eco-morphodynamic years"  # 'hydrodynamic days' or 'eco-morphodynamic years'

# create vegetation species object
veg_1 = VegetationSpecies("veg1.json", "veg1")

# instantiate ANUGA model
HydroModel = ANUGA_hydro.ANUGA(HydroDomain.domain, vegetation=veg_1)

# do timestepping
HydroModel.run_simulation(sim_time, sim_time_unit=time_unit)
