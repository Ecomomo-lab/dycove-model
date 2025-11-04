"""
Simple example script to run an ANUGA hydrodynamic simulation with DYCOVE.

Adapted from ANUGA example "channel1.py"
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from dycove.sim.vegetation import VegetationSpecies
from dycove.sim.engines.ANUGA_hydro import ANUGA
from gen_anuga_domain import RectangSlopeDomainGenerator as RectangDomain


#------------------------------------------------------------------------------
# Create a sloped rectangular ANUGA domain using anuga.rectangular_cross_domain
#------------------------------------------------------------------------------


HydroDomain = RectangDomain("rectang_beach", 
                            rectang_dims=(400, 200), 
                            mesh_spacing=10,
                            min_elev=-0.5,
                            slope=0.0025,
                            mannings_n=0.025,
                            tide_props={
                                'amplitude': 0.5,
                                'period': 12.*3600,
                                'MWL': 0.25},
                                )


#------------------------------------------------------------------------------
# Run ANUGA with DYCOVE
#------------------------------------------------------------------------------

# define simulation time period
sim_time = 3
time_unit = "eco-morphodynamic years"  # 'hydrodynamic days' or 'eco-morphodynamic years'

# create vegetation species object
veg_1 = VegetationSpecies("veg1.txt", "veg1")

# instantiate ANUGA model
HydroModel = ANUGA(HydroDomain, vegetation=veg_1)

# start timestepping
HydroModel.run_simulation(sim_time, sim_time_unit=time_unit)
