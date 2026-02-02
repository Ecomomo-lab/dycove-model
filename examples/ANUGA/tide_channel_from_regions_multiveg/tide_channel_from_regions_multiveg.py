"""
Simple example script to run an ANUGA hydrodynamic simulation with DYCOVE.

Adapted from ANUGA example `channel1.py`.
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

from dycove import VegetationSpecies, ANUGA_hydro
from dycove import MultipleVegetationSpecies as MultiVeg
from gen_anuga_domain_from_regions import RectangSlopeDomainGenerator as RectangDomain

#------------------------------------------------------------------------------
# Create a sloped rectangular ANUGA domain using anuga.rectangular_cross_domain
#------------------------------------------------------------------------------

HydroDomain = RectangDomain("rectang_beach",
                            ["exterior.csv", "interior.csv"],
                            [800, 200],
                            plotting=True
                            )

#------------------------------------------------------------------------------
# Run ANUGA with DYCOVE
#------------------------------------------------------------------------------

# define simulation time period
sim_time = 4
time_unit = "eco-morphodynamic years"  # 'hydrodynamic days' or 'eco-morphodynamic years'

# create vegetation species objects
veg_1 = VegetationSpecies("NelumboLutea.json")
veg_2 = VegetationSpecies("ColocasciaEsculenta.json")

# instantiate ANUGA model
HydroModel = ANUGA_hydro.ANUGA(HydroDomain.domain, vegetation=MultiVeg([veg_1, veg_2]))

# do timestepping
HydroModel.run_simulation(sim_time, sim_time_unit=time_unit)
