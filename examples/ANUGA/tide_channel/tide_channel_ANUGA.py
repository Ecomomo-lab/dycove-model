"""
Simple example script to run an ANUGA hydrodynamic simulation with DYCOVE.

Adapted from ANUGA example "channel1.py"
"""


#------------------------------------------------------------------------------
# Import necessary modules
#------------------------------------------------------------------------------

#import dycove
from dycove import VegetationSpecies, ANUGA
from gen_anuga_domain import RectangSlopeDomainGenerator as RectangDomain

#------------------------------------------------------------------------------
# Create a sloped rectangular ANUGA domain using anuga.rectangular_cross_domain
#------------------------------------------------------------------------------

HydroDomain = RectangDomain("rectang_beach", 
                            rectang_dims=(1000,     # length in x-direction (m)
                                          600),     # length in y-direction (m)
                            mesh_spacing=40,        # exact mesh spacing (regular triangles) (m)
                            min_elev=-0.5,          # elevation at tidal/left boundary (m)
                            dune_elev=0.5,          # elevation of dune crest (m)
                            dune_location=0.65,     # location of dune crest as fraction of domain length
                            channel_bot_elev=0.05,  # bottom elevation of tidal channel
                            lagoon_bot_elev=-0.5,   # bottom elevation of lagoon
                            channel_width_frac=0.1, # width of channel as fraction of domain y
                            tide_props={
                                'amplitude': 0.4,
                                'period': 12.*3600,
                                'MWL': 0},
                            mannings_n=0.025,
                                )

#------------------------------------------------------------------------------
# Run ANUGA with DYCOVE
#------------------------------------------------------------------------------

# define simulation time period
sim_time = 3
time_unit = "eco-morphodynamic years"  # 'hydrodynamic days' or 'eco-morphodynamic years'

# create vegetation species object
veg_1 = VegetationSpecies("veg1.json", "veg1")

# instantiate ANUGA model
HydroModel = ANUGA(HydroDomain.domain, vegetation=veg_1)

# do timestepping
HydroModel.run_simulation(sim_time, sim_time_unit=time_unit)
