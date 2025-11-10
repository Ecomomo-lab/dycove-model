"""
The DYCOVE (DYnamic COastal VEgetation) model provides an interface for the 
tight coupling of dynamic vegetation processes to existing numerical models 
(ANUGA and Delft3D FM for now) of coastal hydrodynamics and morphodynamics 
(in the latter case). The user provides an existing numerical model, such as
an MDU file and supporting files for DFM or a Python domain object for ANUGA,
and one or more vegetation input files (one per species) that contain species
characteristics related to colonization, growth, and mortality. The user 
instantiates the model based on these inputs and a simulation time frame
(including a coupling interval).

This file contains the classes that form the basis of the DYCOVE model, and 
can be imported directly without specifying complete file paths as long 
import statements. The simple import statement below allows the user to 
directly import these classes.

>>> import dycove

Note that users of one numerical model or another may not have installed the
requirements for the other model(s). In the higher-level classes below, the 
try-except statements protect against those kinds of import errors.

"""

# ----------------------------------------------------------------
# Higher-level classes (called by the user in entry-point script)
# ----------------------------------------------------------------

# called by the user in entry-point script
from dycove.sim.vegetation import VegetationSpecies
try:
    from dycove.sim.engines.ANUGA_hydro import ANUGA
    from dycove.utils.baptist_operator import Baptist_operator
except:
    pass
try:
    from dycove.sim.engines.DFM_hydro import DFM
except:
    pass

# ----------------------------------------------------------------
# Lower-level classes (called by other scripts/classes)
# ----------------------------------------------------------------

# Base hydrodynamic simulation classes
from dycove.sim.base import HydroSimulationBase, HydroEngineBase
from dycove.utils.log import Reporter
from dycove.sim.coupler import VegetationCoupler
from dycove.sim.simulation_data import SimulationTimeState, HydrodynamicStats
from dycove.sim.outputs import OutputManager
from dycove.utils.simulation_reporting import print_model_time_info, print_runtime_updates

from dycove.sim.vegetation_data import VegetationAttributes, VegCohort
from dycove.utils.array_math import cell_averaging, sum_product, sum_elementwise
