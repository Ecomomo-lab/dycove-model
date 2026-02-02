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

This file provides import shortcuts for the public-facing classes that form 
the basis of the DYCOVE model, and can be imported directly without specifying 
complete file paths as long import statements. For example, from a working 
model script (entry point), ANUGA users looking to model a single species 
should import these high-level classes to get a DYCOVE model running:

>>> from dycove import VegetationSpecies, ANUGA_hydro

...create an ANUGA domain:

>>> domain = anuga.create_domain_from_regions(...)

...and then instantiate and run the coupled model:

>>> veg_1 = VegetationSpecies("veg1.json")
>>> model = ANUGA_hydro.ANUGA(domain, vegetation=veg_1)
>>> model.run_simulation(3)

Or, for Delft3D FM users looking to model multiple species at once:

>>> from dycove import VegetationSpecies, DFM_hydro
>>> from dycove import MultipleVegetationSpecies as MultiVeg
>>> veg_1 = VegetationSpecies("veg1.json")
>>> veg_2 = VegetationSpecies("veg2.json")

...define paths to DFM executables, DIMR config file, and model MDU files, then run:

>>> model = DFM_hydro.DFM(DFM_DLL_path, 'dimr_config.xml', 'FlowFM.mdu', vegetation=MultiVeg([veg_1, veg_2])
>>> model.run_simulation(3)

Note that users of one numerical model or another may not have installed the
requirements for the other model(s). For the higher-level classes below, the 
_optional_import/_LazyEngine protect against those kinds of import errors.

"""


# modules called by all entry-point scripts 
from dycove.sim.vegetation import VegetationSpecies, MultipleVegetationSpecies

# Optional backends
__all__ = [
    "VegetationSpecies", 
    "MultipleVegetationSpecies",
    "ANUGA_hydro",
    "DFM_hydro",
    "plotting",
    ]

import importlib
# import sys
# import types
import warnings

# ---------------------------------------------------------------------
# Helper: optional imports
# ---------------------------------------------------------------------
def _optional_import(name: str):
    """ Attempt to import an optional backend when accessed. """

    import_paths = {
        "DFM": "dycove.sim.engines.DFM_hydro",
        "ANUGA": "dycove.sim.engines.ANUGA_hydro",
        "plot": "dycove.utils.plotting",
    }

    try:
        return importlib.import_module(import_paths[name])
    except ImportError as e:
        warnings.warn(
            f"Optional dependency for '{name}' not found. "
            f"For dependency checks, install with: pip install dycove[{name.lower()}]",
            ImportWarning,
            stacklevel=2,
        )
        raise e

# ---------------------------------------------------------------------
# Lazy wrapper class
# ---------------------------------------------------------------------
class _LazyEngine:
    """ Lazy loading of an optional coupling engine when first accessed. """

    def __init__(self, name):
        self._name = name
        self._module = None

    def __getattr__(self, attr):
        if self._module is None:
            self._module = _optional_import(self._name)
        return getattr(self._module, attr)


# Expose optional engines lazily
ANUGA_hydro = _LazyEngine("ANUGA")
DFM_hydro = _LazyEngine("DFM")
plotting = _LazyEngine("plot")
