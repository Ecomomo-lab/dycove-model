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
directly import these classes:

>>> import dycove

From a working model script (entry point), ANUGA users looking to model a 
single species should import these high-level classes to get a MWE running:

>>> from dycove import VegetationSpecies, ANUGA

or, for example, for DFM users looking to model multiple species at once:

>>> from dycove import VegetationSpecies, MultipleVegetationSpecies, DFM

Note that users of one numerical model or another may not have installed the
requirements for the other model(s). In the higher-level classes below, the 
try-except statements protect against those kinds of import errors.

"""


# modules called by all entry-point scripts 
from dycove.sim.vegetation import VegetationSpecies, MultipleVegetationSpecies

# Optional backends
__all__ = ["VegetationSpecies", "MultipleVegetationSpecies"]


def _optional_import(name: str):
    """ Attempt to import an optional backend when accessed. """
    import importlib
    import warnings

    try:
        return importlib.import_module(f"dycove.sim.engines.{name}_hydro")
    except ImportError as e:
        warnings.warn(
            f"Optional dependency for '{name.upper()}' coupling not found. "
            f"Install with: pip install dycove[{name}]",
            ImportWarning,
            stacklevel=2,
        )
        raise e


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
ANUGA = _LazyEngine("ANUGA")
DFM = _LazyEngine("DFM")










# for coupling with ANUGA engine
try:
    from dycove.sim.engines.ANUGA_hydro import ANUGA
    __all__.append("ANUGA")
except ImportError:
    import warnings
    warnings.warn(
        "Optional dependency for ANUGA coupling not found. "
        "Install with: pip install dycove[anuga]",
        ImportWarning,
        stacklevel=2,
    )

# for coupling with DFM engine
try:
    from dycove.sim.engines.DFM_hydro import DFM
    __all__.append("DFM")
except ImportError:
    import warnings
    warnings.warn(
        "Optional dependency for Delft3D-FM coupling not found. "
        "Install with: pip install dycove[dfm]",
        ImportWarning,
        stacklevel=2,
    )

try:
    from dycove.utils.plotter import ModelPlotter
    __all__.append("ModelPlotter")
except ImportError:
    import warnings
    warnings.warn(
        "Optional dependency for plotting module not found. "
        "Install with: pip install dycove[plot]",
        ImportWarning,
        stacklevel=2,
    )