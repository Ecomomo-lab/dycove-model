**********
Simulation
**********

This section provides detailed documentation for DYCOVE simulation-related 
classes.


Base Simulation Module
======================

.. currentmodule:: dycove.sim.base

.. autosummary::
    :toctree: ../../_autosummary

    HydroSimulationBase
    HydroEngineBase


Hydrodynamic Model Engines
==========================

ANUGA Hydrodynamic Engine
-------------------------

.. currentmodule:: dycove.sim.engines.ANUGA_hydro

.. autosummary::
    :toctree: ../../_autosummary

    ANUGA
    AnugaEngine


Delft3D-FM Hydrodynamic Engine
------------------------------

.. currentmodule:: dycove.sim.engines.DFM_hydro

.. autosummary::
    :toctree: ../../_autosummary

    DFM
    DFMEngine


Baptist Operator (ANUGA)
------------------------

.. currentmodule:: dycove.sim.engines.ANUGA_baptist

.. autosummary::
    :toctree: ../../_autosummary

    Baptist_operator


Simulation Helper Classes
=========================

Model Coupling
--------------

.. currentmodule:: dycove.sim.coupler

.. autosummary::
    :toctree: ../../_autosummary

    VegetationCoupler

.. currentmodule:: dycove.sim.outputs

Output Management
-----------------

.. autosummary::
    :toctree: ../../_autosummary

    OutputManager

.. currentmodule:: dycove.sim.simulation_data

Simulation Data
---------------

.. autosummary::
    :toctree: ../../_autosummary

    SimulationTimeState
    HydrodynamicStats


Vegetation Modules
==================

Species Classes
---------------

.. currentmodule:: dycove.sim.vegetation

.. autosummary::
    :toctree: ../../_autosummary

    VegetationSpecies
    MultipleVegetationSpecies
    SharedVegMethods

Species Data
------------

.. currentmodule:: dycove.sim.vegetation_data

.. autosummary::
    :toctree: ../../_autosummary

    VegetationAttributes
    VegCohort