.. _user-guide:

User Guide
==========

This user guide provides a quick overview of how to use DYCOVE.
Users can get started by reviewing this guide and running the :ref:`Examples <examples>`, but advanced use of DYCOVE requires an understanding of the concepts discussed in the :ref:`Background <background>`.


General Steps
-------------

1. Create a :class:`~dycove.sim.vegetation.VegetationSpecies` instance using a vegetation attribute `.json` file.
2. Create a :class:`~dycove.sim.base.HydroSimulationBase` instance by combining numerical model inputs/pointers with the :class:`~dycove.sim.vegetation.VegetationSpecies` object. 
   This would be called as either :class:`~dycove.sim.engines.ANUGA_hydro.ANUGA` or :class:`~dycove.sim.engines.DFM_hydro.DFM`, depending on which numerical model is used.
3. Run the simulation using :meth:`~dycove.sim.base.HydroSimulationBase.run_simulation` by providing simulation time parameters.
4. Inspect outputs using built-in plotting functions or external tools.


Imports
-------

Import DYCOVE's public API as follows:

.. code-block:: python

   from dycove import VegetationSpecies                      # main vegetation class
   from dycove import MultipleVegetationSpecies as MultiVeg  # for combining multiple species
   from dycove import ANUGA_hydro                            # ANUGA model interface
   from dycove import DFM_hydro                              # Delft3D FM model interface
   from dycove import plotting                               # model plotting module

Of course, users would only import the numerical model interface they plan to use.


Defining the Vegetation Species
-------------------------------

Creating the vegetation species object is the first step in setting up a DYCOVE simulation.
This object contains all of the vegetation attributes needed to run the simulation, read from a user-provided `.json` file (see the :ref:`Background <input-json>` for details on the required inputs).

.. code-block:: python

   veg = VegetationSpecies("path/to/vegetation_attributes.json", "species_name")

Creating a :class:`~dycove.sim.vegetation.VegetationSpecies` object automatically creates a :class:`~dycove.sim.vegetation_data.VegetationAttributes` instance, which contains all of the vegetation attributes read from the input file.
Following the first colonization of the simulation, a :class:`~dycove.sim.vegetation_data.VegCohort` object is also created.


Creating the Model Interface (Delft3D FM)
-----------------------------------------

Instantiate the model interface by calling :class:`~dycove.sim.engines.DFM_hydro.DFM`.
This class inherits :class:`~dycove.sim.base.HydroSimulationBase`, which actually runs the simulation, and instantiates the Delft3D FM "engine" :class:`~dycove.sim.engines.DFM_hydro.DFMEngine`.

.. code-block:: python

   model = DFM_hydro.DFM("path/to/system/files/", 'dimr_config.xml', 'model.mdu', vegetation=veg)

The first argument ``dfm_path`` points to the Delft3D FM system files (`.dll` library).
The directory must contains `dflowfm/` and `dimr/` directories.
On Windows machines, this directory might look like this:

.. code-block:: python

   dfm_path = 'C:\\Program Files (x86)\\Deltares\\Delft3D Flexible Mesh Suite HM (2021.03)\\plugins\\DeltaShell.Dimr\\kernels\\x64'

The second argument ``dimr_config`` is the path to the DIMR configuration file, which will be located in your working directory, after exporting the model as DIMR from the GUI (see :ref:`tide channel example <tide-channel-dfm>`).
The third argument ``mdu_file`` is the name of the Delft3D FM `.mdu` file.


Creating the Model Interface (ANUGA)
------------------------------------

Instantiate the model interface by calling :class:`~dycove.sim.engines.ANUGA_hydro.ANUGA`.
This class inherits :class:`~dycove.sim.base.HydroSimulationBase`, which actually runs the simulation, and instantiates the ANUGA "engine" :class:`~dycove.sim.engines.ANUGA_hydro.AnugaEngine`.

.. code-block:: python

   model = ANUGA_hydro.ANUGA(anuga_domain, vegetation=veg)

The ``anuga_domain`` object is the domain instance created for ANUGA simulations, after all aspects like elevation, friction, forcings, and boundary conditions have been assigned to it.


Multiple Vegetation Species
---------------------------

Users have the option to model multiple species together, using :class:`~dycove.sim.vegetation.MultipleVegetationSpecies`.
A separate `.json` file and :class:`~dycove.sim.vegetation.VegetationSpecies` object are needed for that species, which can be added to model like so (using the ANUGA engine as an example):

.. code-block:: python

   veg_alt = VegetationSpecies("path/to/vegetation_attributes_alt.json", "species_name_alt")
   model = ANUGA_hydro.ANUGA(anuga_domain, vegetation=MultiVeg[veg, veg_alt])


Running the Simulation
----------------------

Users can run a simulation using the following method of :class:`~dycove.sim.base.HydroSimulationBase`, with only the total desired simulation time as an input:

.. code-block:: python

   model.run_simulation(3, sim_time_unit="eco-morphodynamic years")

Note that there are several optional arguments of :meth:`~dycove.sim.base.HydroSimulationBase.run_simulation` that define the time scaling of DYCOVE models, which are described in the API and in greater detail in the :ref:`Background <eco-time-scale>`.


Outputs
-------

Numerical model outputs can be accessed and inspected as they typically would be without vegetation.
For Delft3D FM, the `_his.nc` file for point output and `_map.nc` file for map output are found in the `/dflowfm/output` directory.
For ANUGA, the output `.sww` file is located in the working directory.

DYCOVE vegetation outputs are found in `/dflowfm/veg_output` for Delft3D FM and in `/veg_output` for ANUGA.
DYCOVE output consists of a series of `.npz` (`numpy`) files, one file per vegetation cohort per ecological time step.
One `.npz` file represents the attribute state of a single vegetation cohort at a point in time.
It contains all of the attributes that are stored in a :class:`~dycove.sim.vegetation_data.VegCohort` object (fractions, stem heights, mortality fractions, etc.).

Plotting
--------

Users can plot numerical model results using their existing plotting tools.
However, DYCOVE provides a built-in plotting class :class:`~dycove.utils.plotting.ModelPlotter` that can be used to plot time series of 2-D quantities from both the numerical model and DYCOVE:

.. code-block:: python

   plotter = plotting.ModelPlotter(
      simdir = Path("."),
      quantity = "Stem Height",  # or: quantity = "Velocity", see examples and ModelPlotter API for list of quantities
      plot_times = {
         "plotHR_0": 0*24.,
         "plotHR_f": 21*24.,
         "mapHR_int": 1,
         "plotHR_int": 1,
      },
   )

   plotter.run()

Attributes in the `.npz` file that are 1-D arrays (fractions, mortality fractions) use the same indexing as the underlying mesh, and can be interpolated to 2-D grids using the 1-D ``X`` and ``Y`` arrays in the numerical model output files.
Users can create an interpolation function directly using :func:`~dycove.utils.plotting.create_nn_interpFunc`.
This function is created and called as shown in :meth:`~dycove.utils.plotting.ModelPlotter.get_quantity_grids`.
Loading of vegetation output files is performed in :class:`~dycove.utils.model_loader.BaseMapLoader`, while loading of numerical model outputs occurs in :class:`~dycove.utils.model_loader.ANUGAMapLoader` or :class:`~dycove.utils.model_loader.DFMMapLoader`.
Users can look through these functions and methods to develop their own plotting codes, if that is of interest.