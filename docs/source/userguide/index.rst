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

   veg = VegetationSpecies("path/to/vegetation_attributes.json")

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

   veg_alt = VegetationSpecies("path/to/vegetation_attributes_alt.json")
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

However, DYCOVE provides a built-in plotting class :class:`~dycove.utils.plotting.ModelPlotter` that can be used to plot time series of 2-D quantities from both the numerical model and DYCOVE:

.. code-block:: python

   plotter = plotting.ModelPlotter(
      simdir = Path("."),
      # or: quantity = "Velocity", etc. See examples and ModelPlotter API for list of quantities
      quantity = "Stem Height",
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


.. _plot-multi-species:

Plotting Results with Multiple Species
--------------------------------------

Plotting results of simulations with multiple species is no different from plotting results of a single-species simulation.
Plotting stem height using :class:`~dycove.utils.plotting.ModelPlotter`, for example, will average all stem heights in each grid cell among the various species and life stages present within that cell.
When plotting vegetation fractions, however, the output images will alternate between fractions of each species.
Example: A simulation of one ecological year and two species consists of one colonization event for each species in the second ecological time step (ETS) of the year.
Plots of stem heights will consist of one plot per ETS, starting at ETS 2.
Plots of fractions will consist of two plots per ETS, one for each species, starting at ETS 2.


Accessing and Plotting Outputs Directly
---------------------------------------

The following explanations and examples are for those users that prefer to develop their own plotting codes.
ANUGA model quantities can be read from the output `.sww` file using the ``xarray`` library (also ``netCDF4``):

.. code-block:: python

   map_vars = xr.open_dataset("path/to/anuga_model.sww")
   x_nodes = map_vars["x"]
   y_nodes = map_vars["y"]
   z_nodes = map_vars["elevation"]
   stage = map_vars["stage"]
   depth = stage - elevation

These are 1-D arrays, that can be interpolated using :func:`~dycove.utils.plotting.create_nn_interpFunc` or other tools.
Note that `stage` and other time-varying quantities have a time axis as the first axis (same with DFM).
Centroid data is also included in the `.sww` file, but `x` and `y` values at centroids would need to be computed as they are in ``~dycove.utils.model_loader.ANUGAMapLoader._load_outputs``.
Similarly, DFM model (centroid) quantities can be read from the output `_map.nc` file:

.. code-block:: python

   map_vars = xarray.open_dataset("path/to/dflowfm/output/FlowFM_map.nc")
   x_nodes = map_vars["x"]
   y_nodes = map_vars["y"]
   z_nodes = map_vars["elevation"]
   stage = map_vars["stage"]
   depth = stage - elevation

Vegetation quantities can be accessed via the files in the `veg_output` directory.
Quantities that are saved in these output files are listed as attributes of the :class:`~dycove.sim.vegetation_data.VegCohort` class.
For example:

.. code-block:: python

   # Loop through all cohort files saved for a given year and ETS, load quantities to running lists
   veg_fractions, veg_quantity = [], []
   for file in self.ecodir.glob(f'cohort*_year1_ets7.npz'):
      c = dict(numpy.load("path/to/veg_output/cohort1_year1_ets7.npz", allow_pickle=True))
      fractions.append(c["fraction"])
      heights.append(c["height"])  # or "density", "applied_mort_flood", etc.

Note that for a given cohort, each ``c["fraction"]`` is an array with one value per grid cell, while other quantities like ``c["height"]`` have only a single value per cohort, which grows over time.
Mortality outputs are arrays, similar to Fractions.
While vegetation fractions and cohorts are tracked individually in DYCOVE, we typically plot them as a weighted average in each grid cell, with weights determined by the fractions present in the cell.
First, we convert the data to a grid using :func:`~dycove.utils.plotting.create_nn_interpFunc`, then we use the helper functions in :py:mod:`~dycove.utils.array_math` to compute weighted averages:

.. code-block:: python

   from dycove.utils.array_math import cell_averaging, sum_product
   grid_data = cell_averaging(fractions, heights)

With the data as a 2-D array, it can be plotted as usual using ``imshow``.
Users can (and should) loop over "years" and "ETS" to view species evolution through time.