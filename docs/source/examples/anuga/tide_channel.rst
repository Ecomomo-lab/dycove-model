.. _tide-channel-anuga:

Example - Simple beach and tide channel (ANUGA)
===============================================

This example corresponds to the Python file `tide_channel_ANUGA.py` located in the `examples/ANUGA/tide_channel/` directory, along with supporting input files.

This example consists of a simple, symmetrical beach slope with a channel that bisects the dune and connects to a lagoon.
A simple mesh is developed using one of ANUGA's built-in functions that creates a rectangular mesh of structured triangles.
A tidal boundary condition is imposed on the left (ocean) side, that propagates through the channel and fills the lagoon behind the dune.
A vegetation species, loosely based on `Salicorniaspp.` (pioneer species, first life stage) and `Spartina anglica` (second life stage), is added to the model via a set of attributes defined in the file `veg1.json`.
After the simulation is finished, we can inspect a variety of 2-D model results using :class:`~dycove.utils.plotting.ModelPlotter`.


Run the Model
-------------

The following import statement forms the basis of all DYCOVE-ANUGA coupling simulations:

.. code-block:: python

   >>> from dycove import VegetationSpecies, ANUGA_hydro

The class ``RectangSlopeDomainGenerator`` in the file `gen_anuga_domain.py` is a utility class that was developed to generate the ANUGA domain for this example.
The new variable ``HydroDomain`` holds the ANUGA domain object required to run ANUGA models:

.. code-block:: python

   >>> from gen_anuga_domain import RectangSlopeDomainGenerator as RectangDomain
   >>> HydroDomain = RectangDomain("rectang_beach")

While this class is not part of DYCOVE, it demonstrates how to set up a simple domain in ANUGA.
The majority of the inputs to this class are related to describing the exact geometry of the problem.
In general, the most critical aspects of the ANUGA domain setup are the following lines (written generally):

.. code-block:: python

   domain = anuga.rectangular_cross_domain(...)  # there are other ways to create domains...
   domain.set_quantity("elevation", topography)  # can also set elevation from a file (e.g., .ASC)
   domain.set_quantity("friction", friction)     # same as above
   Br = anuga.Reflective_boundary(domain)        # define the reflective BC
   Bt = anuga.Time_boundary(domain, function=lambda t: [tidal_func(t), 0.0, 0.0])  # define tidal BC
   domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})  # assign BCs to each domain side

More important than the specifics of the domain geometry is the implementation of vegetation with DYCOVE. 
Instantiate a ``VegetationSpecies`` object using the vegetation attribute file, then pass that object to the :class:`~dycove.sim.engines.ANUGA_hydro.ANUGA` hydrodynamic engine:

.. code-block:: python

   >>> veg_1 = VegetationSpecies("veg1.json", "veg1")
   >>> HydroModel = ANUGA_hydro.ANUGA(HydroDomain.domain, vegetation=veg_1)

The `veg1.json` attribute file found in the working example directory contains a number of parameters related to when/where this species will colonize, how it will grow, and under what conditions it will die off.
The variables in this file map directly to :class:`~dycove.sim.vegetation_data.VegetationAttributes`, a class that contains the documentation for all required input variables.

Finally, run the simulation for a specified number of eco-morphodynamic years (here, 3 years):

.. code-block:: python

   >>> HydroModel.run_simulation(3, sim_time_unit="eco-morphodynamic years")

Note that above command to run the simulation uses default values for ecological time scaling, namely, the number of ecological time steps per year ``n_ets``, hydrodynamic time between ecological time steps ``veg_interval``, and ecological scaling factor ``ecofac`` (see :meth:`~dycove.sim.base.HydroSimulationBase.run_simulation` for details and default values).
For further explanation of the ecological time scaling logic used in DYCOVE, refer to the `background documentation <https://Ecomomo-lab.github.io/dycove-model/background/ecological_time_scaling.html>`_.
Using the default values of ``n_ets`` and ``veg_interval``, ``ecofac`` is equal to 52.
So, three years of ecological time is equivalent to 21 days of hydrodynamic time.
Therefore, the above line can also be written as below, with the same resulting simulation:

.. code-block:: python

   HydroModel.run_simulation(21, sim_time_unit="hydrodynamic days")


Plot the Results
----------------

ANUGA output is saved in netCDF files with the `.sww` extension. 
DYCOVE output is saved in `.npz` files, with one file per vegetation cohort per ecological time step.
For a list of outputs stored with each cohort, see :class:`~dycove.sim.vegetation_data.VegCohort`.

Plot 2-D maps of various model quantities over time using :class:`~dycove.utils.plotting.ModelPlotter`.
For example, we can plot vegetation stem heights over the entire 21-day (hydrodynamic time) simulation period.

.. code-block:: python

   >>> from dycove import plotting
   >>> plotter = plotting.ModelPlotter(
         simdir = Path("."),
         quantity = "Stem Height",
         plot_times = {
           "plotHR_0": 0*24.,
           "plotHR_f": 21*24.,
           "mapHR_int": 1,
           "plotHR_int": 1,
         },
         cmap_lims = {
           "Bathymetry": (-0.5, 0.5),
         },
       )
   >>> plotter.run()

This plotting code is located in the same example directory at `examples/ANUGA/tide_channel/plot_tide_channel_ANUGA.py`, where there are additional comments on how to use the code, as well as more ideas for quantities to plot.
