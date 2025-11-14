Example - Simple beach and tide channel (ANUGA)
===============================================

This example corresponds to the Python file "tide_channel_ANUGA.py" located in the "examples/ANUGA/tide_channel/" directory, along with supporting input files.

This example consists of a simple, symmetrical beach slope with a channel that bisects the dune and connects to a lagoon.
A simple mesh is developed using one of ANUGA's built-in functions that creates a rectangular mesh of structured triangles.
A tidal boundary condition is imposed on the left (ocean) side, that propagates through the channel and fills the lagoon behind the dune.
A vegetation species, loosely based on Nelumbo Lutea (American lotus), is added to the model via a set of parameters defined in the file `veg1.json`.

The following import statement forms the basis of all DYCOVE-ANUGA coupling simulations:

.. doctest::

   >>> from dycove import VegetationSpecies, ANUGA_hydro

The class ``RectangSlopeDomainGenerator`` in the file `gen_anuga_domain.py` is a utility class that was developed to generate the ANUGA domain for this example.
The new variable ``HydroDomain`` holds the ANUGA domain object required to run ANUGA models:

.. doctest::

   >>> from gen_anuga_domain import RectangSlopeDomainGenerator as RectangDomain
   >>> HydroDomain = RectangDomain("rectang_beach")

While this class is not part of DYCOVE, it demonstrates how to set up a simple domain in ANUGA.
The majority of the inputs to this class are related to describing the exact geometry of the problem.
In general, the most critical aspects of the ANUGA domain setup are the following lines (written generally):

.. doctest::

   domain = anuga.rectangular_cross_domain(...)  # there are other ways to create domains...
   domain.set_quantity("elevation", topography)  # can also set elevation from a file (e.g., .ASC)
   domain.set_quantity("friction", friction)     # same as above
   Br = anuga.Reflective_boundary(domain)        # define the reflective BC
   Bt = anuga.Time_boundary(domain, function=lambda t: [tidal_func(t), 0.0, 0.0])  # define tidal BC
   domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})  # assign BCs to each domain side

More important than the specifics of the domain geometry is the implementation of vegetation with DYCOVE. 
Instantiate a ``VegetationSpecies`` object using the vegetation parameter file, then pass that object to the ``ANUGA`` hydrodynamic engine:

.. doctest::

   >>> veg_1 = VegetationSpecies("veg1.json", "veg1")
   >>> HydroModel = ANUGA_hydro.ANUGA(HydroDomain.domain, vegetation=veg_1)

Finally, run the simulation for a specified number of eco-morphodynamic years (here, 3 years):

   >>> HydroModel.run_simulation(3, sim_time_unit="eco-morphodynamic years")

Note that above command to run the simulation uses default values for ecological time scaling, namely, ``n_ets``, ``veg_interval``, and ``ecofac`` (see :meth:`~dycove.sim.base.HydroSimulationBase.run_simulation' for details).
For further explanation of the ecological time scaling logic used in DYCOVE, refer to the `background documentation <https://Ecomomo-lab.github.io/dycove-model/background/ecological_time_scaling.html>`_.