.. _tide-channel-anuga-multiveg:

Example - Simple beach and tide channel with two vegetation species (ANUGA)
===========================================================================

This example corresponds to the Python file `tide_channel_from_regions_multiveg.py` located in the `examples/ANUGA/tide_channel/` directory, along with supporting input files.

This example consists of a simple, symmetrical beach slope with a channel that bisects the dune and connects to a lagoon.
An unstructured mesh is developed using ANUGA's ``create_domain_from_regions`` method that creates a triangular mesh with variable resolution based on input polygons.
A tidal boundary condition is imposed on the left (ocean) side, that propagates through the channel and fills the lagoon behind the dune.
Two vegetation species, loosely based on `Nelumbo Lutea` and `Colocascia Esculenta`, are added to the model via a set of attributes defined in the files `NelumboLutea.json` and `ColocasciaEsculenta.json`.
After the simulation is finished, we can inspect a variety of 2-D model results using :class:`~dycove.utils.plotting.ModelPlotter`.

Running the model and plotting the results are handled mostly as described in the :doc:`tide_channel example <tide_channel>`.
Users running this example can access the associated ANUGA domain creator class (different file name from previous examples), which uses the ``create_domain_from_regions`` method to construct a mesh:

.. code-block:: python

   from gen_anuga_domain_from_regions import RectangSlopeDomainGenerator as RectangDomain
   HydroDomain = RectangDomain("rectang_beach",
                        ["exterior.csv", "interior.csv"],
                        [800, 200],
                        plotting=True
                        )

ANUGA's ``create_domain_from_regions`` method requires input polygons (typically as `.csv` files) that define the exterior domain boundary as well as any number of non-overlapping, interior regions.
Each polygon is assigned a maximum triangle area.
In this example, we supply a file `exterior.csv` that contains relative vertex coordinates mimicking the dimensions of the :doc:`tide_channel example <tide_channel>`.
The file `interior.csv` represents the boundary of the tidal channel, and we supply a smaller triangle size for this region to increase resolution.
Interior boundaries can also be used to form breaklines in the mesh, which is one of its purposes in this domain.
By setting ``plotting=True``, ``RectangSlopeDomainGenerator`` will save images in the working directory showing the resulting mesh and mesh elevation.

For modeling of multiple species, the only change needed is in creating multiple :class:`~dycove.sim.vegetation_data.VegetationSpecies` objects:

.. code-block:: python

   from dycove import MultipleVegetationSpecies as MultiVeg
   veg_1 = VegetationSpecies("NelumboLutea.json")
   veg_2 = VegetationSpecies("ColocasciaEsculenta.json")
   HydroModel = ANUGA_hydro.ANUGA(HydroDomain.domain, vegetation=MultiVeg([veg_1, veg_2]))
   HydroModel.run_simulation(4, sim_time_unit="eco-morphodynamic years")

For notes on plotting results with multiple vegetation species, see the :ref:`User Guide <plot-multi-species>`.