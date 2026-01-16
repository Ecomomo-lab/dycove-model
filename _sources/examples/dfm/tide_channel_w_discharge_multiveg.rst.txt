.. _tide-channel-dfm- multiveg:

Example - Simple beach and tide channel with discharge boundary and morphology, with two vegetation species (Delft3D FM)
========================================================================================================================


This example corresponds to the Python file `tide_channel_w_discharge_multiveg.py` located in the `examples/DFM/tide_channel_w_discharge_multiveg/` directory, along with supporting input files.
This example uses the same bathymetry and vegetation inputs as the :doc:`tide_channel example <tide_channel>`.

This example consists of a simple, symmetrical beach slope with a channel that bisects the dune and connects to a lagoon.
A simple structured mesh is developed using Delft3D's `RGFGRID <https://content.oss.deltares.nl/delft3d4/RGFGRID_User_Manual.pdf>`_.
A tidal boundary condition is imposed on the left (ocean) side, that propagates through the channel and fills the lagoon behind the dune.
Two vegetation species, loosely based on `Nelumbo Lutea` and `Colocascia Esculenta`, are added to the model via a set of attributes defined in the files `NelumboLutea.json` and `ColocasciaEsculenta.json`.
After the simulation is finished, we can inspect a variety of 2-D model results using :class:`~dycove.utils.plotting.ModelPlotter`.

Running the model and plotting the results are handled mostly as described in the :doc:`tide_channel example <tide_channel>`.
The only change needed when modeling multiple species is in creating multiple :class:`~dycove.sim.vegetation_data.VegetationSpecies` objects:

.. code-block:: python

   >>> from dycove import MultipleVegetationSpecies as MultiVeg
   >>> veg_1 = VegetationSpecies("NelumboLutea.json", "nelumbo", mor=1, rand_seed_frac=0.8)
   >>> veg_2 = VegetationSpecies("ColocasciaEsculenta.json", "colocascia", mor=1, rand_seed_frac=0.8)
   >>> HydroModel = DFM_hydro.DFM(DFM_DLLs, config_file, mdu_file, vegetation=MultiVeg([veg_1, veg_2]))
   >>> HydroModel.run_simulation(3, sim_time_unit="eco-morphodynamic years")

In this example, the optional argument `rand_seed_frac` has also been changed from its default value.
By specifying `rand_seed_frac=0.8`, we are telling DYCOVE to only colonize vegetation of that species in 80 percent of eligible grid cells.

Plotting results of simulations with multiple species is no different from plotting results of a single-species simulation.
Plotting stem height using :class:`~dycove.utils.plotting.ModelPlotter`, for example, will average all stem heights in each grid cell among the various species and life stages present within that cell.
When plotting vegetation fractions, however, the output images will alternate between fractions of each species.
Example: A simulation of one ecological year and two species consists of one colonization event for each species in the second ecological time step (ETS) of the year.
Plots of stem heights will consist of one plot per ETS, starting at ETS 2.
Plots of fractions will consist of two plots per ETS, one for each species, starting at ETS 2.