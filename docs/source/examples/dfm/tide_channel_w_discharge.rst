.. _tide-channel-dfm-discharge:

Example - Simple beach and tide channel with discharge boundary and morphology (Delft3D FM)
===========================================================================================

This example corresponds to the Python file `tide_channel_w_discharge.py` located in the `examples/DFM/tide_channel_w_discharge/` directory, along with supporting input files.
This example uses the same bathymetry and vegetation inputs as the :doc:`tide_channel example <tide_channel>`.

This example consists of a simple, symmetrical beach slope with a channel that bisects the dune and connects to a lagoon.
A simple structured mesh is developed using Delft3D's `RGFGRID <https://content.oss.deltares.nl/delft3d4/RGFGRID_User_Manual.pdf>`_.
A tidal boundary condition is imposed on the left (ocean) side, that propagates through the channel and fills the lagoon behind the dune.
A vegetation species, loosely based on `Salicorniaspp.` (pioneer species, first life stage) and `Spartina anglica` (second life stage), is added to the model via a set of attributes defined in the file `veg1.json`.
After the simulation is finished, we can inspect a variety of 2-D model results using :class:`~dycove.utils.plotting.ModelPlotter`.

This example differs from the :doc:`tide_channel example <tide_channel>` in that morphology is active, and to promote substantial bed level changes in the domain, a discharge boundary has been added to the right side that discharges water and sediment to the lagoon.

Running the model and plotting the results are handled mostly as described in the :doc:`tide_channel example <tide_channel>`.
The only change needed (on the DYCOVE side) when morphology is turned on is to add ``mor = 1`` when instantiating :class:`~dycove.sim.vegetation_data.VegetationSpecies`:

.. code-block:: python

   veg_1 = VegetationSpecies("veg1.json", mor=1)
