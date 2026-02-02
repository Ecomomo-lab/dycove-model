.. _tide-channel-dfm:

Example - Simple beach and tide channel (Delft3D FM)
====================================================

This example corresponds to the Python file `tide_channel.py` located in the `examples/DFM/tide_channel/` directory, along with supporting input files.
This example uses the same bathymetry and vegetation inputs as the corresponding :doc:`ANUGA example <../anuga/tide_channel>`.
That example includes some more discussion on DYCOVE's vegetation logic and inputs than this example does, and users are always encouraged to read the `background documentation <https://Ecomomo-lab.github.io/dycove-model/background/ecological_time_scaling.html>`_ for full exposure to DYCOVE concepts.

This example consists of a simple, symmetrical beach slope with a channel that bisects the dune and connects to a lagoon.
A simple structured mesh is developed using Delft3D's `RGFGRID <https://content.oss.deltares.nl/delft3d4/RGFGRID_User_Manual.pdf>`_, one that mimics the mesh developed for the ANUGA version of this example.
A tidal boundary condition is imposed on the left (ocean) side, that propagates through the channel and fills the lagoon behind the dune.
A vegetation species, loosely based on `Salicorniaspp.` (pioneer species, first life stage) and `Spartina anglica` (second life stage), is added to the model via a set of attributes defined in the file `veg1.json`.
After the simulation is finished, we can inspect a variety of 2-D model results using :class:`~dycove.utils.plotting.ModelPlotter`.

This example does not consider morphology. Please see this :doc:`example <tide_channel_w_discharge>` for a model that has morphology turned on.


Run the Model
-------------

The following import statement forms the basis of all DYCOVE-DFM coupling simulations:

.. code-block:: python

   from dycove import VegetationSpecies, DFM_hydro

To run Delft3D FM, you will need to create and export your model outside of the coding environment.
After creating a Delft3D model using Delta Shell, export the model as DIMR:

  - Alt-click on the model name in the Project tree
  - "Export..." → "DIMR configuration"

This will create an `.xml` file at the same level as a directory called `dflowfm/`, which contains the `.mdu` file inside.
Next, identify the path to your Delft3D FM code library, which probably looks something like the path shown below:

.. code-block:: python

   from pathlib import Path
   config_file = 'dimr_config.xml'
   mdu_file = Path('dflowfm/FlowFM.mdu')
   DFM_DLLs = Path('C:/Program Files (x86)/Deltares/Delft3D Flexible Mesh Suite HM (2021.03)/plugins/DeltaShell.Dimr/kernels/x64')

To turn on the internal vegetation module in Delft3D FM, we must create an additional forcing `.ext` file and add a new block in the `.mdu` file with the `[veg]` heading.
DYCOVE handles this automatically, as long as you supply a vegetation species to your DYCOVE model.

After setting up the model and identifying the model file paths, DYCOVE-DFM models are run in much the same way as DYCOVE-ANUGA models.
The only difference is the call to the hydrodynamic engine, which in this case is :class:`~dycove.sim.engines.DFM_hydro.DFM`:

.. code-block:: python

   veg_1 = VegetationSpecies("veg1.json")
   HydroModel = DFM_hydro.DFM(DFM_DLLs, config_file, mdu_file, vegetation=veg_1)
   HydroModel.run_simulation(3, sim_time_unit="eco-morphodynamic years")


Plot the Results
----------------

DFM map output is saved in the `dflowfm/output/FlowFM_map.nc` file (or whatever the model name is, if not `FlowFM`). 
DYCOVE output is saved in `.npz` files, with one file per vegetation cohort per ecological time step.
For a list of outputs stored with each cohort, see :class:`~dycove.sim.vegetation_data.VegCohort`.

Plot 2-D maps of various model quantities over time using :class:`~dycove.utils.plotting.ModelPlotter`.
For example, we can plot vegetation stem heights over the entire 21-day (hydrodynamic time) simulation period.

.. code-block:: python

   from dycove import plotting
   plotter = plotting.ModelPlotter(
      simdir = Path("."),
      #quantity = "Velocity",
      quantity = "Stem Height",
      plot_times = {
         "plotHR_0": 0*24.,
         "plotHR_f": 28*24.,
         "mapHR_int": 1,
         "plotHR_int": 1,
      },
      cmap_lims = {
         "Bathymetry": (-0.5, 0.5),
         "Velocity": (0, 0.3),
      },
      )
   plotter.run()

This plotting code is located in the same example directory at `examples/DFM/tide_channel/plot_tide_channel_DFM.py`, where there are additional comments on how to use the code, as well as more ideas for quantities to plot.
