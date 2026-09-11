.. _organic-accretion-dfm-example:

Example - Tide channel with discharge and organic accretion (Delft3D FM)
===========================================================================

This standalone example is located in
``examples/DFM/tide_channel_w_discharge_organic_accretion``. It extends the
single-species :doc:`tide-channel-with-discharge example
<tide_channel_w_discharge>` with vegetation-driven organic accretion. The
directory contains the Python/DIMR inputs and a complete DeltaShell project.
See the :ref:`organic accretion user guide <organic-accretion>` for the
equations, parameter definitions, validation rules, and output conventions.

Model configuration
-------------------

The example retains the original mixing-layer morphology. ``FlowFM.sed``
contains two sediment fractions:

* ``sedimentsand`` at zero-based sediment index 0;
* ``sedimentorganic`` at zero-based sediment index 1.

The organic fraction is represented as a mud-style sediment and begins with
``IniSedThick=0``. DYCOVE supplies its areal dry-mass increment through BMI.
The organic sediment ``CDryB`` value controls the mass-to-bed-volume
conversion. The sand boundary condition is unchanged, and no organic
sediment boundary concentration is imposed.

``veg1.json`` contains one top-level ``organic_accretion`` block. Its values
are illustrative rather than species-calibrated and act as the fallback for
both vegetation life stages:

.. code-block:: json

   "organic_accretion": {
     "method": "biomass_geometry",
     "refractory_fraction": 0.10,
     "root_to_shoot": 2.0,
     "turnover_rate": 2.0,
     "plant_bulk_density_kg_m3": 100.0
   }

Running the example
-------------------

Navigate to the example directory and run:

.. code-block:: console

   python run_tide_channel_organic_accretion.py

The central additions in that script are:

.. code-block:: python

   from dycove import DFM_hydro, OrganicAccretion, VegetationSpecies

   vegetation = VegetationSpecies(work_dir / "veg1.json", mor=1)
   organic = OrganicAccretion(
       bed_type="mixing_layer",
       target_sediment=1,
   )

   model = DFM_hydro.DFM(
       DFM_DLLs,
       config_file,
       mdu_file,
       vegetation=vegetation,
       organic=organic,
   )
   model.run_simulation(
       4,
       "eco-morphodynamic years",
       n_ets=14,
       veg_interval=43200,
   )

Because morphology is active, DYCOVE reads ``MorFac=52`` from ``FlowFM.mor``.
The script does not supply a conflicting ``ecofac`` value. Four ecological
years with 14 ecological time steps per year produce 56 OM coupling updates
and require 28 hydrodynamic model days for this configuration.

Outputs
-------

The usual DFM results are written under ``dflowfm/output`` and vegetation
results under ``dflowfm/veg_output``. Organic-accretion results are written
under ``dflowfm/om_output``:

* ``om_summary.csv`` stores one species-level summary per ecological step;
* ``om_yearYYYY_etsEEE.nc`` stores ``fveg``, ``agb_estimate``, and
  ``organic_increment`` for every internal cell.

The example has been verified to begin with zero organic sediment, create 56
nonzero OM-step files, and produce positive organic storage in
``mesh2d_bodsed[:, :, 1]``. The final ecological update is recorded in the
dedicated OM output but may not appear in the last DFM map record because it
occurs after the final hydrodynamic advance.

The copied ``FlowFM.sed`` under the DeltaShell ``.dsproj_data`` directory
contains the same two sediment definitions, so users can inspect the complete
model in DeltaShell as well as run it through Python.
