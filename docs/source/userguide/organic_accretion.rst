.. _organic-accretion:

Organic Accretion
=================

DYCOVE can optionally calculate vegetation-driven organic-matter (OM)
accretion and add the resulting areal-mass increment to an existing
Delft3D FM sediment fraction through BMI. Organic accretion currently
requires both a DYCOVE vegetation model and the Delft3D FM engine.
Vegetation-only simulations remain unchanged when no
:class:`~dycove.sim.organic_accretion.OrganicAccretion` object is supplied.

All OM increments use units of ``kg m-2``. DYCOVE calculates an increment
once per ecological time step, using

.. math::

   \Delta t = \frac{1}{N_{ETS}} \quad \mathrm{yr},

where :math:`N_{ETS}` is ``n_ets``. The OM calculation therefore follows
the accelerated ecological clock rather than elapsed hydrodynamic seconds.

Activation
----------

Create an OM coupler and pass it to the Delft3D FM model:

.. code-block:: python

   from dycove import DFM_hydro, OrganicAccretion, VegetationSpecies

   vegetation = VegetationSpecies("veg1.json", mor=1)
   organic = OrganicAccretion(
       bed_type="underlayer",
       target_layer=0,
       target_sediment=2,
   )

   model = DFM_hydro.DFM(
       dfm_path,
       "dimr_config.xml",
       "dflowfm/FlowFM.mdu",
       vegetation=vegetation,
       organic=organic,
   )
   model.run_simulation(
       3,
       "eco-morphodynamic years",
       n_ets=12,
       veg_interval=43200,
       ecofac=30,
   )

Omit ``organic=organic`` to run the original vegetation-only coupling.
``target_sediment`` and ``target_layer`` are zero-based Python indices.
The target sediment must be specified explicitly to prevent accidental
deposition into the wrong sediment fraction.

Sediment storage
----------------

``bed_type="underlayer"`` updates Delft3D FM ``msed[cell, layer, sediment]``.
``target_layer=0`` selects the top underlayer. ``bed_type="mixing_layer"``
updates ``bodsed[cell, sediment]`` and ignores ``target_layer``.

The target must already exist as a Delft3D FM sediment fraction. Configure
the sediment and morphology files consistently before enabling OM. DYCOVE
supplies an areal dry-mass increment; Delft3D's ``CDryB`` value controls the
conversion between dry sediment mass and bed volume. Do not add
``rho_organic`` to the vegetation JSON because it is not used by the OM
production formulations.

For a standard serial DFM grid, DYCOVE uses the active-cell count reported
by the engine and maps vegetation cells to the corresponding prefix of the
BMI storage array. Extra boundary entries are left unchanged. For a custom
ordering, provide an explicitly verified integer ``cell_indices`` array with
one unique storage index per vegetation cell.

Calculation methods
-------------------

Select one of three methods in each JSON ``organic_accretion`` block.
All parameters listed for the selected method are required.

Rate
~~~~

The ``rate`` method applies a prescribed annual areal OM rate:

.. math::

   \Delta M = f_v R_{OM} \Delta t,

where :math:`f_v` is cohort fractional cover and :math:`R_{OM}` is
``om_rate_kg_m2_yr`` in ``kg m-2 yr-1``.

.. code-block:: json

   {
     "method": "rate",
     "om_rate_kg_m2_yr": 0.10
   }

Biomass proxy
~~~~~~~~~~~~~

The ``biomass_proxy`` method begins with a prescribed aboveground biomass
stock:

.. math::

   \Delta M = f_v B_{AGB} r k T \Delta t,

where :math:`B_{AGB}` is ``agb_proxy_kg_m2`` (``kg m-2``), :math:`r` is
``refractory_fraction``, :math:`k` is ``root_to_shoot``, and :math:`T` is
``turnover_rate`` (``yr-1``). The former name ``agb_proxy_kg_m2_yr`` is not
accepted because turnover already provides the inverse-time factor.

.. code-block:: json

   {
     "method": "biomass_proxy",
     "agb_proxy_kg_m2": 0.80,
     "refractory_fraction": 0.10,
     "root_to_shoot": 2.0,
     "turnover_rate": 2.0
   }

Biomass from vegetation geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``biomass_geometry`` method estimates cell-average aboveground biomass
from the current cohort geometry:

.. math::

   B_{AGB} = f_v n \frac{\pi d^2}{4} h \rho_p,

.. math::

   \Delta M = B_{AGB} r k T \Delta t,

where :math:`n`, :math:`d`, and :math:`h` are the cohort stem density
(``m-2``), diameter (``m``), and height (``m``), and :math:`\rho_p` is
``plant_bulk_density_kg_m3`` (``kg m-3``).

.. code-block:: json

   {
     "method": "biomass_geometry",
     "refractory_fraction": 0.10,
     "root_to_shoot": 2.0,
     "turnover_rate": 2.0,
     "plant_bulk_density_kg_m3": 100.0
   }

Species and life stages
-----------------------

OM parameters may be placed inside each entry of ``life_stage_attr`` so that
each vegetation life stage uses its own formulation and values:

.. code-block:: json

   {
     "nls": 2,
     "life_stage_attr": [
       {
         "organic_accretion": {
           "method": "biomass_geometry",
           "refractory_fraction": 0.10,
           "root_to_shoot": 2.0,
           "turnover_rate": 2.0,
           "plant_bulk_density_kg_m3": 100.0
         }
       },
       {
         "organic_accretion": {
           "method": "biomass_geometry",
           "refractory_fraction": 0.20,
           "root_to_shoot": 3.0,
           "turnover_rate": 2.0,
           "plant_bulk_density_kg_m3": 150.0
         }
       }
     ]
   }

The abbreviated JSON above shows only the OM blocks; retain all required
vegetation attributes in a real input file. Alternatively, one top-level
``organic_accretion`` block acts as a fallback for every life stage that does
not define its own block. In a multi-species run, each species reads parameters
from its own vegetation JSON file. Contributions from all cohorts and species
are added to the same explicitly selected sediment fraction. Total fractional
cover across all active cohorts and species must not exceed one in any cell.

Validation and diagnostics
--------------------------

DYCOVE rejects missing or nonfinite required parameters, negative rates and
biomass values, ``refractory_fraction`` outside zero to one, nonpositive plant
bulk density, invalid storage indices, and inconsistent grid lengths.

Set ``verbose=True`` for step-level summaries. For detailed verification of
one vegetation cell, also provide ``verification_cell=<index>``. Leave both at
their defaults for production simulations to avoid large log files.

Outputs and timing
------------------

Organic outputs are written under ``dflowfm/om_output``:

* ``om_summary.csv`` contains species-level summaries for each ecological step.
* ``om_yearYYYY_etsEEE.nc`` contains ``fveg``, ``agb_estimate``, and
  ``organic_increment`` for every internal vegetation cell. Parallel runs use
  rank-specific suffixes to avoid file collisions.

The summary reports both the increment for each species and the cumulative
increment applied since the start of the simulation. The NetCDF global
attributes record the storage variable, target sediment, and target layer
(``not_applicable`` for ``bodsed``).

DYCOVE applies each OM increment after the corresponding hydrodynamic advance
and vegetation update. Consequently, an increment is visible in the standard
DFM map output after the next hydrodynamic advance. The final OM update is
preserved in ``om_output`` but may not appear in the final scheduled DFM map
record when the simulation ends immediately after that coupling step.
