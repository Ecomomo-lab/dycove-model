User Guide
==========

This user guide provides explanations for DYCOVE concepts.
While users can get started quickly by running the examples, advanced use of DYCOVE requires an understanding of these concepts.


The Ecological Time Scale
#########################

DYCOVE simulations run on the idea that vegetation and morphological changes both occur on longer time scales than hydrodynamic changes.
Therefore, we "accelerate" vegetation processes in a manner similar to the Morphological Acceleration Factor (MORFAC) of DFM.
In fact, for DYCOVE-DFM simulations with morphology active, we set the Ecological Acceleration Factor (ECOFAC, ``ecofac``) equal to the MORFAC of the DFM model.

To reconcile these two time scales, we run simulations using an Ecological Time Step, or ``ets``, which represents the vegetation coupling interval ``veg_interval`` in the model.
And because DYCOVE is most useful for modeling vegetation in the intertidal zone, it is convenient to structure simulations around a full tidal cycle.
Therefore, we often set ``veg_interval`` to be equal to one tidal cycle. 
For your typical semi-diurnal tide, one cycle is equal to about 12 hours, rounding down for convenience.
We must decide, then, whether we want to prescribe the ``ecofac`` of our model or the number of ``ets```` per year ``n_ets``, and then the other variable will follow:

   [
   \text{ecofac} \approx \frac{365 \times 86400}{\text{veg_interval} \times \text{n}_{\text{ets}}}
   ]

Standard values of ``veg_interval = 43200`` and ``n_ets = 14`` will yield ``ecofac = 52``. 
In a similar scenario, users can decide to use a more round value of ``ecofac = 50``, which would correspond to 350 days per year.
It is fine, and in fact encouraged, to use a round value in this way, but keep in mind that DYCOVE includes a hardcoded limit of ``DAYS_PER_YEAR`` between 350 and 380 to keep results reasonable.

Note that DYCOVE simulations can always be run using the default values provided in :meth:`~dycove.sim.base.run_simulation`.


The Vegetation Input File
#########################



The Vegetation Input File
#########################