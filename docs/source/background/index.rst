Background
==========

Motivation
----------

DYCOVE originated as a MATLAB code coupled with Delft3D 4 [1]_.
In that study by Brückner et al. (2019), the dynamic vegetation approach was shown to reproduce observed vegetation distributions in salt marshes compared to static approaches.
With DYCOVE written in Python, we can couple to the widely-used and modernized hydro-morphodynamic model `Delft3D FM (DFM) <https://oss.deltares.nl/web/delft3dfm>`_, while also offering a fully open-source avenue with the `ANUGA hydrodynamic model <https://github.com/GeoscienceAustralia/anuga_core>`_.
Furthermore, DYCOVE is written in such a way that will allow for straightforward integration of other hydro-morphodynamic models in the future.


Model Overview
--------------

General
#######

DYCOVE allows users to incorporate vegetation in their numerical models in a dynamic way. 
A set of vegetation species input parameters must be provided that describe species attributes such as initial and maximum stem height and diameter, number of life stages, colonization window, and more.







Vegetation Coupling Processes
#############################

As mentioned in the previous section, vegetation is coupled to the numerical model once every ``veg_interval`` seconds.
During each coupling, the following steps take place before the hydrodynamics continue:

1. **Life-stage update**

   * Increment ecological counters (``ets``, ``eco_year``, etc.)
   * Apply ecological year and life-stage transitions

2. **Growth**

   * Apply height, diameter, and root growth based on linear growth rates computed from vegetation inputs.

3. **Colonization**

   * Establish new vegetation fractions where conditions satisfy:

     * *Dry* at some point during previous ``ets``,
     * *Wet* at some point during previous ``ets``,
     * Random or deterministic seed fraction filtering (optional),
     * Available space in the cell

4. **Mortality**

   * Compute hydrodynamic mortality (flooding, desiccation, uprooting),
   * Compute morphodynamic mortality (burial, scour, if enabled),
   * Reduce vegetation fractions according to mortality fractions

5. **Vegetation model quantity computation**

   * Stem density: compute area-weighted sum of cohort densities × fractions
   * Stem height/diameter: compute area-weighted averages based on vegetation fractions

6. **Injection back into hydrodynamics**

   * Update numerical model with new arrays of stem density, diameter, and height



Colonization
############







### **Hydrodynamic Time Loop**

Hydrodynamics advance in short increments (`hydro_interval`, e.g., 900 seconds).
Within each ecological timestep, the engine performs:

1. **Hydrodynamic stepping** via :meth:`HydroEngineBase.step`.
2. **Extraction of velocity and depth** fields at each substep.
3. **Accumulation of statistics** stored in :class:`HydrodynamicStats`:

   * Minimum and maximum water depth
   * Maximum velocity
   * Flooded/dry timestep counts
   * (Optional) bed‐level change if morphology is active

These aggregated statistics represent the environmental stress experienced by vegetation during that ETS.

---

### **2. Vegetation Update Loop**

The vegetation model updates once per ecological timestep (`veg_interval`, e.g., 12 hours).
At each vegetation step, the coupler performs:


---

## **Time Scaling and Ecological Time Steps**

Vegetation operates on a **slow time scale**, while hydrodynamics update quickly.
To reconcile these, the model uses two key mechanisms:

### **Ecological Time Step (ETS)**

* One ETS represents a discrete vegetation-update period.
* Total vegetation years simulated = `sim_time` (in years) or converted from hydrodynamic days.
* Number of ETS per year = `n_ets` (default: 14).

### **Ecological Acceleration Factor (ecofac)**

The ecological acceleration factor links hydrodynamic time (seconds) to ecological time (years). It determines how fast the “environmental clock” runs relative to real-time seconds of hydrodynamics.

The logic is:

1. If morphology is active, `ecofac = MORFAC`.
2. If the user provides an `ecofac`, it must match `MORFAC` (if morphology exists).
3. Otherwise, the model computes:

   [
   \text{ecofac} \approx \frac{365 \times 86400}{\text{veg_interval} \times n_{\text{ets}}}
   ]

This ensures that vegetation years remain physically meaningful even when using accelerated hydrodynamics.

### **Days Per Ecological Year**

From `ecofac`, `veg_interval`, and `n_ets`, the model derives an “implied days per year,” which must fall within a realistic range (350–380 days). This acts as a consistency check for input time parameters.

---

## **Hydrodynamic Statistics Used by Vegetation**

Vegetation response is controlled by several aggregated hydrodynamic metrics:

* **Minimum depth**: used to detect dry periods (desiccation stress).
* **Maximum depth**: used to detect flooding events.
* **Maximum velocity**: used for uprooting mortality calculations.
* **Flooded/dry fractions**: ratio of hydrodynamic substeps meeting wet/dry criteria.
* **Bed-level change** *(only if morphology enabled)*:

  * Positive values indicate scour,
  * Negative values indicate deposition/burial.

These variables provide a robust, averaged representation of flow conditions for the vegetation lifecycle processes.

---

## **Vegetation Parameter Inputs**

Vegetation traits and life-cycle settings are defined in JSON input files and loaded into a :class:`VegetationAttributes` dataclass. Key parameters include:

### **1. Geometric and Functional Traits**

* Initial shoot height, root length, stem diameter
* Maximum plant size per life stage
* Linear seasonal growth-rate parameters (computed internally)

### **2. Life-Stage Structure**

* Number of life stages (`nls`)
* Maximum years spent in each stage
* Stem density per stage
* Winter retreat height
* Growth window: `start_growth_ets`, `end_growth_ets`, `winter_ets`

### **3. Colonization Parameters**

* Initial colonization fraction (`fraction_0`)
* Colonization period (`start_col_ets` → `end_col_ets`)
* Random/deterministic seed distribution
* Maximum fraction limited by available space in each cell

### **4. Mortality Thresholds**

Hydrodynamic mortality is linear between two thresholds per life stage:

* Flooding: `flood_no_mort`, `flood_all_mort`
* Desiccation: `desic_no_mort`, `desic_all_mort`
* Uprooting: `uproot_no_mort`, `uproot_all_mort`

Morphodynamic mortality (if enabled):

* Burial fraction (relative to vegetation height)
* Scour fraction (relative to root length)

---

## **Multi-Species Capability**

Multiple species can be combined using :class:`MultipleVegetationSpecies`, which:

* Manages all species jointly,
* Ensures consistent use of morphology flag (`mor`),
* Flattens all cohorts into a unified set for weighted averaging,
* Passes inter-species information during colonization to handle space competition.

---

## **Summary**

The DYCOVE vegetation coupling system provides:

* A clear separation between hydrodynamic and vegetation logic,
* Flexible, engine-agnostic hydrodynamic interaction,
* Robust treatment of growth, colonization, and mortality through ETS-based updates,
* Automatic and user-controlled ecological time scaling,
* Support for single and multiple species,
* Parallel-safe vegetation output and merging.

This structure allows users to incorporate vegetation in morphodynamic simulations ranging from small idealized problems to large-scale, fully parallelized riverine or deltaic systems.

---

If you’d like, I can also generate:

✅ A diagram of the coupling loop
✅ A simple “Quick Start” usage example
✅ A glossary of vegetation terms
✅ A version tailored for Sphinx `.rst` with directives

Just let me know!
