###############################################################
#  organic_accretion.py
###############################################################

import numpy as np


class OrganicAccretion:
    """
    Optional DYCOVE organic accretion coupler.

    Supports:
    - bed_type="underlayer"   -> sediment_storage="msed"
    - bed_type="mixing_layer" -> sediment_storage="bodsed"
    - single-species and multi-species vegetation
    - species-specific OM parameters
    - life-stage-specific OM parameters

    Organic matter increments are computed in kg/m2.
    """

    def __init__(
        self,
        method=None,
        om_rate_kg_m2_yr=None,
        agb_proxy_kg_m2=None,
        refractory_fraction=None,
        root_to_shoot=None,
        turnover_rate=None,
        plant_bulk_density_kg_m3=None,
        target_layer=0,
        target_sediment=None,
        bed_type=None,
        sediment_storage=None,
        cell_indices=None,
        verification_cell=None,
        verbose=False,
    ):
        self.method = method
        self.om_rate_kg_m2_yr = om_rate_kg_m2_yr
        self.agb_proxy_kg_m2 = agb_proxy_kg_m2
        self.refractory_fraction = refractory_fraction
        self.root_to_shoot = root_to_shoot
        self.turnover_rate = turnover_rate
        self.plant_bulk_density_kg_m3 = plant_bulk_density_kg_m3

        self.target_layer = target_layer
        self.target_sediment = target_sediment
        self.bed_type = bed_type
        self.verification_cell = verification_cell
        self.verbose = verbose
        self.cell_indices = None if cell_indices is None else np.asarray(cell_indices)
        self.cumulative_added_mass = 0.0

        self.allowed_methods = [
            "rate",
            "biomass_proxy",
            "biomass_geometry",
        ]

        if self.method is not None and self.method not in self.allowed_methods:
            raise ValueError(
                f"Unknown organic accretion method: {self.method}. "
                f"Allowed methods are: {self.allowed_methods}"
            )

        self.sediment_storage = self.resolve_sediment_storage(
            bed_type=bed_type,
            sediment_storage=sediment_storage,
        )

        allowed_storage = ["msed", "bodsed"]

        if self.sediment_storage not in allowed_storage:
            raise ValueError(
                f"Unknown sediment_storage: {self.sediment_storage}. "
                f"Allowed values are: {allowed_storage}"
            )

        if self.verbose and self.bed_type == "mixing_layer" and target_layer != 0:
            print(
                "OrganicAccretion: target_layer is ignored when "
                "bed_type='mixing_layer'."
            )

    def resolve_sediment_storage(
        self,
        bed_type=None,
        sediment_storage=None,
    ):
        """
        Resolve user-facing bed_type into the internal Delft3D-FM
        BMI sediment-storage variable.

        bed_type="underlayer"
            -> msed[cell, layer, sediment]

        bed_type="mixing_layer"
            -> bodsed[cell, sediment]

        sediment_storage is retained for backward compatibility.
        """

        if sediment_storage is not None:
            return sediment_storage

        if bed_type is None:
            return "msed"

        allowed_bed_types = {
            "underlayer": "msed",
            "mixing_layer": "bodsed",
            "mixinglayer": "bodsed",
            "bed": "bodsed",
        }

        if bed_type not in allowed_bed_types:
            raise ValueError(
                f"Unknown bed_type: {bed_type}. "
                f"Allowed values are: {list(allowed_bed_types.keys())}"
            )

        return allowed_bed_types[bed_type]

    def get_dt_years(self, simstate):
        """
        Return the vegetation coupling interval in years.
        """

        if simstate is None:
            raise ValueError(
                "simstate is required for time-dependent "
                "organic accretion methods."
            )

        if simstate.n_ets <= 0:
            raise ValueError("simstate.n_ets must be greater than zero.")

        # Vegetation geometry advances once per ecological time step. OM must
        # use that same accelerated clock, independent of hydrodynamic seconds.
        return 1.0 / simstate.n_ets

    def compute_total_fveg(self, cohorts):
        """
        Compute total vegetation cover from all cohorts of one species.
        """

        if len(cohorts) == 0:
            return None

        fractions = [cohort.fraction for cohort in cohorts]
        fveg = np.sum(fractions, axis=0)

        return np.clip(fveg, 0.0, 1.0)

    def compute_geometry_biomass_estimate(self, cohorts):
        """
        Compute geometry-based biomass using the class-level default
        plant bulk density.

        This method is retained for backward compatibility and general
        diagnostics. Life-stage calculations use compute_cohort_om().
        """

        if len(cohorts) == 0:
            return None

        biomass = None

        for cohort in cohorts:
            stem_area = np.pi * (cohort.diameter ** 2) / 4.0

            volume_index = (
                cohort.fraction
                * cohort.density
                * stem_area
                * cohort.height
            )

            cohort_biomass = (
                volume_index
                * self.plant_bulk_density_kg_m3
            )

            if biomass is None:
                biomass = cohort_biomass.copy()
            else:
                biomass += cohort_biomass

        return biomass

    def get_life_stage_params(self, species, cohort):
        """
        Return the OM parameter dictionary corresponding to the
        cohort's current vegetation life stage.

        Vegetation-only simulations may use original JSON files without
        OM parameters. However, if OrganicAccretion is activated, the
        active cohort must have a valid OM parameter dictionary.
        """

        life_stage = int(cohort.lifestage)

        if hasattr(
            species,
            "organic_accretion_params_by_lifestage",
        ):
            params_by_stage = (
                species.organic_accretion_params_by_lifestage
            )

            stage_index = life_stage - 1

            if stage_index < 0 or stage_index >= len(params_by_stage):
                raise IndexError(
                    "OrganicAccretion: cohort life stage is outside "
                    "the available OM parameter range. "
                    f"species={species.name}, "
                    f"cohort_lifestage={life_stage}, "
                    f"available_stages={len(params_by_stage)}"
                )

            stage_params = params_by_stage[stage_index]

            if not stage_params:
                raise ValueError(
                    "OrganicAccretion is enabled, but no "
                    "organic_accretion parameters were provided for "
                    f"species='{species.name}', life_stage={life_stage}. "
                    "Add an organic_accretion block to this life stage "
                    "in the vegetation JSON, or disable OrganicAccretion."
                )

            return stage_params

        if hasattr(species, "organic_accretion_params"):
            species_params = species.organic_accretion_params

            if not species_params:
                raise ValueError(
                    "OrganicAccretion is enabled, but no "
                    "organic_accretion parameters were provided for "
                    f"species='{species.name}'. Add OM parameters to "
                    "the vegetation JSON, or disable OrganicAccretion."
                )

            return species_params

        raise RuntimeError(
            "OrganicAccretion: vegetation species has neither "
            "'organic_accretion_params_by_lifestage' nor "
            "'organic_accretion_params'."
        )

    def compute_cohort_om(
        self,
        cohort,
        stage_params,
        simstate,
    ):
        """
        Compute one cohort's organic matter contribution using the
        parameters associated with that cohort's current life stage.

        Parameters
        ----------
        cohort
            DYCOVE VegCohort instance.

        stage_params : dict
            Organic matter parameters for the current life stage.

        simstate
            DYCOVE simulation-time state.

        Returns
        -------
        delta_msed : numpy.ndarray
            Organic matter increment per vegetation cell, kg/m2.

        agb_estimate : numpy.ndarray or None
            Geometry-based aboveground biomass estimate, kg/m2.

        method : str
            Organic matter formulation used for this cohort.
        """

        if not isinstance(stage_params, dict):
            raise TypeError("Organic accretion parameters must be a dictionary.")

        method = stage_params.get("method", self.method)

        if method is None:
            raise ValueError(
                "Organic accretion method must be provided in the vegetation "
                "JSON or OrganicAccretion constructor."
            )

        if "rho_organic" in stage_params:
            raise ValueError(
                "rho_organic is not an OM-production parameter. Remove it from "
                "the vegetation JSON; Delft3D bed-volume conversion is controlled "
                "by CDryB in the sediment file."
            )

        if method not in self.allowed_methods:
            raise ValueError(
                f"Unknown organic accretion method '{method}' "
                f"for cohort life stage {cohort.lifestage}. "
                f"Allowed methods are: {self.allowed_methods}"
            )

        om_rate_kg_m2_yr = stage_params.get(
            "om_rate_kg_m2_yr",
            self.om_rate_kg_m2_yr,
        )

        if "agb_proxy_kg_m2_yr" in stage_params:
            raise ValueError(
                "agb_proxy_kg_m2_yr has inconsistent units. Rename it to "
                "agb_proxy_kg_m2; turnover_rate supplies the per-year factor."
            )

        agb_proxy_kg_m2 = stage_params.get(
            "agb_proxy_kg_m2",
            self.agb_proxy_kg_m2,
        )

        refractory_fraction = stage_params.get(
            "refractory_fraction",
            self.refractory_fraction,
        )

        root_to_shoot = stage_params.get(
            "root_to_shoot",
            self.root_to_shoot,
        )

        turnover_rate = stage_params.get(
            "turnover_rate",
            self.turnover_rate,
        )

        plant_bulk_density = stage_params.get(
            "plant_bulk_density_kg_m3",
            self.plant_bulk_density_kg_m3,
        )

        required = {
            "rate": {"om_rate_kg_m2_yr": om_rate_kg_m2_yr},
            "biomass_proxy": {
                "agb_proxy_kg_m2": agb_proxy_kg_m2,
                "refractory_fraction": refractory_fraction,
                "root_to_shoot": root_to_shoot,
                "turnover_rate": turnover_rate,
            },
            "biomass_geometry": {
                "refractory_fraction": refractory_fraction,
                "root_to_shoot": root_to_shoot,
                "turnover_rate": turnover_rate,
                "plant_bulk_density_kg_m3": plant_bulk_density,
            },
        }[method]
        missing = [name for name, value in required.items() if value is None]
        if missing:
            raise ValueError(
                f"Method '{method}' requires parameter(s): {', '.join(missing)}."
            )
        for name, value in required.items():
            if not np.isscalar(value) or not np.isfinite(value):
                raise ValueError(f"{name} must be a finite scalar.")
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative.")
        if refractory_fraction is not None and not 0.0 <= refractory_fraction <= 1.0:
            raise ValueError("refractory_fraction must be between zero and one.")
        if method == "biomass_geometry" and plant_bulk_density <= 0.0:
            raise ValueError("plant_bulk_density_kg_m3 must be greater than zero.")

        for name in ("density", "diameter", "height"):
            values = np.asarray(getattr(cohort, name))
            if not np.all(np.isfinite(values)) or np.any(values < 0.0):
                raise ValueError(f"Cohort {name} must be finite and nonnegative.")

        fraction = np.clip(
            cohort.fraction,
            0.0,
            1.0,
        )

        agb_estimate = None

        if method == "rate":
            dt_years = self.get_dt_years(simstate)

            delta_msed = (
                om_rate_kg_m2_yr
                * fraction
                * dt_years
            )

        elif method == "biomass_proxy":
            dt_years = self.get_dt_years(simstate)

            retained_om_rate = (
                agb_proxy_kg_m2
                * refractory_fraction
                * root_to_shoot
                * turnover_rate
            )

            delta_msed = (
                fraction
                * retained_om_rate
                * dt_years
            )

        elif method == "biomass_geometry":
            dt_years = self.get_dt_years(simstate)

            stem_area = (
                np.pi
                * cohort.diameter ** 2
                / 4.0
            )

            volume_index = (
                fraction
                * cohort.density
                * stem_area
                * cohort.height
            )

            agb_estimate = (
                volume_index
                * plant_bulk_density
            )

            retained_om_rate = (
                agb_estimate
                * refractory_fraction
                * root_to_shoot
                * turnover_rate
            )

            delta_msed = (
                retained_om_rate
                * dt_years
            )

        else:
            raise RuntimeError(
                f"Unhandled organic accretion method: {method}"
            )

        if not np.all(np.isfinite(delta_msed)) or np.any(delta_msed < 0.0):
            raise ValueError("Calculated organic increment must be finite and nonnegative.")

        return delta_msed, agb_estimate, method

    def get_storage_indices(
        self,
        storage,
        delta_msed,
        storage_name,
        n_active_cells=None,
    ):
        """
        Return sediment-storage indices corresponding to vegetation cells.
        """

        n_om = len(delta_msed)
        n_storage = storage.shape[0]

        if self.cell_indices is not None:
            indices = np.asarray(self.cell_indices)
            if indices.ndim != 1 or len(indices) != n_om:
                raise ValueError(
                    "OrganicAccretion: cell_indices must be a one-dimensional "
                    "array with one index per vegetation cell."
                )
            if np.issubdtype(indices.dtype, np.bool_) or not np.issubdtype(
                indices.dtype, np.integer
            ):
                raise TypeError("OrganicAccretion: cell_indices must contain integers.")
            if len(np.unique(indices)) != len(indices):
                raise ValueError("OrganicAccretion: cell_indices must be unique.")
            if np.any(indices < 0) or np.any(indices >= n_storage):
                raise IndexError("OrganicAccretion: cell_indices are outside storage bounds.")
            return indices

        if n_storage == n_om:
            if not self.verbose:
                return np.arange(n_om)
            print(
                f"OrganicAccretion: {storage_name} and vegetation "
                f"grid match exactly: {n_om} cells."
            )

            return np.arange(n_om)

        if n_storage > n_om and n_active_cells == n_om:
            n_extra = n_storage - n_om

            if self.verbose:
                print(
                    f"OrganicAccretion: {storage_name} has "
                    f"{n_storage} entries, vegetation/cohort grid has "
                    f"{n_om} cells. Applying OM to first {n_om} cells "
                    f"only; leaving {n_extra} extra storage entries "
                    "unchanged."
                )

            return np.arange(n_om)

        if n_storage > n_om:
            raise ValueError(
                f"OrganicAccretion: {storage_name} has {n_storage} entries but "
                f"the vegetation grid has {n_om}. Provide verified cell_indices "
                "or ensure engine.get_cell_count() reports the vegetation-grid "
                "length (DFM ndxi)."
            )

        raise ValueError(
            f"OrganicAccretion: {storage_name} has fewer entries "
            f"than the vegetation grid: storage={n_storage}, "
            f"vegetation={n_om}. Cannot safely apply OM."
        )

    def print_life_stage_verification(
        self,
        species,
        cohort_index,
        cohort,
        stage_params,
        method,
        agb_estimate,
        delta_msed,
        simstate,
    ):
        """
        Print a focused verification record for one cohort and one
        selected grid cell.
        """

        cell = self.verification_cell

        if cell is None or cell < 0 or cell >= len(delta_msed):
            return

        print("\n========== LIFE-STAGE OM VERIFICATION ==========")
        print(f"Species: {species.name}")
        print(f"Cohort index: {cohort_index}")
        print(f"Life stage: {cohort.lifestage}")
        print(f"Method: {method}")
        print(f"Vegetation interval (s): {simstate.veg_interval}")
        print(f"Vegetation interval (yr): {self.get_dt_years(simstate)}")

        print("\nLife-stage OM parameters:")

        for key, value in stage_params.items():
            print(f"{key} = {value}")

        print(f"\nCell {cell}")
        print(f"fraction = {cohort.fraction[cell]}")
        print(f"density = {cohort.density}")
        print(f"diameter = {cohort.diameter}")
        print(f"height = {cohort.height}")

        if agb_estimate is not None:
            print(f"agb_estimate = {agb_estimate[cell]}")

        print(f"organic_increment = {delta_msed[cell]}")
        print("================================================\n")

    def update(self, engine, simstate=None):
        """
        Compute OM contributions for every species and cohort, select
        life-stage-specific OM parameters, and update Delft3D-FM bed
        sediment storage through BMI.

        Returns
        -------
        organic_records : list[dict]
            Species-level summary records for CSV output.

        organic_cell_records : list[dict]
            Species-by-cell arrays for NetCDF output.
        """

        storage_name = self.sediment_storage

        organic_records = []
        organic_cell_records = []

        storage_raw = engine.dflowfm.get_var(storage_name)

        if storage_raw is None:
            raise RuntimeError(
                "OrganicAccretion: BMI variable "
                f"'{storage_name}' is not available."
            )

        storage = storage_raw.copy()

        if storage_name == "msed":
            if storage.ndim != 3:
                raise ValueError(
                    f"OrganicAccretion: msed must be 3-D; got shape {storage.shape}."
                )
            if not isinstance(self.target_layer, (int, np.integer)):
                raise TypeError("target_layer must be an integer for msed storage.")
            if not 0 <= self.target_layer < storage.shape[1]:
                raise IndexError(
                    f"target_layer={self.target_layer} is outside msed layer "
                    f"bounds [0, {storage.shape[1] - 1}]."
                )
        elif storage_name == "bodsed" and storage.ndim != 2:
            raise ValueError(
                f"OrganicAccretion: bodsed must be 2-D; got shape {storage.shape}."
            )

        if self.target_sediment is None:
            raise ValueError(
                "target_sediment must be explicitly specified when organic "
                "accretion is enabled."
            )
        if not isinstance(self.target_sediment, (int, np.integer)):
            raise TypeError("target_sediment must be an integer.")
        if not 0 <= self.target_sediment < storage.shape[-1]:
            raise IndexError(
                f"target_sediment={self.target_sediment} is outside sediment "
                f"bounds [0, {storage.shape[-1] - 1}]."
            )

        if storage_name == "msed":
            before_total = np.nansum(
                storage[:, :, self.target_sediment]
            )

        elif storage_name == "bodsed":
            before_total = np.nansum(
                storage[:, self.target_sediment]
            )

        else:
            raise RuntimeError(
                f"Unhandled sediment_storage: {storage_name}"
            )

        if hasattr(engine.veg, "species_list"):
            species_units = engine.veg.species_list

        elif hasattr(engine.veg, "cohorts"):
            species_units = [engine.veg]

        else:
            raise RuntimeError(
                "OrganicAccretion: vegetation object has neither "
                "'species_list' nor 'cohorts'."
            )

        total_fraction = None
        expected_shape = None
        for species in species_units:
            for cohort in species.cohorts:
                fraction = np.asarray(cohort.fraction)
                if fraction.ndim != 1 or not np.all(np.isfinite(fraction)):
                    raise ValueError("Cohort fractions must be finite one-dimensional arrays.")
                if np.any(fraction < 0.0) or np.any(fraction > 1.0):
                    raise ValueError("Every cohort fraction must be between zero and one.")
                if expected_shape is None:
                    expected_shape = fraction.shape
                    total_fraction = np.zeros_like(fraction, dtype=float)
                elif fraction.shape != expected_shape:
                    raise ValueError("All cohort fraction arrays must have the same shape.")
                total_fraction += fraction

        if total_fraction is not None and np.any(total_fraction > 1.0 + 1e-12):
            raise ValueError(
                "Summed vegetation cover across all species and cohorts exceeds one."
            )

        n_active_cells = engine.get_cell_count()
        if expected_shape is not None and expected_shape[0] != n_active_cells:
            raise ValueError(
                "Vegetation grid length does not match the engine's active-cell "
                f"count: vegetation={expected_shape[0]}, active={n_active_cells}."
            )

        total_added_mass = 0.0

        for species in species_units:
            cohorts = species.cohorts

            if len(cohorts) == 0:
                if self.verbose:
                    print(
                        f"OrganicAccretion: species={species.name}, "
                        "no active cohorts; no OM added."
                    )

                continue

            fveg = self.compute_total_fveg(cohorts)

            species_delta_msed = np.zeros_like(
                fveg,
                dtype=float,
            )

            species_agb_estimate = None
            methods_used = []
            life_stages_used = []

            for cohort_index, cohort in enumerate(cohorts):
                stage_params = self.get_life_stage_params(
                    species,
                    cohort,
                )

                (
                    cohort_delta_msed,
                    cohort_agb_estimate,
                    cohort_method,
                ) = self.compute_cohort_om(
                    cohort=cohort,
                    stage_params=stage_params,
                    simstate=simstate,
                )

                species_delta_msed += cohort_delta_msed

                if cohort_agb_estimate is not None:
                    if species_agb_estimate is None:
                        species_agb_estimate = (
                            cohort_agb_estimate.copy()
                        )
                    else:
                        species_agb_estimate += (
                            cohort_agb_estimate
                        )

                methods_used.append(cohort_method)
                life_stages_used.append(
                    int(cohort.lifestage)
                )

                self.print_life_stage_verification(
                    species=species,
                    cohort_index=cohort_index,
                    cohort=cohort,
                    stage_params=stage_params,
                    method=cohort_method,
                    agb_estimate=cohort_agb_estimate,
                    delta_msed=cohort_delta_msed,
                    simstate=simstate,
                )

            storage_indices = self.get_storage_indices(
                storage=storage,
                delta_msed=species_delta_msed,
                storage_name=storage_name,
                n_active_cells=n_active_cells,
            )

            if storage_name == "msed":
                species_before = np.nansum(
                    storage[:, :, self.target_sediment]
                )

                storage[
                    storage_indices,
                    self.target_layer,
                    self.target_sediment,
                ] += species_delta_msed

                species_after = np.nansum(
                    storage[:, :, self.target_sediment]
                )

            elif storage_name == "bodsed":
                species_before = np.nansum(
                    storage[:, self.target_sediment]
                )

                storage[
                    storage_indices,
                    self.target_sediment,
                ] += species_delta_msed

                species_after = np.nansum(
                    storage[:, self.target_sediment]
                )

            species_added_mass = (
                species_after
                - species_before
            )

            total_added_mass += species_added_mass
            self.cumulative_added_mass += species_added_mass

            unique_methods = sorted(set(methods_used))

            if len(unique_methods) == 1:
                species_method = unique_methods[0]
            else:
                species_method = "mixed"

            if species_agb_estimate is None:
                agb_max = np.nan
                agb_mean = np.nan

                agb_for_output = np.full_like(
                    fveg,
                    np.nan,
                    dtype=float,
                )

            else:
                agb_max = np.nanmax(
                    species_agb_estimate
                )

                agb_mean = np.nanmean(
                    species_agb_estimate
                )

                agb_for_output = (
                    species_agb_estimate.copy()
                )

            organic_records.append(
                {
                    "species": species.name,
                    "method": species_method,

                    # Life-stage diagnostics
                    "life_stages": ",".join(
                        str(stage)
                        for stage in sorted(set(life_stages_used))
                    ),
                    "cohort_count": len(cohorts),

                    "storage": storage_name,
                    "target_layer": (
                        self.target_layer
                        if storage_name == "msed"
                        else ""
                    ),
                    "target_sediment": self.target_sediment,
                    "fveg_max": np.nanmax(fveg),
                    "fveg_mean": np.nanmean(fveg),
                    "agb_max": agb_max,
                    "agb_mean": agb_mean,
                    "organic_increment_max": np.nanmax(
                        species_delta_msed
                    ),
                    "organic_increment_mean": np.nanmean(
                        species_delta_msed
                    ),
                    "summed_areal_increment_kg_m2": species_added_mass,
                    "cumulative_summed_areal_increment_kg_m2": self.cumulative_added_mass,
                }
            )

            organic_cell_records.append(
                {
                    "species": species.name,
                    "fveg": fveg.copy(),
                    "agb_estimate": agb_for_output,
                    "delta_msed": (
                        species_delta_msed.copy()
                    ),
                }
            )

            if self.verbose:
                print(
                    "OrganicAccretion: "
                    f"species={species.name}, "
                    f"storage={storage_name}, "
                    f"methods={unique_methods}, "
                    f"cohorts={len(cohorts)}, "
                    f"life_stages={life_stages_used}, "
                    f"fveg max={np.nanmax(fveg)}, "
                    f"fveg mean={np.nanmean(fveg)}, "
                    f"organic increment max="
                    f"{np.nanmax(species_delta_msed)}, "
                    f"species summed areal increment={species_added_mass}, "
                    f"agb_estimate max={agb_max}, "
                    f"agb_estimate mean={agb_mean}"
                )

        engine.dflowfm.set_var(
            storage_name,
            storage,
        )

        if storage_name == "msed":
            after_total = np.nansum(
                storage[:, :, self.target_sediment]
            )

        elif storage_name == "bodsed":
            after_total = np.nansum(
                storage[:, self.target_sediment]
            )

        added_mass = after_total - before_total

        if self.verbose:
            print(
                "OrganicAccretion: "
                f"storage={storage_name}, "
                f"species_count={len(species_units)}, "
                f"total summed areal increment={added_mass}"
            )

        return organic_records, organic_cell_records
