###############################################################
#  vegetation.py
###############################################################

from pathlib import Path
import numpy as np
import json

from dycove.sim.vegetation_data import VegetationAttributes, VegCohort
from dycove.utils.simulation_reporting import Reporter
from dycove.utils.array_math import cell_averaging, sum_product, sum_elementwise

r = Reporter()


##### ---------- Shared methods ---------- #####
class SharedVegMethods:

    def compute_veg_model_quantities(self):
        """
        Function to compute weighted area averages of model vegetation variables.

        Handles the various fractions that are tracked within each grid cell.

        This function works the same no matter the number of species, and so it is shared 
        by both classes.
        """

        # for single species, self.cohorts is a list of VegCohort objects
        # for multiple species, self.cohorts is a flattened list of VegCohort objects for all species
        cohorts = self.cohorts
        if not self.cohorts:
            self.stemdensity  = None
            self.stemdiameter = None
            self.stemheight   = None
            return
        
        # get list of cohort fractions and other attributes
        fractions = [c.fraction for c in cohorts]
        densities = [c.density for c in cohorts]
        diameters = [c.diameter for c in cohorts]
        heights   = [c.height   for c in cohorts]        
        # stem density in a single cell: sum of n_i*frac_i over all fractions frac_i in cell
        self.stemdensity  = sum_product(fractions, densities)
        # for veg height/diameter we just need a weighted average value in each cell based on cell fractions.
        # cells may have a sum of fractions less than 1, so we need to safely divide by the sum of fractions
        self.stemdiameter = cell_averaging(fractions, diameters)
        self.stemheight   = cell_averaging(fractions, heights)


class VegetationSpecies(SharedVegMethods):
    """
    A single vegetation species. Can (and will) handle multiple cohorts of the species.

    Parameters
    ----------
    input_veg_filename : str or Path
        Path to JSON file containing vegetation attributes.
    species_name : str
        String representing name of species (can be anything).
    mor : int
        1 if we want to consider morphology (burial/scour) for this vegetation unit, else 0
    rand_seed_frac : float
        Fraction of potentially colonized cells that will actually colonize, distributed randomly.
    rand_seed_method : str
        Options are 'random' or 'deterministic'. If rand_seed_frac is between 0 and 1, 
        'deterministic' will use a seed to evaluate the SAME set of cells for each colonization, 
        whereas 'random' will be truly random each time.
    """

    def __init__(self, 
                 input_veg_filename,
                 species_name="veg1",
                 mor=0,
                 rand_seed_frac=1.0,
                 rand_seed_method="random",
                 ):

        self.attrs       = self.load_vegetation_attributes(Path(input_veg_filename))
        self.name        = species_name
        self.mor         = mor
        self.seed_frac   = rand_seed_frac
        self.seed_method = rand_seed_method

        # will be come list of VegCohort objects, one for each cohort/colonization that occurs
        self.cohorts: list[VegCohort] = []


    @staticmethod
    def load_vegetation_attributes(filename: Path) -> VegetationAttributes:
        """ Parse vegetation input json file and store attributes in a dataclass. """
        with open(filename, "r") as f:
            data = json.load(f)

        ls_data = data.pop("life_stage_attr")
        for key in ls_data[0].keys():
            data[key] = [stage[key] for stage in ls_data]

        return VegetationAttributes(**data)


    def colonization(self, ets, min_depths, max_depths, fl_dr, combined_cohorts=None):
        """
        Colonize a :class:`~dycove.sim.vegetation_data.VegCohort` by adding a new 
        fraction to cells.

        This method only adds new fractions. Effective stem diameters and heights are 
        computed elsewhere.

        Parameters
        ----------
        ets : int
            Current ecological timestep within current ecological year.
        min_depths : numpy.ndarray
            Array of minimum water depths [m] at each cell over the previous period.
        max_depths : numpy.ndarray
            Array of maximum water depths [m] at each cell over the previous period.
        fl_dr : float
            Wet/dry threshold [m]; cells with depth above this value are considered wet.
        combined_cohorts : list of VegetationSpecies or None, optional
            Relevant only if multiple species are present. Provides information about
            other species occupying space in cells. See colonization method of 
            :class:`~dycove.sim.vegetation.MultipleVegetationSpecies` for the use-case.
        """

        if self.attrs.start_col_ets <= ets < self.attrs.end_col_ets:
            r.report(f"Doing colonization for species \"{self.name}\", at eco step {ets}")

            # get conditions for colonization
            dry_cond = (min_depths <  fl_dr)
            fld_cond = (max_depths >= fl_dr)
            # get indices for potential colonization
            potential_inds = np.where(dry_cond & fld_cond)[0]
            # create mask for potential colonization based on input fraction and method
            rand_mask = self.create_seed_fraction_mask(len(min_depths))
            # apply potential inds to mask
            potential_inds_mask = np.zeros_like(rand_mask, dtype=bool)
            potential_inds_mask[potential_inds] = True
            selected_inds = potential_inds_mask & rand_mask

            existing_cohorts = combined_cohorts if combined_cohorts is not None else self.cohorts
            new_fraction = self.compute_new_fraction(existing_cohorts, selected_inds, len(min_depths))

            # create new cohort for latest colonization
            # new_veg_fraction is an array, other quantities are scalars
            new_cohort = VegCohort(
                fraction=new_fraction,
                density=self.attrs.stemdens[0],
                diameter=self.attrs.stemdiam_0,
                height=self.attrs.stemht_0,
                rootlength=self.attrs.rootlength_0,
                lifestage=1,
                lifestage_year=1,
            )
            self.cohorts.append(new_cohort)


    def create_seed_fraction_mask(self, array_len):
        # number of indices where we will allow colonization
        n_inds = int(np.floor(self.seed_frac*array_len))
        # initialize mask with all False
        rand_mask = np.zeros(array_len, dtype=bool)
        if self.seed_method == "deterministic":
            np.random.seed(0)            
        inds = np.random.choice(array_len, n_inds, replace=False)
        rand_mask[inds] = True
        return rand_mask


    def compute_new_fraction(self, existing_cohorts, inds2colonize, array_len):
        # compute a new cohort’s fractional cover based on available space
        new_fraction = np.zeros(array_len)
        existing_fractions = [c.fraction for c in existing_cohorts]

        if not existing_fractions:
            # create vegetation fraction array based on the hydro conditions and random filtering
            new_fraction[inds2colonize] = self.attrs.fraction_0
        else:
            # TODO: Maybe find a better way of avoiding continuous appends with each colonization. 
            #       E.g., avoid adding a new fraction if no space anywhere, etc. But that would mean that ZERO cells 
            #       anywhere have favorable hydro conditions AND space available (not many cases).
            # determine available space based on sum of fractions in each cell
            fraction_uncovered = 1 - sum_elementwise(existing_fractions)
            # compute new fraction to be colonized; min: space available, max: initial veg fraction
            candidate_fraction = np.minimum(fraction_uncovered, self.attrs.fraction_0)
            # we apply new fractions where partial flooding condition is met, and within cells selected via seed_frac
            new_fraction[inds2colonize] = candidate_fraction[inds2colonize]

        return new_fraction


    def mortality(self, hydro_vars, morpho_vars):
        """ Delegate to internal methods. """
        self.mortality_hydrodynamic(**hydro_vars)
        self.mortality_morphodynamic(**morpho_vars)
        #self.apply_mortality()
        self.apply_mortality_using_initial_fractions()

    def mortality_hydrodynamic(self, fld_frac, dry_frac, vel_max):
        """
        Compute linear mortality functions for each hydrodynamic stressor (if activated).

        Hydrodynamic mortality is computed independently from vegetation characteristics
        like stem height, root length, etc. Mortality parameters can be dependent on life 
        stage though, so we compute mortality for each life stage and apply to the 
        fractions later (self.apply_mortality).

        Parameters
        ----------
        fld_frac : numpy.ndarray
            Array containing fractions of time that each cell was flooded during the 
            previous period.
        dry_frac : numpy.ndarray
            Array containing fractions of time that each cell was dry during the 
            previous period.
        vel_max : numpy.ndarray
            Array of maximum water velocities [m/s] at each cell over the previous period.
        """

        # these arrays of potential mortality get saved as cohort attributes, along with applied mortality below
        self.potential_mort_flood  = [None]*self.attrs.nls
        self.potential_mort_desic  = [None]*self.attrs.nls
        self.potential_mort_uproot = [None]*self.attrs.nls
        # loop over life stages present in vegetation file, thresholds depend on lifestage
        for n in range(self.attrs.nls):
            if self.attrs.flood_no_mort[n] == 0:  # if flood mortality is turned off, flag = 0 (TODO: create a better flag)
                self.potential_mort_flood[n] = np.zeros_like(fld_frac)
            else:
                self.potential_mort_flood[n] = self.linear_mortality_func(fld_frac, 
                                                                          self.attrs.flood_no_mort[n], 
                                                                          self.attrs.flood_all_mort[n])
            if self.attrs.desic_no_mort[n] == 0:  # if dessication is turned off, flag = 0 (TODO: create a better flag)
                self.potential_mort_desic[n] = np.zeros_like(dry_frac)
            else:
                self.potential_mort_desic[n] = self.linear_mortality_func(dry_frac, 
                                                                          self.attrs.desic_no_mort[n], 
                                                                          self.attrs.desic_all_mort[n])
            if self.attrs.uproot_no_mort[n] == 0:  # if uprooting is turned off, flag = 0 (TODO: create a better flag)
                self.potential_mort_uproot[n] = np.zeros_like(vel_max)
            else:
                self.potential_mort_uproot[n] = self.linear_mortality_func(vel_max,
                                                                           self.attrs.uproot_no_mort[n],
                                                                           self.attrs.uproot_all_mort[n])
                

    # TODO: verify that we want the default scour_frac to be 10%, previous codes have just used 100% same as stem burial
    def mortality_morphodynamic(self, bl_diff=None, burial_frac=1.0, scour_frac=0.1):
        """
        Compute linear mortality functions for each morphodynamic stressor (if activated).

        Morphodynamic stressors (burial, scour) are not currently functions of life stage,
        and they are binary (e.g., stem is either buried or not).

        Creates arrays of zeros if using a non-morphology model or with morphology turned 
        off.

        Parameters
        ----------
        bl_diff : numpy.ndarray
            Array of cell-wise differences in bed level [m], from beginning to end of the 
            previous period, where positive values signify burial.
        burial_frac : float
            Fraction of stem height above which vegetation is considered buried.
        scour_frac : float
            Fraction of root length above which vegetation is considered scoured.

        Notes
        -----
        - If morphology is turned on (mor=1), but morphology is not active in the model (e.g., ANUGA)
          bl_diff will be passed as zero-difference arrays, and scour/burial will be set to zero.
        """

        # loop over fractions and their current shoot/root lengths and compare to erosion/sedimentation
        for c in self.cohorts:
            # if morphology is off, we still need these zero arrays for calculations in apply_mortality()
            if self.mor == 0:
                c.potential_mort_burial = np.zeros_like(c.fraction)
                c.potential_mort_scour  = np.zeros_like(c.fraction)
            else:
                c.potential_mort_burial = np.where(bl_diff >= burial_frac*c.height, 1, 0)
                c.potential_mort_scour  = np.where(bl_diff <= -scour_frac*c.rootlength, 1, 0)


    def apply_mortality(self):
        """
        Multiply potential mortality by vegetation fractions to determine actual mortality.

        Populate mortality-related fields of each active VegCohort object.
        """

        # TODO: we track the fractions lost via each cause, but if all causes yield mortality fractions of 1,
        #       how should we track causes when, for instance, a cell only had a single fraction covering 40%?
        #       Is there an order of operations?
        for c in self.cohorts:
            # vegetation fractions lost to flooding
            c.potential_mort_flood = self.potential_mort_flood[c.lifestage-1]
            c.applied_mort_flood  = c.fraction*c.potential_mort_flood
            # vegetation fractions lost to dessication
            c.potential_mort_desic = self.potential_mort_desic[c.lifestage-1]
            c.applied_mort_desic  = c.fraction*c.potential_mort_desic 
            # vegetation fractions lost to uprooting
            c.potential_mort_uproot = self.potential_mort_uproot[c.lifestage-1]
            c.applied_mort_uproot = c.fraction*c.potential_mort_uproot
            # vegetation fractions lost to deposition
            c.applied_mort_burial = c.fraction*c.potential_mort_burial
            # vegetation fractions lost to erosion
            c.applied_mort_scour  = c.fraction*c.potential_mort_scour
            # subtract all mortality fractions from actual fractions, but maintain minimum fraction of zero
            c.applied_mort_total = c.applied_mort_flood + c.applied_mort_desic + c.applied_mort_uproot + \
                                   c.applied_mort_burial + c.applied_mort_scour
            fractions_left = c.fraction - c.applied_mort_total
            fractions_left = np.maximum(fractions_left, 0)  # no negative fractions

            # update fractions in cohort
            # for fractions that decay slowly over time, round down to zero when they get small enough
            c.fraction = np.where(fractions_left > 0.025, fractions_left, 0.)


    def apply_mortality_using_initial_fractions(self):
        """ Replacing `c.fraction` with `self.attrs.fraction_0`, always a function of initial colonization """
        for c in self.cohorts:
            # vegetation fractions lost to flooding
            c.potential_mort_flood = self.potential_mort_flood[c.lifestage-1]
            # for accounting purposes, no mortality if there is no fraction
            c.applied_mort_flood  = np.where(c.fraction > 0.01, self.attrs.fraction_0*c.potential_mort_flood, 0.)
            # vegetation fractions lost to dessication
            c.potential_mort_desic = self.potential_mort_desic[c.lifestage-1]
            c.applied_mort_desic  = np.where(c.fraction > 0.01, self.attrs.fraction_0*c.potential_mort_desic, 0.)
            # vegetation fractions lost to uprooting
            c.potential_mort_uproot = self.potential_mort_uproot[c.lifestage-1]
            c.applied_mort_uproot = np.where(c.fraction > 0.01, self.attrs.fraction_0*c.potential_mort_uproot, 0.)
            # vegetation fractions lost to deposition
            c.applied_mort_burial = np.where(c.fraction > 0.01, self.attrs.fraction_0*c.potential_mort_burial, 0.)
            # vegetation fractions lost to erosion
            c.applied_mort_scour  = np.where(c.fraction > 0.01, self.attrs.fraction_0*c.potential_mort_scour, 0.)
            # subtract all mortality fractions from actual fractions, but maintain minimum fraction of zero
            c.applied_mort_total = c.applied_mort_flood + c.applied_mort_desic + c.applied_mort_uproot + \
                                    c.applied_mort_burial + c.applied_mort_scour
            fractions_left = c.fraction - c.applied_mort_total
            c.fraction = np.maximum(fractions_left, 0)  # no negative fractions


    @staticmethod
    def linear_mortality_func(stressor, th_min, th_max):
        """
        Generic, hydrodynamic stressor mortality function (linear).

        Parameters
        ----------
        stressor : numpy.ndarray
            Array of relevant stressor magnitude for each cell. For flooding/dessication,
            it is an array of fraction of time where cells are wet/dry. For uprooting, 
            it is an array of maximum velocities.
        th_min : float
            Stressor value below which there is no mortality. Comes from JSON input file.
        th_max : float
            Stressor value above which there is total mortality. Comes from JSON input file.
        """

        # compute fractional mortality based on linear interpolation
        mort_frac = (stressor - th_min)/(th_max - th_min)
        
        # updated, cleaner method to just bring all fractions outside the 0-1 range to 0 and 1
        mort_frac[mort_frac > 1] = 1
        mort_frac[mort_frac < 0] = 0

        return mort_frac


    # TODO: re-implement this, removed once the multi-species capability was added
    # the likelihood of this being called is fairly low
    def clean_out_old_fractions(self):
        # if all fractions within a given cohort are zero, remove it from the list
        cohorts_to_remove = []
        for i, c in enumerate(self.cohorts):
            if np.sum(c.fraction) < 0.001:
                cohorts_to_remove.append(i)
                r.report("Removed a cohort that has disappeared.")
        if cohorts_to_remove:
            self.remove_old_cohorts(cohorts_to_remove)


    def remove_old_cohorts(self, list_of_inds):
        self.cohorts = [x for i, x in enumerate(self.cohorts) if i not in list_of_inds]

      
    def stemheight_growth(self, ets):
        """ Computes stem height based on growth functions in VegAttributes. """
        # during first growth ets, height starts at previous winter value, so no need to loop
        if ets > self.attrs.start_growth_ets:
            for c in self.cohorts:
                # if none of these conditions are met, height stays the same
                if ets <= self.attrs.end_growth_ets and c.height < self.attrs.stemht_max[c.lifestage-1]:
                    c.height += self.attrs.ht_growth_rates[c.lifestage-1]
                elif ets >= self.attrs.winter_ets:
                    # drop down to winter height, but by some chance if we are already below winter max, do nothing
                    c.height = min(c.height, self.attrs.stemht_winter_max[c.lifestage-1])


    def stemdiam_growth(self, ets):
        """ Computes stem diameter based on growth functions in VegAttributes. """
        # during first growth ets, height starts at previous winter value, so no need to loop
        if ets > self.attrs.start_growth_ets:
            for c in self.cohorts:
                if c.diameter < self.attrs.stemdiam_max[c.lifestage-1] and ets < self.attrs.winter_ets:
                    c.diameter += self.attrs.diam_growth_rates[c.lifestage-1]  # else, remains constant


    def root_growth(self, ets):  
        """ Computes root length based on growth functions in VegAttributes. """
        # during first growth ets, height starts at previous winter value, so no need to loop
        if ets > self.attrs.start_growth_ets:
            for c in self.cohorts:
                if c.rootlength < self.attrs.rootlength_max[c.lifestage-1] and ets < self.attrs.winter_ets:
                    c.rootlength += self.attrs.root_growth_rates[c.lifestage-1]  # else, remain constant


    def update_lifestage_and_stemdensity(self):
        """ 
        Function to be called at the end of every eco year. 

        Updates the life stage and the year within the life stage for each cohort.

        Stem density update is done here because it only depends on eco year and does not 
        change during the year.
        """

        cohorts_to_remove = []
        for i, c in enumerate(self.cohorts):
            # if we have reached final year in the life stage
            if c.lifestage_year == self.attrs.years_max[c.lifestage-1]:
                # if unit has not yet reached final life stage
                if c.lifestage != self.attrs.nls:
                    c.lifestage   += 1  # move this fraction to next life stage
                    c.lifestage_year = 1  # restart counter for years within life stage
                    # TODO: consider moving this elsewhere?
                    c.density = self.attrs.stemdens[c.lifestage-1]  # apply stemdens from the NEXT lifestage (lifestage was just increased by 1)
                else:
                    # REMOVE FRACTION FROM LIST OF FRACTIONS, need to be careful with this
                    cohorts_to_remove.append(i)
            else:
                c.lifestage_year += 1
        if cohorts_to_remove:
            self.remove_old_cohorts(cohorts_to_remove)


class MultipleVegetationSpecies(SharedVegMethods):
    """
    A class for handling multiple vegetation species.

    We want to keep the hydrodynamics code clean, avoiding extra if statements regarding 
    the number of species. This class handles the extra steps required, mostly containing 
    wrappers of functions from VegetationSpecies that distribute the tasks to each species
    present.

    Parameters
    ----------
    species_list : list[VegetationSpecies]
        A list of VegetationSpecies objects.

    """

    def __init__(self, species_list: list[VegetationSpecies]):

        # list of individual vegetation species objects
        self.species_list = species_list
        self.check_species_consistency()

    @property
    def cohorts(self) -> list[VegCohort]:
        # flatten vegetation cohorts across species into single list
        return [c for sp in self.species_list for c in sp.cohorts]

    def colonization(self, ets, min_depths, max_depths, fl_dr):
        for sp in self.species_list:
            # passing flattened list of all cohorts for all species for colonization calculation
            sp.colonization(ets, min_depths, max_depths, fl_dr, combined_cohorts=self.cohorts)
        r.report(f"Number of veg fractions total: {len(self.cohorts)}")

    def stemheight_growth(self, ets):
        for sp in self.species_list:
            sp.stemheight_growth(ets)

    def stemdiam_growth(self, ets):
        for sp in self.species_list:
            sp.stemdiam_growth(ets)

    def root_growth(self, ets):
        for sp in self.species_list:
            sp.root_growth(ets)

    def mortality(self, hydro_vars, morpho_vars):
        for sp in self.species_list:
            sp.mortality(hydro_vars, morpho_vars)

    ### Not called outside of vegetation module, so doesn't need a wrapper method
    # def apply_mortality(self):
    #     for sp in self.species_list:
    #         sp.apply_mortality()

    def update_lifestage_and_stemdensity(self):
        for sp in self.species_list:
            sp.update_lifestage_and_stemdensity()

    def check_species_consistency(self):
        # inherit the "mor" input parameter from the individual species object
        # we cannot have some species with mor=1 and others with mor=0
        if len(set([sp.mor for sp in self.species_list])) != 1:
            msg = "All vegetation species inputs must have the same value for 'mor'"
            r.report(msg, level="ERROR")
            raise ValueError(msg)
        self.mor = self.species_list[0].mor
