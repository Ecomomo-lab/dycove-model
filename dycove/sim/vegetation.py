###############################################################
#  vegetation.py
###############################################################

from pathlib import Path
import numpy as np

from dycove.utils.log import Reporter
from dycove.utils.array_math import cell_averaging, sum_product, sum_elementwise
from dycove.sim.vegetation_data import VegetationAttributes, VegCohort

r = Reporter()


##### ---------- Shared methods ---------- #####
class SharedVegMethods:

    def compute_veg_model_quantities(self):
        """
        Function to compute weighted area averages of model vegetation variables.
        Handles the various fractions that are tracked within each grid cell.
        This function works the same no matter the number of species, and so it is shared by both classes
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
    
    def __init__(self, 
                 input_veg_filename,
                 species_name="veg1",
                 mor=0,
                 rand_seed_frac=1.0,
                 rand_seed_method="deterministic",
                 ):
        """
        input_veg_filename: text file containing vegetation attributes
        species_name: string representing name of species (can be anything)
        mor: 1 if we want to consider morphology (burial/scour) for this vegetation unit, else 0
        rand_seed_frac: fraction of potentially colonized cells that will actually colonize, distributed randomly
        rand_seed_method: 'random' or 'deterministic'. If rand_seed_frac is between 0 and 1, 
                          'deterministic' will use a seed to evaluate the SAME set of cells for each colonization, 
                          whereas 'random' will be truly random each time
        """
        self.attrs       = self.parse_input_file(Path(input_veg_filename))
        self.name        = species_name
        self.mor         = mor
        self.seed_frac   = rand_seed_frac
        self.seed_method = rand_seed_method

        # will be come list of VegCohort objects, one for each cohort/colonization that occurs
        self.cohorts: list[VegCohort] = []


    # parse vegetation input file and store attributes in a dataclass
    @staticmethod
    def parse_input_file(filename: Path) -> VegetationAttributes:
        if not filename.exists():
            msg = f"Vegetation input file not found: {filename}"
            r.report(msg, level="ERROR")
            raise FileNotFoundError(msg)
        with open(filename, "r") as f:
            lines = f.readlines()[41:]

        # get general vegetation attributes, convert to ints and floats
        gen = [float(v) if "." in v else int(v) for v in lines[0].split()]
        # get ecological timesteps where colonisation starts and stops
        start_col, end_col = map(int, lines[1].split())
        # get vegetation life stage attributes, convert to ints and floats
        nls = int(gen[5])
        ls = [[float(v) if "." in v else int(v) for v in lines[n+2].split()] for n in range(nls)]

        return VegetationAttributes(  # wow thats a lot of attributes!
            age_max=gen[0],
            seed_dispersal=gen[1],
            col_method=gen[2],
            growth_method=gen[3],
            veg_formula=gen[4],
            nls=nls,
            rootlength_0=gen[6],
            shootlength_0=gen[7],
            stemdiam_0=gen[8],
            veg_type=gen[9],
            start_growth_ets=gen[13],
            end_growth_ets=gen[14],
            winter_ets=gen[15],
            start_col_ets=start_col,
            end_col_ets=end_col,
            ht_max=[row[0] for row in ls],
            diam_max=[row[1] for row in ls],
            rootlength_max=[row[2] for row in ls],
            years_max=[row[3] for row in ls],
            stemdens=[row[4] for row in ls],
            fraction_0=ls[0][5],
            drag=[row[6] for row in ls],
            mud_percent=[row[7] for row in ls],
            desic_no_mort=[row[8] for row in ls],
            desic_all_mort=[row[9] for row in ls],
            flood_no_mort=[row[10] for row in ls],
            flood_all_mort=[row[11] for row in ls],
            uproot_no_mort=[row[12] for row in ls],
            uproot_all_mort=[row[13] for row in ls],
            ht_winter_max=[row[14] for row in ls],
        )


    def colonization(self, ets, min_depths, max_depths, fl_dr, combined_cohorts=None):
        """
        This function only adds new fractions to cells, effective stem diams, heights are computed elsewhere

        min_depths      : array of minimum depths at each cell over the previous period
        max_depths      : array of maximum depths at each cell over the previous period
        fl_dr           : scalar threshold [meters] above which a cell is considered wet
        combined_cohorts: relevant only if multiple species, 
                          each species needs info on other species taking up room in cells
        """
        if self.attrs.start_col_ets <= ets < self.attrs.end_col_ets:
            r.report(f"Doing colonization for species \"{self.name}\", at eco step {ets}")
            # get conditions for colonization
            dry_cond = (min_depths <  fl_dr)
            fld_cond = (max_depths >= fl_dr)
            # get indices for potential colonization
            potential_inds = np.where(dry_cond & fld_cond) #[0]
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
                height=self.attrs.shootlength_0,
                rootlength=self.attrs.rootlength_0,
                lifestage=1,  # life stage always starts at 1 with colonization. TODO: check this
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
        """ Delegate to internal methods """
        self.mortality_hydrodynamic(**hydro_vars)
        self.mortality_morphodynamic(**morpho_vars)
        self.apply_mortality()

    def mortality_hydrodynamic(self, fld_frac, dry_frac, vel_max):
        """
        Compute linear mortality function for hydrodynamic stressors.
        Hydrodynamic mortality is computed independently from vegetation characteristics like stem height, root length, etc.
        Mortality parameters can be dependent on life stage though, so we compute mortality for each life stage and apply 
        to the fractions later (self.apply_mortality).
        Flooding/dessication generally do not depend on life stage, but we currently have the option to have them vary, so 
        we treat them here as if they do. TODO: rearrange input variables to remove this choice?
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
        bl_diff    : difference in bed level [m], from beginning to end of ets, where positive values signify scour
        burial_frac: fraction of shoot length used as threshold for burial
        scour_frac : fraction of root length used as threshold for scour
        """
        # loop over fractions and their current shoot/root lengths and compare to erosion/sedimentation
        # morphodynamic mortality modes are not functions of life stage, and they are binary
        for c in self.cohorts:
            # if morphology is off, we still need these zero arrays for calculations in apply_mortality()
            if self.mor == 0:
                c.potential_mort_burial = np.zeros_like(c.fraction)
                c.potential_mort_scour  = np.zeros_like(c.fraction)
            else:
                c.potential_mort_burial = np.where(bl_diff < -burial_frac*c.height, 1, 0)
                c.potential_mort_scour  = np.where(bl_diff > scour_frac*c.rootlength, 1, 0)


    def apply_mortality(self):
        # multiply potential mortality by vegetation fractions to determine actual mortality, which we track
        # TODO: we track the fractions lost via each cause, but if all causes yield mortality fractions of 1,
        #       how should we track causes when, for instance, a cell only had a single fraction covering 40%?
        #       Is there an order of operations?
        # append empty lists (of length n_fractions) to running list of mortality fractions
        # empty slots will be filled with arrays of mortality fractions corresponding to each fraction that 
        #   currently exists

        # loop over cohorts, add mortality values/causes as attribute of the cohort (both potential and applied)
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
            # TODO: verify that this works the same regardless of the number of species
            c.applied_mort_total = c.applied_mort_flood + c.applied_mort_desic + c.applied_mort_uproot + \
                                   c.applied_mort_burial + c.applied_mort_scour
            fractions_left = c.fraction - c.applied_mort_total
            fractions_left = np.maximum(fractions_left, 0)  # no negative fractions

            # update fractions in cohort
            # for fractions that decay slowly over time, round down to zero when they get small enough
            c.fraction = np.where(fractions_left > 0.025, fractions_left, 0.)


    @staticmethod
    def linear_mortality_func(stressor, th_min, th_max):
        """
        Hydrodynamic stressor mortality function (linear)

        stressor: for flooding/dessication, it is 1-D array of fraction of time within ets where cells are wet/dry
                  for uprooting, it is 1-D array of maximum velocities
        th_min: stressor value below which there is no mortality
        th_max: stressor value above which there is total mortality
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
        """
        Function to compute vegetation height
        (linear growth, constant during summer, constant during winter at winter value)
        """
        # during first growth ets, height starts at previous winter value, so no need to loop
        if ets > self.attrs.start_growth_ets:
            for c in self.cohorts:
                # if none of these conditions are met, height stays the same
                if ets <= self.attrs.end_growth_ets and c.height < self.attrs.ht_max[c.lifestage-1]:
                    c.height += self.attrs.ht_growth_rates[c.lifestage-1]
                elif ets >= self.attrs.winter_ets:
                    # drop down to winter height, but by some chance if we are already below winter max, do nothing
                    c.height = min(c.height, self.attrs.ht_winter_max[c.lifestage-1])


    def stemdiam_growth(self, ets):
        """
        Function to compute vegetation stem diameter growth.
        Linear growth until winter, then constant, then continued linear growth the following year.
        Process continues until max diameter is reached.
        """
        # during first growth ets, height starts at previous winter value, so no need to loop
        if ets > self.attrs.start_growth_ets:
            for c in self.cohorts:
                if c.diameter < self.attrs.diam_max[c.lifestage-1] and ets < self.attrs.winter_ets:
                    c.diameter += self.attrs.diam_growth_rates[c.lifestage-1]  # else, remains constant


    def root_growth(self, ets):  
        """
        Function to compute vegetation root growth.
        Linear growth until winter, then constant, then continued linear growth the following year.
        Process continues until max root length is reached.
        """
        # during first growth ets, height starts at previous winter value, so no need to loop
        if ets > self.attrs.start_growth_ets:
            for c in self.cohorts:
                if c.rootlength < self.attrs.rootlength_max[c.lifestage-1] and ets < self.attrs.winter_ets:
                    c.rootlength += self.attrs.root_growth_rates[c.lifestage-1]  # else, remain constant


    def update_lifestage_and_stemdensity(self):
        """ Function to be called at the end of every eco year """
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
    - A class for handling multiple vegetation species.
    - We want to keep the hydrodynamics code clean, avoiding extra if statements regarding the number of species.
    - This class handles the extra steps required, mostly containing wrappers of functions from VegetationSpecies 
      that distribute the tasks to each species present.
    - The functions below are in the order they are called in the hydrodynamics.run_simulation module
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

    def apply_mortality(self):
        for sp in self.species_list:
            sp.apply_mortality()

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











    #     # inherit the "mor" input parameter from the individual species object
    #     # we cannot have some species with mor=1 and others with mor=0
    #     if len(set([veg.mor for veg in self.species_list])) != 1:
    #         msg = "All vegetation species inputs must have the same value for 'mor'"
    #         logger.error(msg)
    #         raise ValueError(msg)
    #     self.mor = self.species_list[0].mor

    #     # a dictionary for tracking fractions across multiple species. As more fractions are colonized, a unique 
    #     #   ID (veg_fraction_ID) is appended to the list corresponding to each species so that we can access them 
    #     #   later when tracking mortality. See colonization method for more explanation.
    #     self.veg_fraction_inds = {}
    #     self.veg_fraction_ind = 0
    #     for veg in self.species_list:
    #         self.veg_fraction_inds[veg.name] = []

    #     # master list of all veg fractions, identical in shape to veg_fractions lists for individual species,
    #     #   but it includes cell-wise fractions for all species for accounting purposes. For example, for one
    #     #   species, the arrays within the veg_fractions list are arrays of fractions for each consecutive 
    #     #   cohort/colonization. For two species, each species colonizes during the colonization time step,
    #     #   and their fractions are appended to this list. So, every other array within veg_fractions will
    #     #   correspond to the same species. For three species, every third array, and so on.
    #     self.veg_fractions = []
    #     # also keep lists of other vegetation parameters. By keeping these lists and using the same variable 
    #     #   names as for the individual species, we can run shared functions (e.g., compute_veg_model_quantities) 
    #     #   without needing to pass arguments or make adjustments
    #     self.veg_densities = []
    #     self.veg_diameters = []
    #     self.veg_heights   = []
    #     self.mort_from_flood  = []
    #     self.mort_from_desic  = []
    #     self.mort_from_uproot = []
    #     self.mort_from_burial = []
    #     self.mort_from_scour  = []
    #     self.mort_total       = []

    # def update_lifestage_and_stemdensity(self):
    #     for veg in self.species_list:
    #         # update lifestages and associated stem densities for individual species and update the corresponding lists for that species
    #         veg.update_lifestage_and_stemdensity()
    #         # iterate through once more, taking the latest stem density values and updating the global list
    #         for i, vdens in enumerate(veg.veg_densities):
    #             corresponding_ind = self.veg_fraction_inds[veg.name][i]
    #             self.veg_densities[corresponding_ind] = vdens


    # def stemheight_growth(self, ets):
    #     for veg in self.species_list:
    #         # compute stem height growth for individual species and update the corresponding list for that species
    #         veg.stemheight_growth(ets)
    #         # iterate through once more, taking the latest stem height values and updating the global list
    #         for i, vht in enumerate(veg.veg_heights):
    #             corresponding_ind = self.veg_fraction_inds[veg.name][i]
    #             self.veg_heights[corresponding_ind] = vht

    # def stemdiam_growth(self, ets):
    #     for veg in self.species_list:
    #         # compute stem diameter growth for individual species and update the corresponding list for that species
    #         veg.stemdiam_growth(ets)
    #         # iterate through once more, taking the latest stem diameter values and updating the global list
    #         for i, vd in enumerate(veg.veg_diameters):
    #             corresponding_ind = self.veg_fraction_inds[veg.name][i]
    #             self.veg_diameters[corresponding_ind] = vd            

    # def root_growth(self, ets):
    #     # we do not track root_growth globally because we don't export this at the end of the run (as of now)
    #     for veg in self.species_list:
    #         veg.root_growth(ets)

    # def colonization(self, ets, min_depths, max_depths, fl_dr):
    #     for veg in self.species_list:
    #         # for the first colonization of the first species, self.veg_fractions_combined will be an empty list,
    #         #   and colonization will proceed normally
    #         veg.colonization(ets, min_depths, max_depths, fl_dr, combined_fractions_list=self.veg_fractions)
    #         # update the global fractions list
    #         self.update_fractions_with_new_colonization(ets, veg)
    #     logger.info(f"Number of veg fractions total: {len(self.veg_fractions)}")

    # def update_fractions_with_new_colonization(self, ets, veg):
    #     # ETS checks happen here and under colonization() because they depend on inputs for individual species (start/end)
    #     if veg.start_col_ets <= ets < veg.end_col_ets:
    #         # after colonization of an individual species, we append that most recent array of cell-wise fractions
    #         #   and corresponding parameters to the master lists, so when we do colonization for the next species, 
    #         #   we have a global account of cell-wise fractions and therefore available space in each cell. We use
    #         #   veg_fraction_inds to track fractions for each species within that master list, even though the
    #         #   fractions arrays will be spaced regularly within the master list. 
    #         self.veg_fractions.append(veg.veg_fractions[-1])
    #         # this means, for the first colonization of species 1: veg_fraction_inds = {"species 1": [0], "species 2" = []}
    #         # and for the first colonization of species 2:         veg_fraction_inds = {"species 1" = [0], "species 2" = [1]}
    #         # and for the second colonization of species 1:        veg_fraction_inds = {"species 1" = [0, 2], "species 2" = [1]}
    #         #   and so on. Clearly, the arrays in veg_fractions will alternate between species, but the structure of
    #         #   the list is the same as if there were only one species, which makes the cell-averaging calculations
    #         #   in compute_veg_model_quantities simple
    #         self.veg_fraction_inds[veg.name].append(self.veg_fraction_ind)
    #         self.veg_densities.append(veg.veg_densities[-1])
    #         self.veg_diameters.append(veg.veg_diameters[-1])
    #         self.veg_heights.append(veg.veg_heights[-1])
    #         self.veg_fraction_ind += 1

    # def mortality_hydrodynamic(self, fld_fraction, dry_fraction, vel_maximum):
    #     for veg in self.species_list:
    #         # Compute the potential mortality for each species, for each stressor, individually. 
    #         # This will yield mortality from hydrodynamic causes as fractions between 0 and 1 for each grid cell, REGARDLESS of 
    #         #   vegetation present in the cell. But potential mortality is dependent on species properties, so they are computed 
    #         #   separately.
    #         veg.mortality_hydrodynamic(fld_fraction, dry_fraction, vel_maximum)

    # # TODO: find a better way to share these function definitions, so that the identical arguments are not written out in two places
    # def mortality_morphodynamic(self, bl_diff, burial_frac=1.0, scour_frac=0.1):
    #     for veg in self.species_list:
    #         # See comments under mortality_hydrodynamic()
    #         veg.mortality_morphodynamic(bl_diff, burial_frac, scour_frac)

    # def apply_mortality(self):
    #     for veg in self.species_list:
    #         veg.apply_mortality()
    #         # Below, we redo the loop that occurs in apply_mortality() and apply updated fractions to the global list.
    #         # veg.veg_fractions is the individual species' fractions list
    #         # self.veg_fractions is the global fractions list
    #         for i in range(len(veg.veg_fractions)):
    #             # NT verified that this works correctly 7/17/25
    #             # i = 0, 1, 2, ...
    #             # vf_ID should be = 0, 2, 4, .... if this is the first of two species
    #             vf_ind = self.veg_fraction_inds[veg.name][i]
    #             # Global list updated with fractions from individual veg objects that were previously updated in apply_mortality()
    #             self.veg_fractions[vf_ind] = veg.veg_fractions[i]
    #         # Add mortality arrays from individual lists to global lists
    #         # This keeps the mortality arrays in the same order as fractions, so they can be indexed with veg_fractions_inds
    #         for attr in ["mort_from_flood", "mort_from_desic", "mort_from_uproot", 
    #                      "mort_from_burial", "mort_from_scour", "mort_total"]:
    #             getattr(self, attr).append(getattr(veg, attr)[-1])
