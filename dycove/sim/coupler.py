###############################################################
#  coupler.py
###############################################################

class VegetationCoupler:
    """
    Class to handle coupling between hydrodynamic model and vegetation model.

    This class orchestrates the vegetation update cycle by coordinating between
    the hydrodynamic engine and vegetation model at each ecological timestep.

    Vegetation object (``self.veg``) may be an instance of 
    :class:`~dycove.sim.vegetation.VegetationSpecies` or 
    :class:`~dycove.sim.vegetation.MultipleVegetationSpecies`. If it is the latter, 
    the methods in :meth:`~dycove.sim.coupler.VegetationCoupler.update` are 
    distrubuted to each individual species under the methods of 
    :class:`~dycove.sim.vegetation.MultipleVegetationSpecies`.

    """

    def __init__(self, engine):
        self.engine = engine
        self.veg = engine.veg

    def update(self, simstate, hydrostats):
        """
        Advance vegetation by one eco timestep.

        These are all the steps that need to be performed within one ecological 
        time step.
        """
        self.lifestage_update(simstate)
        self.apply_growth(simstate)
        self.do_colonization(simstate, hydrostats)
        self.compute_mortality(simstate, hydrostats)
        self.compute_veg_model_quantities()
        self.push_veg_to_hydro()

    def lifestage_update(self, simstate):
        """
        Update ETS and eco_year counter, then update vegetation life stage and stem density 
        if we started a new eco year.
        """
        simstate.update_ets()
        if simstate.ets == 1:
            self.veg.update_lifestage_and_stemdensity()

    def apply_growth(self, simstate):
        """ Apply precomputed growth rates for stem height, diameter, and root depth. """
        self.veg.stemheight_growth(simstate.ets)
        self.veg.stemdiam_growth(simstate.ets)
        self.veg.root_growth(simstate.ets)

    def do_colonization(self, simstate, hydrostats):
        """ Compute colonization based on hydrodynamic conditions. """
        self.veg.colonization(simstate.ets, hydrostats.h_min, hydrostats.h_max, hydrostats.fl_dr)

    def compute_mortality(self, simstate, hydrostats):
        """ Compute mortality based on hydrodynamic statistics. """
        n_substeps = simstate.n_hydro_steps
        hydro_vars  = {"fld_frac"   : hydrostats.flood_frac(n_substeps),
                       "dry_frac"   : hydrostats.dry_frac(n_substeps),
                       "vel_max"    : hydrostats.v_max*0.95}
        morpho_vars = {"bl_diff"    : hydrostats.bedlevel_diff,
                       "burial_frac": 1.0,
                       "scour_frac" : 0.1}
        self.veg.mortality(hydro_vars, morpho_vars)

    def compute_veg_model_quantities(self):
        """ 
        Compute vegetation quantities (stem density, diameter, height), to be pushed to 
        hydro model as single values at each grid cell
        """
        self.veg.compute_veg_model_quantities()

    def push_veg_to_hydro(self):
        """
        Inject updated vegetation properties into hydrodynamic model arrays.
        
        Replaces the current vegetation arrays in the hydro model with values from the 
        vegetation module. Skips injection if no vegetation exists (due to delayed 
        colonization, total mortality, or absence on this subdomain).
        """
        # veg.stemdensity is given a value in compute_veg_model_quantities()
        if self.veg.stemdensity is None:
            return

        # get vegetation variable arrays directly from hydro model
        stemdensity, stemdiameter, stemheight = self.engine.get_vegetation()
        
        # replace the "active" part of the vegetation arrays with values from vegetation module (required for DFM, not ANUGA)
        stemdensity[:len(self.veg.stemdensity)] = self.veg.stemdensity
        stemdiameter[:len(self.veg.stemdiameter)] = self.veg.stemdiameter
        stemheight[:len(self.veg.stemheight)] = self.veg.stemheight

        # replace the vegetation attributes in the hydro model
        self.engine.set_vegetation(stemdensity, stemdiameter, stemheight)
