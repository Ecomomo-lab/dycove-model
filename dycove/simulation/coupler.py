""" For handling coupling between numerical model and vegetation model """
class VegetationCoupler:
    def __init__(self, engine):
        self.engine = engine
        self.veg = engine.veg

    def update(self, simstate, hydrostats):
        """Advance vegetation by one eco timestep and inject into hydro model."""
        self.lifestage_update(simstate)
        self.apply_growth(simstate)
        self.do_colonization(simstate, hydrostats)
        self.compute_mortality(simstate, hydrostats)
        self.compute_veg_model_quantities()
        self.push_veg_to_hydro()

    def lifestage_update(self, simstate):
        simstate.update_ets()
        if simstate.ets == 1:
            self.veg.update_lifestage_and_stemdensity()

    def apply_growth(self, simstate):
        self.veg.stemheight_growth(simstate.ets)
        self.veg.stemdiam_growth(simstate.ets)
        self.veg.root_growth(simstate.ets)

    def do_colonization(self, simstate, hydrostats):
        self.veg.colonization(simstate.ets, hydrostats.h_min, hydrostats.h_max, hydrostats.fl_dr)

    def compute_mortality(self, simstate, hydrostats):
        substeps = simstate.n_hydro_steps
        hydro_vars  = {"fld_frac"   : hydrostats.flood_frac(substeps),
                       "dry_frac"   : hydrostats.dry_frac(substeps),
                       "vel_max"    : hydrostats.v_max}
        morpho_vars = {"bl_diff"    : hydrostats.bedlevel_diff,
                       "burial_frac": 1.0,
                       "scour_frac" : 0.1}
        self.veg.mortality(hydro_vars, morpho_vars)

    def compute_veg_model_quantities(self):
        self.veg.compute_veg_model_quantities()

    def push_veg_to_hydro(self):
        # if no vegetation in model, skip
        # could be due to delayed colonization, total mortality, or no vegetation on this specific subdomain
        # veg.stemdensity is given a value in veg.compute_veg_model_quantities()
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
