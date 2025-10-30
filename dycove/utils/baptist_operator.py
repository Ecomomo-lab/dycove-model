"""
Vegetation operators - downloaded on 04/22/25
"""
from __future__ import division, absolute_import, print_function
import numpy as np
from numpy import sqrt, minimum, maximum, log

# from anuga import Domain
from anuga import Quantity
from anuga.operators.base_operator import Operator

# from math import sqrt, log
# from anuga.config import epsilon, g

# import anuga.utilities.log as log
# from anuga.config import netcdf_mode_r, netcdf_mode_w, netcdf_mode_a, \
                            # netcdf_float

# import os
# from scipy.interpolate import NearestNDInterpolator

#===============================================================================
# Vegetation operator applying drag to the flow
#===============================================================================
class Baptist_operator(Operator):
    """
    Baptist operator that applies a drag on the flow due to the presence of veg.
    Methodology is based on Baptist et al, 2007 for emergent and submerged veg.
    Modified Chezy coefficient is computed as:
    Cv = (Cb^-2 + (Cd*m*D/(2g))*min(h,hv))^-0.5 + (g^0.5/k)*ln(max(h,hv)/hv)
    
    (Cv is vegetated Chezy, Cb is bed Chezy, Cd is drag coefficient, m is stem
    density [m^-2], D is stem diameter [m], g is gravity, h is flow depth [m], 
    hv is stem height [m], k is von karman constant)
    
    However, we precompute some of these terms together in the init for speed.
    The result looks like:
    Cv = (a1 + a2*CD*min(h,hv))^-0.5 + a3*ln(max(h,hv)/hv)
    
    Operator uses explicit form to update velocities, according to:
    d/dt(uh) = -g*uh*sqrt(uh^2 + vh^2)/(Cv^2*h^2)
    """

    def __init__(self, 
                 domain, 
                 veg_diameter=None,
                 veg_density=None,
                 veg_height=None,
                 bed_friction_const=65.0,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):
        """

        Initialize vegetation characteristics

        **Inputs** :

            domain : `object`
                ANUGA domain instance

            veg_diameter : `float or np.ndarray`
                Vegetation stem diameters, given in meters. Input can be either
                one constant value to be applied everywhere, or an array giving
                one value per cell centroid. If None is given, we check to see
                if values were previously defined in domain.quantities()

            veg_density : `float or np.ndarray`
                Vegetation stem density, given in number per unit area [#/m^2].
                Input can be either one constant value to be applied everywhere,
                or an array giving one value per cell centroid. If None is
                given, we check to see if values were previously defined in
                domain.quantities()

            veg_height : `float or np.ndarray`
                Vegetation stem height, given in meters. Input can be either
                one constant value to be applied everywhere, or an array giving
                one value per cell centroid. If None is given, we check to see
                if values were previously defined in domain.quantities()

            bed_friction_const : `float or np.ndarray`, optional
                Bed friction Chezy coefficient. Default is 65. Input can be either
                one constant value to be applied everywhere, or an array giving
                one value per cell centroid.

        """
        super().__init__(domain, description, label, logging, verbose)

        #-----------------------------------------------------
        # Pull domain information
        self.depth = self.stage_c - self.elev_c  # see base Operator constructor for attributes that get transferred from domain
        self.g = self.domain.g  # g value from ANUGA domain/config

        #-----------------------------------------------------
        # Set constants
        self.bed_friction = np.ones_like(self.depth, dtype=float)*bed_friction_const
        self.K = 0.41  # von Karman constant
        self.Cd = 1.0  # drag coefficient

        #-----------------------------------------------------
        # Initialize vegetation characteristics
        self.set_vegetation(veg_diameter, veg_density, veg_height)
            

    def __call__(self):
        """
        Apply vegetation drag according to veg_diameter and veg_density quantities
        """
        # Get the timestep for explicit update
        self.dt = self.get_timestep()
        
        self.depth = self.stage_c - self.elev_c

        # Only apply where vegetation exists and depth is positive
        self.inds = (self.veg_density.centroid_values > 0) & (self.depth > 0.01) # Update active indices

        self.update_quantities()


    def set_vegetation(self, veg_diameter=None, veg_density=None, veg_height=None):
        """
        Set vegetation characteristics, either for the first time or as an update
        """
        for name, values in zip(['veg_diameter', 'veg_density', 'veg_height'],
                                 [veg_diameter, veg_density, veg_height]):
            self.set_veg_quantity(name, values)

        self.update_coefficients()


    def set_veg_quantity(self, name, values):
        """Register vegetation quantity and set values in the domain."""
        if name not in self.domain.quantities:
            Quantity(self.domain, name=name, register=True)

        q = self.domain.quantities[name]

        if values is not None:
            q.set_values(values, location="centroids")

        setattr(self, name, q)  # always keep the Quantity, not raw array


    def update_coefficients(self):
        """Recompute coefficients after vegetation update."""
        vdiam = self.veg_diameter.centroid_values
        vdens = self.veg_density.centroid_values

        self.a1 = self.bed_friction**-2   # First lumped coefficient
        self.a2 = vdiam*vdens/(2*self.g)  # Second lumped coefficient
        self.a3 = sqrt(self.g)/self.K     # Third lumped coefficient


    def update_quantities(self):
        """
        Calculate the drag that vegetation imparts on the flow
        and update momentum quantities
        """
        if np.any(self.inds):
            # Cut down some variables to just vegetated areas
            depth_w = self.depth[self.inds]
            hv_w = self.veg_height.centroid_values[self.inds]
            Cv_w = np.zeros_like(depth_w) # Vegetated Chezy
            a1_w = self.a1[self.inds]
            a2_w = self.a2[self.inds]
            xmom_w = self.xmom_c[self.inds]
            ymom_w = self.ymom_c[self.inds]

            # calculate discharge in the cell
            qcell_w = sqrt(xmom_w**2 + ymom_w**2)

            # Calculate Chezy
            Cv_w = (1./sqrt(a1_w + a2_w*self.Cd*minimum(depth_w, hv_w)) 
                    + self.a3*log(maximum(depth_w, hv_w)/hv_w))

            # Compute friction slope
            Sf_x = self.g*xmom_w*qcell_w/(Cv_w**2*depth_w**2 + 1e-6)
            Sf_y = self.g*ymom_w*qcell_w/(Cv_w**2*depth_w**2 + 1e-6)

            self.xmom_c[self.inds] = xmom_w - Sf_x*self.dt
            self.ymom_c[self.inds] = ymom_w - Sf_y*self.dt


    def parallel_safe(self):
        """If Operator is applied independently on each cell and
        so is parallel safe.
        """
        return True


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Veg_operator, time '
        message += str(self.get_time())
        return message
