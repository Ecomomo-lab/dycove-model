###############################################################
#  ANUGA_baptist.py
###############################################################

"""
The Baptist_operator class for ANUGA was originally developed by Kyle Wright:

- Wright, K., Passalacqua, P., Simard, M., & Jones, C. E. (2022). Integrating 
  connectivity into hydrodynamic models: An automated open-source method to 
  refine an unstructured mesh using remote sensing. Journal of Advances in 
  Modeling Earth Systems, 14, e2022MS003025. 
  https://doi.org/10.1029/2022MS003025

The Baptist_operator used in DYCOVE has been modified from Wright et al. (2022)
in the following ways:

- Quantity checks have been removed (this class is not part of public API).
- Logic has been distributed into more self-contained methods.
- New methods set_vegetation and set_veg_quantity have been added to allow
  DYCOVE to update vegetation parameters after the Operator has been
  instantiated.

The Baptist formulation itself comes from:

- Baptist, M. J., Babovic, C., Uthurburu, J. R., Uittenbogaard, R. E., Mynett, 
  A., & Verwey, A. (2007). "On inducing equations for vegetation resistance." 
  Journal of Hydraulic Research, 45(4), 435–450.
  https://doi.org/10.1080/00221686.2007.9521778

"""

#from __future__ import division, absolute_import, print_function
import numpy as np


class _OperatorBase:
    """Fallback base class when ANUGA is unavailable."""
    pass

def _import_anuga():
    try:
        from anuga import Quantity
        from anuga.operators.base_operator import Operator        
        return Quantity, Operator
    except ImportError:
        msg = ("The `anuga` package is not installed. "
               "Refer to the documentation for installation instructions.")
        raise ImportError(msg)
    

#===============================================================================
# Vegetation operator applying drag to the flow
#===============================================================================
class Baptist_operator(_OperatorBase):
    """
    Applies vegetation-induced flow drag based on the Baptist formulation.

    This operator modifies the local flow field in ANUGA simulations to account for
    emergent and submerged vegetation effects. The vegetation-adjusted Chezy coefficient
    :math:`C_v` is computed as:

    .. math::

        C_v = \\left( C_b^{-2} + \\frac{C_d m D}{2 g} \\min(h, h_v) \\right)^{-0.5}
              + \\frac{\\sqrt{g}}{k} \\ln\\left( \\frac{\\max(h, h_v)}{h_v} \\right)

    where:

    - :math:`C_v` — vegetated Chezy coefficient  
    - :math:`C_b` — bed Chezy coefficient  
    - :math:`C_d` — drag coefficient  
    - :math:`m` — stem density :math:`[m^{-2}]`  
    - :math:`D` — stem diameter :math:`[m]`  
    - :math:`g` — gravitational acceleration :math:`[m/s^2]`  
    - :math:`h` — flow depth :math:`[m]`  
    - :math:`h_v` — vegetation height :math:`[m]`  
    - :math:`k` — von Kármán constant  

    For computational efficiency, some terms are precomputed in the initializer, giving:

    .. math::

        C_v = (a_1 + a_2 C_d \\min(h, h_v))^{-0.5} + a_3 \\ln\\left( \\frac{\\max(h, h_v)}{h_v} \\right)

    The operator updates momentum explicitly using:

    .. math::

        \\frac{\\partial (uh)}{\\partial t} = 
        - g \\frac{uh \\sqrt{u^2 + v^2}}{C_v^2 h^2}

    Parameters
    ----------
    domain : anuga.Domain
        ANUGA domain instance.

    veg_diameter : float or numpy.ndarray
        Vegetation stem diameter [m]. Either a single constant applied everywhere or
        an array with one value per cell centroid.

    veg_density : float or numpy.ndarray
        Vegetation stem density [#/m²]. Either a single constant applied everywhere or
        an array with one value per cell centroid.

    veg_height : float or numpy.ndarray
        Vegetation stem height [m]. Either a single constant applied everywhere or
        an array with one value per cell centroid.

    bed_friction_const : float or numpy.ndarray, optional
        Bed friction Chezy coefficient. Default is 65. Either a single constant applied
        everywhere or an array with one value per cell centroid.
    """

    def __init__(self, 
                 domain, 
                 veg_diameter=None,
                 veg_density=None,
                 veg_height=None,
                 drag=1.0,
                 bed_friction_const=65.0,
                 description = None,
                 label = None,
                 logging = False,
                 verbose = False):

        # ensure lazy import the first time the object is created
        if not hasattr(self, "_anuga_loaded"):
            self.Quantity, Operator = _import_anuga()
            # Patch the base class of Baptist_operator AFTER module load
            Baptist_operator.__bases__ = (Operator,)
            self._anuga_loaded = True

        super().__init__(domain, description, label, logging, verbose)

        #-----------------------------------------------------
        # Pull domain information
        # See base Operator constructor for attributes that get transferred from domain
        self.depth = self.stage_c - self.elev_c
        self.g = self.domain.g  # g value from ANUGA domain/config

        #-----------------------------------------------------
        # Set constants
        self.bed_friction = np.ones_like(self.depth, dtype=float)*bed_friction_const
        self.K = 0.41   # von Karman constant
        self.Cd = drag  # drag coefficient C_d

        #-----------------------------------------------------
        # Initialize vegetation characteristics
        self.set_vegetation(veg_diameter, veg_density, veg_height)
            

    def __call__(self):
        """
        Get the current water depth and vegetated cell indicies, 
        then update Chezy values accordingly.
        """
        # Get the timestep for explicit update
        self.dt = self.get_timestep()
        
        self.depth = self.stage_c - self.elev_c

        # Only apply where vegetation exists and depth is positive
        self.inds = (self.veg_density.centroid_values > 0) & (self.depth > 0.01) # Update active indices

        self.update_quantities()


    def set_vegetation(self, veg_diameter=None, veg_density=None, veg_height=None):
        """ Set vegetation characteristics, either for the first time or as an update. """
        for name, values in zip(['veg_diameter', 'veg_density', 'veg_height'],
                                 [veg_diameter, veg_density, veg_height]):
            self.set_veg_quantity(name, values)

        self.update_coefficients()


    def set_veg_quantity(self, name, values):
        """ Register vegetation quantity and set values in the domain. """
        if name not in self.domain.quantities:
            self.Quantity(self.domain, name=name, register=True)

        q = self.domain.quantities[name]

        if values is not None:
            q.set_values(values, location="centroids")

        setattr(self, name, q)  # always keep the Quantity, not raw array


    def update_coefficients(self):
        """ Recompute Baptist coefficients after vegetation update. """
        vdiam = self.veg_diameter.centroid_values
        vdens = self.veg_density.centroid_values

        self.a1 = self.bed_friction**-2   # First lumped coefficient
        self.a2 = vdiam*vdens/(2*self.g)  # Second lumped coefficient
        self.a3 = np.sqrt(self.g)/self.K     # Third lumped coefficient


    def update_quantities(self):
        """
        Calculate drag that vegetation imparts on flow, then update momentum quantities.
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
            qcell_w = np.sqrt(xmom_w**2 + ymom_w**2)

            # Calculate Chezy
            Cv_w = (1./np.sqrt(a1_w + a2_w*self.Cd*np.minimum(depth_w, hv_w)) 
                    + self.a3*np.log(np.maximum(depth_w, hv_w)/hv_w))

            # Compute friction slope
            Sf_x = self.g*xmom_w*qcell_w/(Cv_w**2*depth_w**2 + 1e-6)
            Sf_y = self.g*ymom_w*qcell_w/(Cv_w**2*depth_w**2 + 1e-6)

            self.xmom_c[self.inds] = xmom_w - Sf_x*self.dt
            self.ymom_c[self.inds] = ymom_w - Sf_y*self.dt


    def parallel_safe(self):
        """ If Operator is applied independently on each cell and so is parallel safe. """
        return True


    def timestepping_statistics(self):
        from anuga import indent

        message  = indent + self.label + ': Veg_operator, time '
        message += str(self.get_time())
        return message
