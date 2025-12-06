"""
A wrapper class for developing a simple rectangular ANUGA domain.

Used for `tide_channel` example.

Creating a domain in parallel requires a few additional lines of code, namely:

- myid: anuga processor ID reference, for limiting domain creation to core 0
- distribute: anuga function to decompose the domain across all processors

"""

import numpy as np
import anuga
from anuga import myid, distribute


class RectangSlopeDomainGenerator:
    """
    Generate an ANUGA rectangular coastal domain with a simple sloped beach,
    a tidal channel, and a lagoon behind a dune.

    This utility class wraps the ANUGA mesh generator and assigns
    topography, boundary conditions, friction, and initial stage.  
    All topography-related parameters are optional and default values used for
    the `tide_channel` example.

    Parameters
    ----------
    model_name : str
        Base name for the ANUGA model output file.
    rectang_dims : tuple of float, optional
        Domain dimensions ``(Lx, Ly)`` in meters.
    mesh_spacing : float, optional
        Target cell size (m) used in `rectangular_cross_domain`.
    min_elev : float, optional
        Elevation at the offshore boundary (left boundary).
    dune_elev : float, optional
        Elevation of the dune crest in meters.
    dune_location : float, optional
        Location of the dune crest as a fraction of domain length
        (e.g., 0.65 = 65% of the domain in x).
    channel_bot_elev : float, optional
        Bottom elevation of the tidal channel.
    lagoon_bot_elev : float, optional
        Minimum elevation in the lagoon.
    channel_width_frac : float, optional
        Fraction of the *y*-dimension occupied by the channel width.
    tide_props : dict, optional
        Dictionary containing tidal forcing parameters:
        ``{"amplitude": A, "period": T, "phase": deg, "MWL": offset}``.
    WL_0 : float, optional
        Initial water level on the seaward side.
    mannings_n : float, optional
        Spatially constant Manning's `n` friction coefficient.
    plotting : bool, optional
        Whether to plot mesh and mesh elevation using ANUGA's 
        ``Domain_plotter``.

    Notes
    -----
    The topography is constructed as a 1D slope in ``x`` with a local horizontal
    channel depression applied. The channel is formed by setting elevations 
    within a given ``channel_width_frac`` band to a constant ``channel_bot_elev``, 
    overriding the sloped surface.

    """

    def __init__(
        self,
        model_name,
        rectang_dims=(1000.0, 600.0),
        mesh_spacing=40.0,
        min_elev=-0.5,
        dune_elev=0.5,
        dune_location=0.65,
        channel_bot_elev=0.05,
        lagoon_bot_elev=-0.5,
        channel_width_frac=0.1,
        tide_props=None,
        WL_0=0.0,
        mannings_n=0.025,
        plotting=False,
    ):

        # If tide_props wasn't supplied, use defaults
        if tide_props is None:
            tide_props = {
                "amplitude": 0.4,
                "period": 12.0 * 3600,
                "MWL": 0.0,
            }

        self.model_name = model_name
        self.dims = rectang_dims
        self.dx = mesh_spacing
        self.min_z = min_elev
        self.dune_z = dune_elev
        self.dune_loc = dune_location
        self.chan_bot_z = channel_bot_elev
        self.lago_bot_z = lagoon_bot_elev
        self.chan_frac = channel_width_frac
        self.tide = tide_props
        self.friction = mannings_n
        self.WL_0 = WL_0
        self.plotting = plotting

        if myid == 0:
            self.create_domain()
            self.set_initial_quantities()
        self.parallelize_domain()
        self.set_boundary_conditions()

    
    def create_domain(self):
        # Create a rectangular domain with named boundaries "left", "right", "top" and "bottom"
        self.domain = anuga.rectangular_cross_domain(self.dims[0]/self.dx,  # no. of cells in x-direction
                                                     self.dims[1]/self.dx,  # no. of cells in y-direction
                                                     len1=self.dims[0],     # length in x-direction
                                                     len2=self.dims[1])     # length in y-direction
        self.domain.set_name(self.model_name)  # ANUGA output name for .sww file
        self.domain.set_low_froude(1)          # Always set to 1 for low-froude flows
        self.domain.set_flow_algorithm('DE0')  # For stable solution, can always try 'DE1' but is twice as slow
        self.domain.set_CFL(0.9)               # A slightly lower CFL number combined with DE0 is usually sufficient
        self.domain.set_minimum_storable_height(0.01)  # Only store heights > 1 cm

        print(self.domain.statistics(), flush=True)
        
        
    def set_initial_quantities(self):        
        self.domain.set_quantity('elevation', self.topography)
        self.domain.set_quantity('friction', self.friction)

            
    def set_boundary_conditions(self):
        # Set initial stage at all cells, but only on the ocean side
        x = self.domain.quantities['x'].centroid_values
        topo = self.domain.quantities['elevation'].centroid_values
        stage = np.where(x < self.dims[0]*self.dune_loc, self.WL_0, topo)  # left of the dune only
        self.domain.set_quantity('stage', stage, location='centroids')

        # Solid reflective wall boundary
        Br = anuga.Reflective_boundary(self.domain)

        # Tidal wave function
        def tidal_func(
                t, 
                A=self.tide["amplitude"], 
                T=self.tide["period"],
                # Phase specifed in deg, converted to rad below
                phase=self.tide["phase"] if "phase" in self.tide else 0,
                offset=self.tide["MWL"]
                ):
            return A*np.sin(2*np.pi * t/T + phase*np.pi/180) + offset
        # Time varying (tidal) boundary condition
        Bt = anuga.Time_boundary(self.domain, function=lambda t: [tidal_func(t), 0.0, 0.0])

        # Map boundary conditions to domain boundaries (names come from rectangular_cross_domain)
        self.domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})


    def topography(self, x, y):
        """
        Returns bed elevation (z) for a simple sloped beach with a central tidal channel 
        connecting to a lagoon behind a dune.
        """

        # position of dune crest
        x_dune = self.dune_loc * self.dims[0]

        # compute slope
        beach_slope = (self.dune_z - self.min_z) / x_dune
        back_slope = beach_slope * 2  # steeper slope behind dune
        # base profile (1D)
        z = np.where(
            x <= x_dune,
            self.min_z + beach_slope * x,  # beach rising toward dune
            np.maximum(self.dune_z - back_slope * (x - x_dune), self.lago_bot_z),  # sloping down toward lagoon
        )

        # define channel depression (Gaussian in y)
        centerline_y = self.dims[1] / 2.
        chan_halfwidth = self.chan_frac*self.dims[1] / 2.
        channel = np.abs(y - centerline_y) <= chan_halfwidth  # boolean mask

        # lower to channel bottom only where terrain is above z_channel
        z[channel & (z > self.chan_bot_z)] = self.chan_bot_z

        return z


    def parallelize_domain(self):
        # set domain to None for all cores except the first
        if myid != 0:
            self.domain = None
        # distribute domain fro the first core to all cores
        self.domain = distribute(self.domain)