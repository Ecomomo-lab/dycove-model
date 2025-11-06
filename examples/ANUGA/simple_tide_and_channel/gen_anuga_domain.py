"""
A wrapper class for developing a simple rectangular ANUGA domain.

Used for `simple_tide` example.
"""

import numpy as np
import matplotlib.pyplot as plt
import anuga
from dycove.utils.log import Reporter

r = Reporter()

class RectangSlopeDomainGenerator:
    def __init__(self, 
                 model_name,
                 rectang_dims,
                 mesh_spacing,
                 min_elev,
                 dune_elev,
                 dune_location,
                 channel_bot_elev,
                 lagoon_bot_elev,
                 channel_width_frac,
                 tide_props,
                 WL_0=0.0,
                 mannings_n=0.025,
                 plotting=False,
                 ):
        
        self.model_name = model_name
        self.dims = rectang_dims
        self.dx = mesh_spacing
        self.min_z = min_elev
        self.dune_z = dune_elev
        self.dune_loc = dune_location
        self.chan_bot_z = channel_bot_elev
        self.lago_bot_z = lagoon_bot_elev
        self.chan_frac = channel_width_frac
        self.friction = mannings_n
        self.tide = tide_props
        self.WL_0 = WL_0

        self.plotting = plotting

        self.create_domain()
        self.set_initial_quantities()
        self.set_boundary_conditions()
    
    
    def create_domain(self):

        r.report("Creating ANUGA domain and mesh...")

        # Create a domain with named boundaries "left", "right", "top" and "bottom"
        self.domain = anuga.rectangular_cross_domain(self.dims[0]/self.dx, self.dims[1]/self.dx, 
                                                     len1=self.dims[0], len2=self.dims[1])
        self.domain.set_name(self.model_name)  # ANUGA output name
        self.domain.set_low_froude(1)          # always set to 1 for low-froude flows
        self.domain.set_flow_algorithm('DE1')      # for stable solution, can always try 'DE0' for faster results
        self.domain.set_minimum_allowed_height(0.01)  # Only store heights > 1 cm

        r.report(self.domain.statistics())

        if self.plotting:
            plt.figure()
            self.dplotter = anuga.Domain_plotter(self.domain)  
            plt.triplot(self.dplotter.triang, linewidth = 0.1)
            plt.axis('scaled')
            plt.savefig('mesh.png', bbox_inches='tight')
            plt.close()
        
        
    def set_initial_quantities(self):        
        r.report("Assigning elevation and friction data to mesh...\n")

        self.domain.set_quantity('elevation', self.topography)
        self.domain.set_quantity('friction', self.friction)

        if self.plotting:
            plt.figure()
            plt.tripcolor(self.dplotter.triang, 
                          facecolors = self.dplotter.elev,
                          vmin=-0.5, vmax=0.5,
                          cmap='bone')

            cbar = plt.colorbar()
            cbar.set_label('Mesh topography [m, NAVD88]', rotation=270, labelpad=20)
            plt.axis('scaled')

            plt.savefig('mesh_elev.png', bbox_inches='tight')
            plt.close()

            
    def set_boundary_conditions(self):
        r.report("Setting boundary conditions...\n")

        # set initial stage, but only on the ocean side
        x = self.domain.quantities['x'].centroid_values
        topo = self.domain.quantities['elevation'].centroid_values
        stage = np.where(x < self.dims[0]*self.dune_loc, self.WL_0, topo)
        self.domain.set_quantity('stage', stage, location='centroids')

        Br = anuga.Reflective_boundary(self.domain)  # Solid reflective wall

        # tidal wave function
        def tidal_func(
                t, 
                A=self.tide["amplitude"], 
                T=self.tide["period"], 
                phase=self.tide["phase"] if "phase" in self.tide else 0, # phase specifed in deg, converted to rad below
                offset=self.tide["MWL"]
                ):
            return A*np.sin(2*np.pi * t/T + phase*np.pi/180) + offset
        Bt = anuga.Time_boundary(self.domain, function=lambda t: [tidal_func(t), 0.0, 0.0])


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
