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
                 slope, 
                 mannings_n,
                 tide_props,
                 WL_0=0.0,
                 ):
        
        self.model_name = model_name
        self.dims = rectang_dims
        self.dx = mesh_spacing
        self.min_z = min_elev
        self.slope = slope
        self.friction = mannings_n
        self.tide = tide_props
        self.WL_0 = WL_0

        self.plotting = True

        self.create_anuga_domain()
        self.set_initial_quantities()
        self.set_boundary_conditions()
    

    def create_anuga_domain(self):

        r.report("Creating ANUGA domain and mesh...")

        # Create a domain with named boundaries "left", "right", "top" and "bottom"
        self.domain = anuga.rectangular_cross_domain(self.dims[0]/self.dx, self.dims[1]/self.dx, 
                                                     len1=self.dims[0], len2=self.dims[1])
        self.domain.set_name(self.model_name)  # ANUGA output name
        self.domain.set_low_froude(1)          # always set to 1 for low-froude flows

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

        def topography(x, y):
            return -self.min_z + x*self.slope  # linear bed slope between -0.5 and 0.5

        self.domain.set_quantity('elevation', topography)
        self.domain.set_quantity('friction', self.friction)

        if self.plotting:
            plt.figure()
            plt.tripcolor(self.dplotter.triang, 
                          facecolors = self.dplotter.elev,
                          vmin=-1, vmax=1,
                          cmap='bone')

            cbar = plt.colorbar()
            cbar.set_label('Mesh topography [m, NAVD88]', rotation=270, labelpad=20)
            plt.axis('scaled')

            plt.savefig('mesh_elev.png', bbox_inches='tight')
            plt.close()

            
    def set_boundary_conditions(self):
        r.report("Setting boundary conditions...\n")

        # set initial stage
        self.domain.set_quantity("stage", self.WL_0)

        Br = anuga.Reflective_boundary(self.domain)  # Solid reflective wall

        # tidal wave function
        def tidal_func(
                t, 
                A=self.tide["amplitude"], 
                T=self.tide["period"], 
                # phase to be specifed in degrees, converted to radians below
                phase=self.tide["phase"] if "phase" in self.tide else 0,
                offset=self.tide["MWL"]
                ):
            return A*np.sin(2*np.pi * t/T + phase*np.pi/180) + offset
        Bt = anuga.Time_boundary(self.domain, function=lambda t: [tidal_func(t), 0.0, 0.0])


        self.domain.set_boundary({'left': Bt, 'right': Br, 'top': Br, 'bottom': Br})

