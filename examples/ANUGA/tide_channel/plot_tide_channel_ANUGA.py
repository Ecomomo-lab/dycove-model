from pathlib import Path
from dycove import plotting


"""
Quantites to plot from numerical model (use exact names):
'Bathymetry', 
'WSE', 
'Depth', 
'Velocity', 
'Max Shear Stress',

Vegetation quantites to plot from DYCOVE (use exact names):
'Stem Height', 
'Stem Diameter', 
'Stem Density', 
'Fractions',
'Potential Mortality -- Flooding',
'Potential Mortality -- Desiccation',
'Potential Mortality -- Uprooting',
'Potential Mortality -- Burial',
'Potential Mortality -- Scour',
'Mortality -- Flooding',
'Mortality -- Desiccation',
'Mortality -- Uprooting',
'Mortality -- Burial',
'Mortality -- Scour',
'Mortality -- Total'
"""

plotter = plotting.ModelPlotter(
    simdir = Path('.'),
    #quantity = 'Velocity',
    quantity = 'Stem Height',
    #quantity = 'Mortality -- Total',
    #quantity = 'Fractions',
    plot_times = {  # times specified here are hydrodynamic time, not eco-morpho time
        # sim hr to start plotting
        'plotHR_0': 0*24.,
        # sim hr to stop plotting, not to exceed total sim length.
        'plotHR_f': 7*24.,  # 21 hydro days ~ 3 eco-morpho years when vegfac ~ 50
        # sim hrs between map outputs, default for ANUGA, value for DFM given in MDU file
        'mapHR_int': 1,
        # hrs between consecutive plots, cannot be less than map_output, unused if plotting vegetation
        'plotHR_int': 1,
        },
    cmap_lims = {
        'Bathymetry': (-0.5, 0.5),
        'Velocity': (0, 0.3),
        },
    animate=False,
)

plotter.run()