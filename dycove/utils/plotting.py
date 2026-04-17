"""
Class for plotting various output quantites from DYCOVE-ANUGA/DFM
"""

from typing import Optional, Any, Union
import numpy as np
import math
from pathlib import Path
import json
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, ListedColormap
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib import path as mpath
from scipy.spatial import cKDTree

from dycove.utils.array_math import cell_averaging, sum_product
from dycove.utils.model_loader import DFMMapLoader, ANUGAMapLoader
from dycove.utils.simulation_reporting import Reporter

r = Reporter()

    
def create_nn_interpFunc(x_coords, y_coords, grid_size, k_nn=1, 
                         polygon_csv=None, extents=None,
                        ):
    """
    Creates a nearest-neighbor interpolator function for a fixed grid.
    
    Parameters
    ----------
    coords : list
        List of (x, y) coordinate pairs.
    grid_size : int
        Desired grid cell size for interpolation.
    k_nn : int
        Number of nearest neighbors to use (typically 1 or 3).
    polygon_csv : Path or str
        Path to csv file containing polygon vertices (columns: x, y). If provided, 
        values outside the polygon will be masked as NaN.
    extents : list or tuple
        Manual x- and y-extent limits for the interpolation grid.
        Example: extents=(10, 50, 20, 80).
        If None, defaults to the min/max of the input coords.

    Returns
    -------
    Function that can be used to interpolate new z values efficiently.

    """

    # Convert input coordinates to NumPy array
    x_coords, y_coords = np.asarray(x_coords), np.asarray(y_coords)

    # Set extents to full mesh extent if smaller window not provided
    if extents is None:
        extents = (x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max())

    ext_inds = np.where((extents[0] <= x_coords) & (x_coords <= extents[1]) &
                        (extents[2] <= y_coords) & (y_coords <= extents[3]))[0]
    
    # Filter coordinates to those within extents
    x_coords, y_coords = x_coords[ext_inds], y_coords[ext_inds]

    # Derive grid limits from coords, which have already been filtered to the desired extent
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()

    # Create the grid
    x_grid = np.arange(x_min, x_max + grid_size, grid_size)
    y_grid = np.arange(y_min, y_max + grid_size, grid_size)
    X, Y = np.meshgrid(x_grid, y_grid, indexing="xy")
    grid_points = np.column_stack([X.ravel(), Y.ravel()])  # Flatten grid

    # Build KDTree and find nearest-neighbor indices
    tree = cKDTree([(x, y) for x, y in zip(x_coords, y_coords)])
    distances, tree_inds = tree.query(grid_points, k=k_nn, workers=-1)  # Fast parallel query

    # ---- Polygon mask (optional) ---- #
    if polygon_csv is None:
        inside_mask = np.ones(X.shape, dtype=bool)
    else:
        # Load polygon vertices (supports headers)
        poly_data = np.loadtxt(polygon_csv, delimiter=",")
        poly_x, poly_y = poly_data[:, 0], poly_data[:, 1]

        poly_path = mpath.Path(np.column_stack((poly_x, poly_y)))
        inside_mask = poly_path.contains_points(grid_points).reshape(X.shape)
        

    def nn_interpolator(z_values):
        """
        Interpolates z-values using k-nearest neighbors.

        Parameters
        ----------
        z_values: np.ndarray or list
            Array of shape (L,) where L = len(coords).

        Returns
        -------
        2D array of interpolated values matching the grid shape.

        """

        # Filter z_values to match extents
        z_values = np.asarray(z_values)[ext_inds]
        assert len(z_values) == len(x_coords), "Mismatch in z-values length"

        if k_nn == 1:
            Z_interp = z_values[tree_inds]
        else:
            # Weight by inverse distance, avoiding divide-by-zero
            with np.errstate(divide="ignore"):
                weights = 1./distances
            weights[np.isinf(weights)] = 1e10  # handle exact matches
            weights /= np.sum(weights, axis=1, keepdims=True)
            Z_interp = np.sum(z_values[tree_inds]*weights, axis=1)

        Z = Z_interp.reshape(Y.shape)

        # Apply mask
        Z[~inside_mask] = np.nan

        # Flip vertically to match typical image orientation
        Z = np.flipud(Z)

        return Z

    return nn_interpolator


class ModelPlotter:
    """
    Plot hydrodynamic and vegetation quantities from DYCOVE model runs (ANUGA or DFM).

    This class handles loading, interpolating, and plotting 2-D spatial quantities
    (e.g., depth, velocity, vegetation fractions, mortality, etc.) from DYCOVE simulations.

    Parameters
    ----------
    simdir : str or Path
        Path to the root simulation directory containing either a `dflowfm/` directory 
        (DFM) or a ``.sww`` model file (ANUGA).
    quantity : str
        Quantity to plot. Must be one of the supported quantities listed below.
    plot_times : dict[str, int]
        Dictionary of plot timing parameters with the following keys:
        
        - ``'plotHR_0'``: Sim. hour to start plotting
        - ``'plotHR_f'``: Sim. hour to end plotting
        - ``'mapHR_int'``: Sim. hours between hydrodynamic model outputs
        - ``'plotHR_int'``: Sim. hours between consecutive plots, cannot be less than 
          map_output, unused if plotting vegetation.

    plot_separate_species : bool, optional
        If True, vegetation quantities are plotted separately for each species (multi-species only).
        Vegetation fractions are plotted by species, rather than by cohort.
    plot_method : str, optional
        ``'interp'`` (default) for nearest-neighbor interpolation to a regular grid, or
        ``'exact'`` for plotting on the original mesh nodes/cells (must be regular grid). NOT WORKING CURRENTLY.
    plot_interp_vars : dict, optional
        Options for interpolating to 2D grid, see keys and details under __init__ below.
    plot_specs : dict, optional
        Plot appearance options, e.g., see keys and details under __init__ below.
    cmap_lims : dict, optional
        Colorbar limits for each quantity, see keys and details under __init__ below.
    cmaps : dict, optional
        Colormaps for each quantity. Accepts Matplotlib ``Colormap`` objects. See keys and details 
        under __init__ below.
    quantity_units : dict, optional
        Unit labels, scaling factors, and minimum plotting thresholds for each quantity, 
        see keys and details under __init__ below.
    plot_vectors : bool, optional
        If True, plots normalized velocity vectors (quantity=="Velocity" only).
    vector_props : dict, optional
        Scaling and spacing options for velocity vectors, used if quantity == "Velocity",
        see keys and details under __init__ below.
    show_title : bool, optional
        If True, include auto-generated plot title at top of plot.
    show_topo_cbar : bool, optional
        If True, include elevation color bar on the left-hand side.
    show_axis_ticks : bool, optional
        If True, axis ticks will be shown on left and bottom sides of image, in data units.
    scalebar : bool, optional
        If True, include a scalebar object on the plot.
    scalebar_props : dict, optional
        Scalebar properties, e.g., see keys and details under __init__ below.
    mask_bndy_file : str or Path, optional
        Path to polygon CSV file for masking interpolation outside the domain. Probably
        required if model domain is not rectangular.
    extents : tuple[float, float, float, float], optional
        Geographic extents for interpolation (xmin, xmax, ymin, ymax). Useful for
        zooming in on a particular location. NOTE: ymin and ymax are in MODEL coordinates, not 
        numpy coordinates.
    animate : bool, optional
        If True, creates an animated GIF after generating all plots. Requires ``imageio`` library.
    delete_static_imgs : bool, optional
        If True, deletes individual PNGs after animation is built.
    save_grids : bool, optional
        If True, saves interpolated grids (``.npz``) alongside figures. Can be helpful for
        external plotting of, say, difference maps, which are not supported here (yet).

    Supported Quantities (use exact names)
    -------------------------------------
    Numerical model quantities:

    - ``'Bathymetry'``
    - ``'WSE'``
    - ``'Depth'``
    - ``'Velocity'``
    - ``'Max Shear Stress'``

    DYCOVE quantities:

    - ``'Stem Height'``
    - ``'Stem Diameter'``
    - ``'Stem Density'``
    - ``'Fractions'``
    - ``'Potential Mortality -- Flooding'``
    - ``'Potential Mortality -- Desiccation'``
    - ``'Potential Mortality -- Uprooting'``
    - ``'Potential Mortality -- Burial'``
    - ``'Potential Mortality -- Scour'``
    - ``'Mortality -- Flooding'``
    - ``'Mortality -- Desiccation'``
    - ``'Mortality -- Uprooting'``
    - ``'Mortality -- Burial'``
    - ``'Mortality -- Scour'``
    - ``'Mortality -- Total'``

    Notes
    -----
    - All optional arguments have default values, but some dictionary arguments
      will almost definitely need to be changed or provided by the user,
      such as ``cmap_lims``, for which appropriate values will depend on model data.
    - The class automatically detects the model type (``'DFM'`` or ``'ANUGA'``)
      based on directory contents.
    - Uses :class:`~dycove.utils.model_loader.DFMMapLoader` or 
      :class:`~dycove.utils.model_loader.ANUGAMapLoader` for reading model-specific 
      outputs.
    - Vegetation quantities (e.g. fractions, stem attributes, mortality)
      are handled differently than hydrodynamic quantities, using the ``ecofac`` 
      time scaling factor.
    - Quantities that include ``'Mortality'`` may be plotted either as potential
      (model-estimated stress exposure) or applied (post-threshold mortality)
      fields. See :class:`~dycove.sim.vegetation_data.VegCohort` documentation and 
      ``mortality*`` methods of
      :class:`~dycove.sim.vegetation.VegetationSpecies` for details.

    """


    def __init__(self, 
                 simdir, 
                 quantity, 
                 plot_times,
                 plot_separate_species = True,
                 plot_method = "interp",  # "exact" not working currently
                 plot_interp_vars: Optional[dict[str, int]] = None,
                 plot_specs: Optional[dict[str, Any]] = None, 
                 cmap_lims: Optional[dict[str, tuple[float, float]]] = None, 
                 cmaps: Optional[dict[str, Colormap]] = None,
                 quantity_units: Optional[dict[str, tuple[str, Union[int, float]]]] = None, 
                 plot_vectors = False,
                 vector_props: Optional[dict[str, Any]] = None,
                 show_title = True,
                 show_topo_cbar = False,
                 show_axis_ticks = False,
                 scalebar = False,
                 scalebar_props: Optional[dict[str, Any]] = None,
                 mask_bndy_file: Optional[Union[Path, str]] = None,
                 extents: Optional[tuple[float, float, float, float]] = None,
                 animate = False, 
                 delete_static_imgs = False,
                 save_grids = False,
                 ):
        
        self.simdir = Path(simdir)
        self.quantity = quantity
        self.plot_times = plot_times
        self.plot_separate_species = plot_separate_species
        self.plot_method = plot_method
        self.plot_vectors = plot_vectors
        self.show_title = show_title
        self.show_topo_cbar = show_topo_cbar
        self.show_axis_ticks = show_axis_ticks
        self.scalebar = scalebar
        self.mask_bndy_file = mask_bndy_file
        self.extents = extents
        self.animate = animate
        self.delete_static_imgs = delete_static_imgs
        self.save_grids = save_grids

        # Simplify quantity name if "Mortality"
        self.full_quantity_name = self.quantity
        if "Mortality" in self.full_quantity_name.split():
            self.quantity = "Mortality"

        self.output_plot_dir = self.simdir / "figures" / self.full_quantity_name.replace(" ", "")
        self.model_type = self.get_model_type()
        self.modeldir = self.get_modeldir()
        self.model_name = self.get_model_name()
        self.eco_plot = self.quantity in ["Fractions", "Stem Height", "Stem Density", "Stem Diameter", "Mortality"]

        # Change how we label time in our plots depending on if we are plotting vegetation (uses ecofac) or hydrodynamics (no morfac/ecofac)
        self.plot_label_time = "eco-morphodynamic" if self.eco_plot else "hydrodynamic"
            
        default_plot_interp_vars = {
            "cell_size": 5,    # inteprolated grid size, in meters
            "n_neighbors": 1,  # neighbors used for interpolation, typically 1 or 3
        }
        self.plot_interp_vars = {**default_plot_interp_vars, **(plot_interp_vars or {})}  # merge provided custom values with default values

        default_plot_specs = {
            "figsize": (6, 6),  # figure image dimensions
            "fontsize": 12,     # font size for all major text
            "output_dpi": 100,  # resolution of output PNG
            "aspect": "equal",  # "set_aspect()" arg, set to 'auto' to squish grid to fit figsize, 
        }                       # or provide float: e.g., 2 means height is stretched by factor of 2, 0.5 is the opposite
        self.plot_specs = {**default_plot_specs, **(plot_specs or {})}  # merge provided custom values with default values

        # Colorbar plotting limits, check quantity_units dictionary for unit consistency
        default_cmap_lims = { # (min, max)
            "Bathymetry"      : (0, 1),
            "WSE"             : (0, 1),
            "Depth"           : (0, 3),
            "Velocity"        : (0, 0.5),
            "Max Shear Stress": (0, 0.1),
            "Fractions"       : (0, 1),
            "Stem Height"     : (0, 2),
            "Stem Density"    : (0, 300),
            "Stem Diameter"   : (0, 50),
            "Mortality"       : (0, 100),
        }
        self.cmap_lims = {**default_cmap_lims, **(cmap_lims or {})}  # merge provided custom values with default values

        default_cmaps = {
            "Bathymetry"      : ListedColormap(plt.colormaps["Greys_r"](np.linspace(0, 1.0, 256))),
            #"Bathymetry"      : ListedColormap(plt.colormaps["bone"](np.linspace(0, 1.0, 256))),
            "WSE"             : ListedColormap(plt.colormaps["Blues_r"](np.linspace(0.1, 1.0, 256))),
            "Depth"           : ListedColormap(plt.colormaps["Blues"](np.linspace(0.1, 1.0, 256))),
            "Velocity"        : ListedColormap(plt.colormaps["viridis"](np.linspace(0, 1.0, 256))),
            "Max Shear Stress": ListedColormap(plt.colormaps["plasma"](np.linspace(0, 1.0, 256))),
            "Fractions"       : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            "Stem Height"     : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            "Stem Density"    : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            "Stem Diameter"   : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            "Mortality" : ListedColormap(plt.colormaps["Oranges"](np.linspace(0.2, 1.0, 256))),
        }
        self.cmaps = {**default_cmaps, **(cmaps or {})}  # merge provided custom values with default values

        # First value: units to use for plot title/colorbar
        # Second value: multiplier based on most standard units (meters, seconds, Newtons, fractions, etc)
        # Third value: minimum value for plotting, anything below this value is masked out
        # examples given for some common conversions
        # Note that the third value for Depth is used to mask out other quantities, e.g. Velocity and WSE
        # check cmap_lims dictionary to ensure units are consistent with colorbar plotting limits
        default_quantity_units = {
            "Bathymetry"      : ("Elevation [m]", 1, -99),  # use ("[ft]", 3.281, -99) for elevations in feet
            "WSE"             : ("WSE [m]", 1, -99),
            "Depth"           : ("Depth [m]", 1, 0.01),
            "Velocity"        : ("Velocity [m/s]", 1, 0),  # use ("[cm/s]", 100, 0) for velocities in cm/s
            "Max Shear Stress": ("Shear Stress [Pa]", 1, 0.001),
            "Fractions"       : ("Fraction [--]", 1, 0.01),
            "Stem Height"     : ("Height [m]", 1, 0.05),
            "Stem Density"    : ("Density [m$^{-2}$]", 1, 5),
            "Stem Diameter"   : ("Diameter [mm]", 1000, 1),  # use ("[m]", 1, 0.001) for diameter in meters
            "Mortality"       : ("Mortality [%]", 100, 1),  # use ("[--]", 1, 0.01) for fraction between 0 and 1
        }
        self.quantity_units = {**default_quantity_units, **(quantity_units or {})}  # merge provided custom values with default values

        default_vector_props = {
            "vect_spacing": 50,  # absolute spacing in meters. If x-y spacing should be different, this must be a tuple giving (x, y) spacing
            "color": "red",      
            "scale": 30,         # see matplotlib.axes.Axes.quiver
            "pivot": "mid",      # see matplotlib.axes.Axes.quiver
            "width": 0.003,      # see matplotlib.axes.Axes.quiver
        }
        self.vector_props = {**default_vector_props, **(vector_props or {})}  # merge provided custom values with default values

        default_scalebar_props = {
            "distance": "200 m",   # scalebar distance to display
            "loc": "upper left",   # position
            "frameon": True,       # add frame (usually looks better)
        }
        self.scalebar_props = {**default_scalebar_props, **(scalebar_props or {})}  # merge provided custom values with default values

        self.read_assign_eco_vars()
        self.setup_model_loaders()
        self.check_inputs()

        # If plot_method is "interp", will build interpolation function for converting 1-D arrays to a grid
        self.interp_func = None
        # If plot_method is "exact", will build empty grid and x/y indices for mapping 1-D arrays to grid
        self.empty_grid, self.x_idx, self.y_idx = None, None, None

        # To save for building animations later
        self.timestrings = []
        self.img_paths = []


    def read_assign_eco_vars(self):
        metadata_file = self.modeldir / "veg_output" / "_eco_time_vars.json"
        # These variables are assigned if the file exists, and are protected by "if self.eco_plot" statements
        # TODO: This may be brittle, consider refactoring.
        if metadata_file.exists():
            with open(metadata_file, "r") as f:
                eco_vars = json.load(f)
            self.n_ets = eco_vars["n_ets"]
            self.veg_int_hr = int(eco_vars["veg_interval"]/3600.)
            self.ecofac = eco_vars["ecofac"]
            if self.ecofac is None:
                if self.eco_plot:
                    msg = ("Specified quantity to plot is ecological, but model ECOFAC is None. "
                           "Was vegetation turned on in the simulation?")
                    raise ValueError(msg)
                self.ecofac = 1
                self.days_per_year = 365.
            else:
                self.days_per_year = (self.ecofac * self.veg_int_hr * self.n_ets) / 24.

            self.save_freq = eco_vars["save_frequency"]

            # TODO: finalize this (remove if statement) once fully transitioned
            if "save_mortality" in eco_vars:
                self.save_mort = eco_vars["save_mortality"]
            else:
                self.save_mort = True
        elif self.eco_plot:
            raise FileNotFoundError(
                f"The _eco_time_vars.json was not found in the veg_output directory, but an ecological "
                 "variable is being plotted. Make sure vegetation was turned on in the simulation, or "
                 "choose a hydrodynamic quantity to plot."
            )


    def setup_model_loaders(self):
        args = (self.modeldir, self.model_name, self.full_quantity_name, self.eco_plot)
        if self.model_type == "DFM": 
            self.map_loader = DFMMapLoader(*args)
        elif self.model_type == "ANUGA": 
            self.map_loader = ANUGAMapLoader(*args)


    def check_inputs(self):
        self.check_quantity_input()
        self.check_time_inputs()
        

    def check_quantity_input(self):
        if self.plot_vectors and self.quantity != "Velocity":
            raise ValueError("plot_vectors option is only supported for quantity = 'Velocity'.")
        supported_quantities = [
            "Bathymetry",
            "WSE",
            "Depth",
            "Velocity",
            "Max Shear Stress",
            "Stem Height",
            "Stem Diameter",
            "Stem Density",
            "Fractions",
            "Potential Mortality -- Flooding",
            "Potential Mortality -- Desiccation",
            "Potential Mortality -- Uprooting",
            "Potential Mortality -- Burial",
            "Potential Mortality -- Scour",
            "Mortality -- Flooding",
            "Mortality -- Desiccation",
            "Mortality -- Uprooting",
            "Mortality -- Burial",
            "Mortality -- Scour",
            "Mortality -- Total",
        ]
        if self.full_quantity_name not in supported_quantities:
            raise ValueError("Input quantity not supported. Please use a quantity from the "
                                f"following list: {supported_quantities}")
        if self.quantity == "Mortality" and not self.save_mort:
            raise ValueError("Mortality arrays were not saved for this simulation. Please "
                             "specify a different quantity or rerun the simulation with "
                             "mortality output turned on (see run_simulation).")    


    def check_time_inputs(self):
        # veg_plot or not, need to check times against length of actual hydrodynamic simulation
        final_ind = int(math.ceil(self.plot_times["plotHR_f"] / self.plot_times["mapHR_int"]))
        self.map_loader.check_final_index(final_ind)


    def run(self):
        from tqdm import tqdm

        ti_0, ti_f, dti = self.get_loop_indices()
        for i in tqdm(range(ti_0, ti_f, dti)):
            self.create_timestrings(i)
            hydro_i, ets, eco_year = self.get_map_indices(i)
            map_vars = self.map_loader.load(hydro_i, ets, eco_year)
            z_grid, grid2plot = self.get_quantity_grids(map_vars)
            self.plot_quantity(z_grid, grid2plot)

        if self.animate: 
            self.create_gif()


    def create_timestrings(self, i):
        time_parts = self.get_time_breakdown(i)
        fname_timestring = self.format_fname_time(time_parts)
        title_timestring = self.format_title_time(time_parts)
        self.timestrings.append((fname_timestring, title_timestring))
        
        
    def get_time_breakdown(self, i):
        """ Returns a dict time components for formatting titles and filenames """
        sim_hours = float(i*self.veg_int_hr if self.eco_plot else \
                          i*self.plot_times["mapHR_int"])

        time_parts = {
            "sim_days": int(sim_hours // 24),
            "sim_hrs_rem": int(sim_hours % 24),
            "sim_mins_rem": int(sim_hours*60 % 60),  # for if we have plot_int < 1 hr, to avoid duplicate file names
        }

        if self.eco_plot:
            veg_years = sim_hours*self.ecofac/24/self.days_per_year
            time_parts["veg_years"] = int(veg_years)
            time_parts["veg_days_rem"] = int(round((veg_years % 1)*self.days_per_year))  # days after full years
            time_parts["veg_days_tot"] = int(round(sim_hours*self.ecofac/24.))

        return time_parts
    

    def format_fname_time(self, time_parts):
        # Always use hydrodynamic time for the filenames
        sim_days = time_parts["sim_days"]
        sim_hrs_rem = time_parts["sim_hrs_rem"]
        sim_mins_rem = time_parts["sim_mins_rem"]
        if self.plot_times["plotHR_int"] >= 1:
            return f"{sim_days}days_{sim_hrs_rem}hrs"
        else:
            return f"{sim_days}days_{sim_hrs_rem}hrs_{sim_mins_rem}mins"


    def format_title_time(self, time_parts):
        sim_days = time_parts["sim_days"]
        sim_hrs_rem = time_parts["sim_hrs_rem"]
        sim_mins_rem = time_parts["sim_mins_rem"]
        if self.eco_plot and self.plot_label_time == "eco-morphodynamic":
            veg_years = time_parts["veg_years"]
            veg_days_rem = time_parts["veg_days_rem"]
            veg_days_tot = time_parts["veg_days_tot"]
            if veg_years > 0:
                return f"{veg_years} years, {veg_days_rem} days"
            else:
                return f"{veg_days_tot} days"
        elif self.plot_times["plotHR_int"] >= 1:
            return f"{sim_days} days, {sim_hrs_rem} hrs"
        else:
            return f"{sim_days} days, {sim_hrs_rem} hrs, {sim_mins_rem} mins"    
    

    def get_loop_indices(self):
        if self.eco_plot:
            return self.get_eco_loop_indices()
        else:
            return self.get_hydro_loop_indices()
        

    def get_eco_loop_indices(self):
        plotHR_0 = max(self.plot_times["plotHR_0"], 
                       self.veg_int_hr * self.save_freq)

        ti_0 = int(plotHR_0 / self.veg_int_hr)
        ti_f = int(self.plot_times["plotHR_f"] / self.veg_int_hr) + 1
        dti = self.save_freq

        return ti_0, ti_f, dti
    

    def get_hydro_loop_indices(self):
        ti_0 = int(self.plot_times["plotHR_0"] / self.plot_times["mapHR_int"])
        ti_f = int(self.plot_times["plotHR_f"] / self.plot_times["mapHR_int"]) + 1
        dti = int(self.plot_times["plotHR_int"] / self.plot_times["mapHR_int"])

        return ti_0, ti_f, dti


    def get_map_indices(self, i):
        # How we loop over the model time slices will depend on whether we are plotting hydro or veg quantities
        # First need to determine the year, ets, and hydrodynamic model index hydro_ind based on the value of i
        # hydro_ind is simple counter for hydro quantities, but must be multiplied by veg_interval for veg quantities
        hydro_ind, ets, eco_year = i, None, None

        if self.eco_plot:
            ets, eco_year = self.get_eco_times(i)
            hydro_ind = i*self.veg_int_hr
        
        return hydro_ind, ets, eco_year


    def get_eco_times(self, i):
        ets = ((i-1) % self.n_ets) + 1
        eco_year = ((i-1) // self.n_ets) + 1
        return ets, eco_year


    def get_model_type(self):
        # Determine from files in model directory whether ththe model is DFM or ANUGA
        return "DFM" if (self.simdir/"dflowfm").exists() else "ANUGA"


    def get_modeldir(self):
        return self.simdir/"dflowfm" if self.model_type=="DFM" else self.simdir


    def get_model_name(self):
        # Find file with the model file extension, this could be ".mdu" or "_orig.mdu" for DFM, ".mdf" for D3D4, ".sww" for ANUGA
        if self.model_type == "DFM":
            model_files = list(self.modeldir.glob(f"*.mdu"))
        elif self.model_type == "ANUGA":
            model_files = list(self.modeldir.glob(f"*.sww"))
                
        if len(model_files) == 0:
            raise FileNotFoundError(
                f"No {self.model_type} model file found in {self.modeldir}"
            )
        elif len(model_files) > 1:
            raise ValueError(
                f"Multiple {self.model_type} model files found in {self.modeldir}. "
                f"Please specify which model to use or remove extra files."
        )
        return model_files[0].stem


    def get_quantity_grids(self, map_vars):
        # Create interpolation function the first time -> quick interpolations for all other time steps
        if self.plot_method == "interp" and self.interp_func is None:
            self.interp_func = self.create_interp_func(map_vars)
        # Get grid plotting scheme if plotting a regular grid exactly, no interpolation
        if self.plot_method == "exact" and self.empty_grid is None:
            self.empty_grid, self.x_idx, self.y_idx = self.get_exact_grid_idx(map_vars)

        # For DFM, z_grid should be recomputed every step in case morphology is turned on. For ANUGA, oh well it's fast enough
        z_grid = self.create_grid(map_vars["Bathymetry"]) * self.quantity_units["Bathymetry"][1]

        if self.quantity == "Bathymetry":
            return z_grid, None
        elif self.eco_plot:
            eco_grid = self.get_eco_quantity_grids(map_vars)
            if eco_grid is None:
                # Return empty, not None, because we want a blank plot but also want to retain colorbar properties
                eco_grid = np.ma.masked_all(z_grid.shape)
            return z_grid, eco_grid
        else:
            return z_grid, self.get_hydro_quantity_grids(map_vars)


    def get_eco_quantity_grids(self, map_vars):
        # Adjust "mortality" key if necessary, from full name to generic name
        if self.quantity == "Mortality":
            map_vars["Mortality"] = map_vars.pop(self.full_quantity_name)

        if len(map_vars[self.quantity]) > 0:
            if self.quantity in ["Fractions", "Mortality"]:
                frac_grid_list = []
                # --- Plot fractional quantities by species (summation) --- #
                if self.plot_separate_species:
                    species_list = list(set(spec[0] for spec in map_vars["Cohort Names"]))  # exhaustive set of species names
                    for species in species_list:
                        fractions = [frac for frac, c in zip(map_vars[self.quantity], map_vars["Cohort Names"]) if c[0]==species]
                        fractions_sum = np.sum(fractions, axis=0)
                        eco_grid = self.create_grid(fractions_sum) * self.quantity_units[self.quantity][1]
                        frac_grid_list.append((np.ma.masked_where(eco_grid < self.quantity_units[self.quantity][2], eco_grid),
                                               species))
                # --- Plot fractional quantities by cohort --- #
                else:
                    for frac, cohort in zip(map_vars[self.quantity], map_vars["Cohort Names"]):
                        eco_grid = self.create_grid(frac) * self.quantity_units[self.quantity][1]
                        frac_grid_list.append((np.ma.masked_where(eco_grid < self.quantity_units[self.quantity][2], eco_grid),
                                               cohort))
                return frac_grid_list
            # --- Separate quantities (stem height, diameter, etc) by species --- #
            elif self.plot_separate_species:
                eco_grid_list = []
                species_list = list(set(spec[0] for spec in map_vars["Cohort Names"]))  # exhaustive set of species names
                for species in species_list:
                    fractions = [frac for frac, c in zip(map_vars["Fractions"], map_vars["Cohort Names"]) if c[0]==species]
                    quantities = [q for q, c in zip(map_vars[self.quantity], map_vars["Cohort Names"]) if c[0]==species]
                    eco_grid = self.compute_stem_quantity_grid(fractions, quantities)
                    eco_grid_list.append((np.ma.masked_where(eco_grid < self.quantity_units[self.quantity][2], eco_grid),
                                          species))
                return eco_grid_list
            # --- Combine quantities (stem height, diameter, etc), single value per grid cell --- #
            else:
                eco_grid = self.compute_stem_quantity_grid(map_vars["Fractions"], map_vars[self.quantity])
                return np.ma.masked_where(eco_grid < self.quantity_units[self.quantity][2], eco_grid)
        else:
            return None  # if no veg data at all in this time step


    def get_hydro_quantity_grids(self, map_vars):
        depth_grid = self.create_grid(map_vars["Depth"]) * self.quantity_units[self.quantity][1]            
        show_inds = depth_grid < self.quantity_units["Depth"][2]
        if self.quantity == "Depth":
            return np.ma.masked_where(show_inds, depth_grid)
        elif self.quantity == "WSE":
            wse_grid = self.create_grid(map_vars[self.quantity]) * self.quantity_units[self.quantity][1]
            return np.ma.masked_where(show_inds, wse_grid)
        elif self.quantity == "Velocity":
            V_grid = self.create_grid(map_vars[self.quantity]) * self.quantity_units[self.quantity][1]
            if self.plot_vectors:
                Vx_grid = self.create_grid(map_vars["Vel_x"]) * self.quantity_units[self.quantity][1]
                Vy_grid = self.create_grid(map_vars["Vel_y"]) * self.quantity_units[self.quantity][1]
                self.vector_comps = (np.ma.masked_where(show_inds, Vx_grid),
                                     np.ma.masked_where(show_inds, Vy_grid),
                                     np.ma.masked_where(show_inds, V_grid))
            return np.ma.masked_where(show_inds, V_grid)
        elif self.quantity == "Max Shear Stress":
            tau_grid = self.create_grid(map_vars[self.quantity]) * self.quantity_units[self.quantity][1]
            return np.ma.masked_where(tau_grid < self.quantity_units["Max Shear Stress"][2], tau_grid)
        else:
            raise ValueError(f"Unexpected quantity name: {self.quantity}")
        

    def create_interp_func(self, map_vars):
        # Create interpolation function the first time -> quick interpolations for all other time steps
        interp_func = create_nn_interpFunc(map_vars["X"], map_vars["Y"], 
                                           grid_size=self.plot_interp_vars["cell_size"], 
                                           k_nn=self.plot_interp_vars["n_neighbors"],
                                           polygon_csv=self.mask_bndy_file, extents=self.extents)
        return interp_func


    def get_exact_grid_idx(self, map_vars):
        # Get number of unique x/y values 
        x_unique = np.unique(map_vars["X"])
        y_unique = np.unique(map_vars["Y"])
        # Get number of rows and columns
        nx = x_unique.size
        ny = y_unique.size
        # Create empty template grid
        empty_grid = np.full((ny, nx), np.nan)
        # Map 1-D points to grid indices
        x_idx = np.searchsorted(x_unique, map_vars["X"])
        y_idx = np.searchsorted(y_unique, map_vars["Y"])

        return empty_grid, x_idx, y_idx
    
    
    def create_grid(self, var):
        # if plot_method is "interp"
        if self.interp_func is not None:
            grid = self.interp_func(var)
        # if plot_method is "exact"
        elif self.empty_grid is not None:
            grid = self.empty_grid.copy()
            grid[self.y_idx, self.x_idx] = var
        return grid


    def compute_stem_quantity_grid(self, frac, q):
        if self.quantity == "Stem Density":  # sum_product for "Stem Density"
            veg_data = sum_product(frac, q)
        else:  # cell_averaging for "Stem Diameter", "Stem Face Factor", and "Stem Height"
            veg_data = cell_averaging(frac, q)
        return self.create_grid(veg_data) * self.quantity_units[self.quantity][1]
    

    def plot_quantity(self, base_grid, main_grid):
        # For veg fractions and mortality, or separate species plotting
        if type(main_grid) is list:
            for cohort in main_grid:
                grid = cohort[0]
                c_id = cohort[1]  # species name, or (species name, cohort ID)
                if self.plot_separate_species:
                    c_str = c_id
                    c_fstr = c_id
                else:
                    c_str = f"{c_id[0]} (cohort {c_id[1]+1})"
                    c_fstr = f"{c_id[0]}{c_id[1]}"                   
                self.plot_single_quantity(
                    base_grid, 
                    grid, 
                    # E.g.: title = "Stem Height -- 1 year, 50 days \n SpartinaAnglica"
                    title=f"{self.full_quantity_name} -- {self.timestrings[-1][1]}\n{c_str}",
                    fname=f"{self.full_quantity_name.replace(' ', '')}_{c_fstr}_{self.timestrings[-1][0]}"
                    )
        # Hydrodynamic plots or averaged vegetation quantity plots
        elif self.quantity == "Bathymetry" or not main_grid.mask.all():  # don't plot empty vegetation plots
            self.plot_single_quantity(
                base_grid, 
                main_grid, 
                title=f"{self.full_quantity_name} -- {self.timestrings[-1][1]}",
                fname=f"{self.full_quantity_name.replace(' ', '')}_{self.timestrings[-1][0]}"
                )
            
            
    def plot_single_quantity(self, base_grid, main_grid, title, fname):

        fig, ax = plt.subplots(figsize=self.plot_specs["figsize"])

        base = self.imshow(ax, base_grid,
                           cmap=self.cmaps["Bathymetry"],
                           vlims=self.cmap_lims["Bathymetry"],
                           )
        
        if main_grid is not None:
            main = self.imshow(ax, main_grid,
                               cmap=self.cmaps[self.quantity],
                               vlims=self.cmap_lims[self.quantity],
                               )

            if self.plot_vectors:
                self.make_quiver(ax, main_grid)

        if self.show_title:
            ax.set_title(title, fontsize=self.plot_specs["fontsize"])
        if not self.show_axis_ticks:
            ax.set_xticks([])
            ax.set_yticks([])

        cbar = plt.colorbar(main if main_grid is not None else base, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(f"{self.quantity_units[self.quantity][0]}", rotation=270, labelpad=25, fontsize=self.plot_specs["fontsize"])
        cbar.ax.tick_params(labelsize=self.plot_specs["fontsize"])

        if self.show_topo_cbar:
            cbar_z = plt.colorbar(base, ax=ax, location='left', fraction=0.046, pad=0.04)
            cbar_z.set_label("Bathymetry", rotation=90, labelpad=5)   
            cbar_z.ax.tick_params(labelsize=self.plot_specs["fontsize"]) 

        if self.scalebar:
            self.make_scalebar(ax)

        fig.tight_layout()

        self.output_plot_dir.mkdir(parents=True, exist_ok=True)
        fig_path = self.output_plot_dir / fname
        plt.savefig(f"{fig_path}.png", dpi=self.plot_specs["output_dpi"], bbox_inches="tight")
        self.img_paths.append(f"{fig_path}.png")
        if self.save_grids:
            np.savez_compressed(f"{fig_path}.npz", data=main_grid.data, mask=main_grid.mask)
        plt.close()


    def imshow(self, ax, grid, cmap, vlims):
        extent = None
        if self.show_axis_ticks:
            gridX, gridY = self.get_grid_coords(grid)
            extent = [gridX.min(), gridX.max(), gridY.max(), gridY.min()]

        return ax.imshow(grid,
                         cmap=cmap,
                         vmin=vlims[0],
                         vmax=vlims[1],
                         extent=extent,
                         origin="upper",
                         aspect=self.plot_specs["aspect"],
                         )
    

    def make_quiver(self, ax, grid):  # See pyplot.quiver manual for more info
        # Get x and y components of velocity
        vx, vy, V = self.vector_comps  # created only when plot_vectors and quantity is "Velocity"
        # Normalize to get unit vectors
        with np.errstate(divide="ignore", invalid="ignore"):
            vx, vy = vx/V, vy/V

        gridX, gridY = self.get_grid_coords(grid)

        # Identify mesh array indices where we want to show vectors (cut down number of arrows)
        cell_size = self.plot_interp_vars["cell_size"]
        if type(self.vector_props["vect_spacing"]) in [tuple, list]:
            dx = int(self.vector_props["vect_spacing"][0] / cell_size)
            dy = int(self.vector_props["vect_spacing"][1] / cell_size)
        else:
            dx = dy = int(self.vector_props["vect_spacing"] / cell_size)
        
        ax.quiver(gridX[::dy, ::dx], gridY[::dy, ::dx], 
                  vx[::dy, ::dx], vy[::dy, ::dx],
                  color=self.vector_props["color"],
                  scale=self.vector_props["scale"],
                  pivot=self.vector_props["pivot"],
                  width=self.vector_props["width"],
                  )

    
    def get_grid_coords(self, grid):
        """ Return meshgrid coordinates in the same space as imshow extent """
        ny, nx = grid.shape

        if self.show_axis_ticks:
            xx = np.arange(nx) * self.plot_interp_vars["cell_size"]
            yy = np.arange(ny) * self.plot_interp_vars["cell_size"]
        else:
            xx = np.arange(nx)
            yy = np.arange(ny)

        return np.meshgrid(xx, yy)
        

    def make_scalebar(self, ax):
        unit_conv = 1.
        bar_text = self.scalebar_props["distance"]
        if bar_text.split()[1] == "km": unit_conv = 1000.
        elif bar_text.split()[1] == "mi": unit_conv = 1609.
        scalebar = AnchoredSizeBar(
            ax.transData,
            float(bar_text.split()[0])/self.plot_interp_vars["cell_size"]*unit_conv, 
            bar_text, 
            self.scalebar_props["loc"],
            pad=0.5,
            color="k",
            frameon=self.scalebar_props["frameon"],
            sep=4,
            size_vertical=2,
        )
        ax.add_artist(scalebar)
    

    def create_gif(self):
        import imageio.v2 as imageio
        print("Creating animation...")
        fps = 5 if self.eco_plot else 10
        images = [imageio.imread(p) for p in self.img_paths]
        gif_path = self.output_plot_dir / "animation.gif"
        try:
            imageio.mimsave(str(gif_path), images, fps=fps, loop=0)  # type: ignore[arg-type]
        except ValueError as e:
            if "all input arrays must have the same shape" in str(e):
                raise ValueError(
                    "Animation failed to be created because frames have inconsistent dimensions. "
                    "This is often caused by one of more plot elements (titles, axis labels, "
                    "colorbar elements, etc.) being too long relative to one of the figure "
                    "dimensions. One fix could be to change the 'figsize' parameter by providing "
                    "a non-default argument for 'plot_specs'. Otherwise, inspect the figures that "
                    "were created to see why some of them may be different sizes."
                ) from e
            raise  # re-raise unrelated ValueErrors unchanged
        if self.delete_static_imgs:
            for img in self.img_paths:
                Path.unlink(img)
