"""
Class for plotting various output quantites from DYCOVE-ANUGA/DFM
"""

from typing import Optional, Any, Union
import numpy as np
import math
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from matplotlib import path as mpath
from tqdm import tqdm
from scipy.spatial import cKDTree  # type: ignore

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

    #if tree_inds is None and inside_mask is None:
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
            with np.errstate(divide='ignore'):
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

    hydro_plot_label_time : str, optional
        Time scale to use for plot labels for hydro quantities (e.g., depth, velocity). 
        Default is 'hydrodynamic'.
    veg_plot_label_time : str, optional
        Time scale to use for plot labels for veg quantities (e.g., stem height, fractions).
        Default is 'eco-morphodynamic'.
    plot_label_time : str, optional
        Time scale to use in plot titles, either ``'eco-morphodynamic'`` (default) or
        ``'hydrodynamic'``.
    n_ets : int, optional
        Number of ecological time steps (ETS) per ecological year (default 14).
    veg_interval : int, optional
        Hydrodynamic time per ETS, in seconds (default 43200).
    ecofac : int, optional
        Multiplicative factor converting hydro to eco time. If None, computed from
        ``n_ets`` and ``veg_interval``.
    plot_method : str, optional
        ``'interp'`` (default) for nearest-neighbor interpolation to a regular grid, or
        ``'exact'`` for plotting on the original mesh nodes/cells (must be regular grid).
    plot_interp_vars : dict, optional
        Options for interpolating to 2D grid, e.g. ``{'cell_size': 5, 'n_neighbors': 3}``.
    plot_specs : dict, optional
        Plot appearance options, e.g., 
        ``{'figsize': (6, 6), 'fontsize': 12, 'output_dpi': 150}``.
    cmap_lims : dict, optional
        Colorbar limits for each quantity (min, max), e.g. ``{'Bathymetry': (0, 1)}``.
    cmaps : dict, optional
        Colormaps for each quantity. Accepts Matplotlib ``Colormap`` objects.
    quantity_units : dict, optional
        Unit labels, scaling factors, and minimum plotting thresholds for each quantity, 
        e.g. ``{'Depth': ('[m]', 1, 0.01), 'Stem Diameter': ('[mm]', 1000, 1)}``.
    plot_vectors : bool, optional
        If True, plots normalized velocity vectors (quantity=="Velocity" only).
    vector_props : dict, optional
        Scaling and spacing options for velocity vectors, used if quantity == "Velocity",
        e.g., ``{'vect_spacing': 50, 'color': 'red', 'scale': 30, 'pivot': 'mid', 'width': 0.003}``.
    show_title : bool, optional
        If True, include auto-generated plot title at top of plot.
    show_topo_cbar : bool, optional
        If True, include elevation color bar on the left-hand side.
    scalebar : bool, optional
        If True, include a scalebar object on the plot.
    scalebar_props : dict, optional
        Scalebar properties, e.g., ``{'distance': '200 m', 'loc': 'upper left', 'frameon': True}``.
    mask_bndy_file : str or Path, optional
        Path to polygon CSV file for masking interpolation outside the domain. Probably
        required if model domain is not rectangular.
    extents : tuple[float, float, float, float], optional
        Geographic extents for interpolation (xmin, xmax, ymin, ymax). Useful for
        zooming in on a particular location.
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
                 hydro_plot_label_time='hydrodynamic',
                 veg_plot_label_time='eco-morphodynamic',
                 n_ets = 14,
                 veg_interval = 43200,
                 ecofac = None,
                 plot_method = 'interp',
                 plot_interp_vars: Optional[dict[str, int]] = None,
                 plot_specs: Optional[dict[str, Any]] = None, 
                 cmap_lims: Optional[dict[str, tuple[float, float]]] = None, 
                 cmaps: Optional[dict[str, Colormap]] = None,
                 quantity_units: Optional[dict[str, tuple[str, Union[int, float]]]] = None, 
                 plot_vectors = False,
                 vector_props: Optional[dict[str, Any]] = None,
                 show_title = True,
                 show_topo_cbar = False,
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
        self.plot_method = plot_method
        self.plot_vectors = plot_vectors
        self.show_title = show_title
        self.show_topo_cbar = show_topo_cbar
        self.scalebar = scalebar
        self.mask_bndy_file = mask_bndy_file
        self.extents = extents
        self.animate = animate
        self.delete_static_imgs = delete_static_imgs
        self.save_grids = save_grids

        if self.plot_vectors and self.quantity != "Velocity":
            raise ValueError("plot_vectors option is only supported for quantity = 'Velocity'.")
        
        # Simplify quantity name if "Mortality"
        self.full_quantity_name = self.quantity
        if "Mortality" in self.full_quantity_name.split():
            self.quantity = "Mortality"

        # Ecological time inputs
        self.n_ets = n_ets
        self.veg_int_hr = int(veg_interval/3600.)
        self.ecofac, self.days_per_year = self.compute_ecofac(ecofac)

        # If plot_method is 'interp', will build interpolation function for converting 1-D arrays to a grid
        self.interp_func = None
        # If plot_method is 'exact', will build empty grid and x/y indices for mapping 1-D arrays to grid
        self.empty_grid, self.x_idx, self.y_idx = None, None, None
            
        default_plot_interp_vars = {
            'cell_size': 5,    # in meters
            'n_neighbors': 1,  # typically 1 or 3
        }
        self.plot_interp_vars = {**default_plot_interp_vars, **(plot_interp_vars or {})}  # merge provided custom values with default values

        default_plot_specs = {
            'figsize': (6, 6),
            'fontsize': 12,
            'output_dpi': 100,
        }
        self.plot_specs = {**default_plot_specs, **(plot_specs or {})}  # merge provided custom values with default values

        # Colorbar plotting limits, check quantity_units dictionary for unit consistency
        default_cmap_lims = {
            'Bathymetry'      : (0, 1),
            'WSE'             : (0, 1),
            'Depth'           : (0, 3),
            'Velocity'        : (0, 0.5),
            'Max Shear Stress': (0, 0.1),
            'Fractions'       : (0, 1),
            'Stem Height'     : (0, 2),
            'Stem Density'    : (0, 300),
            'Stem Diameter'   : (0, 50),
            'Mortality'       : (0, 100),
        }
        self.cmap_lims = {**default_cmap_lims, **(cmap_lims or {})}  # merge provided custom values with default values

        default_cmaps = {
            'Bathymetry'      : ListedColormap(plt.colormaps["Greys_r"](np.linspace(0, 1.0, 256))),
            #'Bathymetry'      : ListedColormap(plt.colormaps["bone"](np.linspace(0, 1.0, 256))),
            'WSE'             : ListedColormap(plt.colormaps["Blues_r"](np.linspace(0.1, 1.0, 256))),
            'Depth'           : ListedColormap(plt.colormaps["Blues"](np.linspace(0.1, 1.0, 256))),
            'Velocity'        : ListedColormap(plt.colormaps["viridis"](np.linspace(0, 1.0, 256))),
            'Max Shear Stress': ListedColormap(plt.colormaps["plasma"](np.linspace(0, 1.0, 256))),
            'Fractions'       : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Height'     : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Density'    : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Diameter'   : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Mortality' : ListedColormap(plt.colormaps["Oranges"](np.linspace(0.2, 1.0, 256))),
        }
        self.cmaps = {**default_cmaps, **(cmaps or {})}  # merge provided custom values with default values

        # First value: units to use for plot title/colorbar
        # Second value: multiplier based on most standard units (meters, seconds, Newtons, fractions, etc)
        # Third value: minimum value for plotting, anything below this value is masked out
        # examples given for some common conversions
        # Note that the third value for Depth is used to mask out other quantities, e.g. Velocity and WSE
        # check cmap_lims dictionary to ensure units are consistent with colorbar plotting limits
        default_quantity_units = {
            'Bathymetry'      : ('[m]', 1, -99),  # use ('[ft]', 3.281, -99) for elevations in feet
            'WSE'             : ('[m]', 1, -99),
            'Depth'           : ('[m]', 1, 0.01),
            'Velocity'        : ('[m/s]', 1, 0),  # use ('[cm/s]', 100, 0) for velocities in cm/s
            'Max Shear Stress': ('[N/m$^2$]', 1, 0.001),
            'Fractions'       : ('[--]', 1, 0.01),
            'Stem Height'     : ('[m]', 1, 0.1),
            'Stem Density'    : ('[stems/m$^2$]', 1, 5),
            'Stem Diameter'   : ('[mm]', 1000, 1),  # use ('[m]', 1, 0.001) for diameter in meters
            'Mortality'       : ('[%]', 100, 1),  # use ('[--]', 1, 0.01) for fraction between 0 and 1
        }
        self.quantity_units = {**default_quantity_units, **(quantity_units or {})}  # merge provided custom values with default values

        default_vector_props = {
            'vect_spacing': 50,
            'color': 'red',
            'scale': 30,
            'pivot': 'mid',
            'width': 0.003,
        }
        self.vector_props = {**default_vector_props, **(vector_props or {})}  # merge provided custom values with default values

        default_scalebar_props = {
            'distance': '200 m',
            'loc': 'upper left',
            'frameon': True,
        }
        self.scalebar_props = {**default_scalebar_props, **(scalebar_props or {})}  # merge provided custom values with default values

        self.output_plot_dir = self.simdir / 'figures' / self.full_quantity_name.replace(' ', '')
        self.model_type = self.get_model_type()
        self.modeldir = self.get_modeldir()
        self.model_name = self.get_model_name()
        self.eco_plot = self.quantity in ['Fractions', 'Stem Height', 'Stem Density', 'Stem Diameter', 'Mortality']

        # Change how we label time in our plots depending on if we are plotting vegetation (uses ecofac) or hydrodynamics (no morfac/ecofac)
        self.plot_label_time = veg_plot_label_time if self.eco_plot else hydro_plot_label_time

        # Load model output files using the appropriate loader class
        args = (self.modeldir, self.model_name, self.full_quantity_name, self.eco_plot, self.n_ets)
        if self.model_type == 'DFM': 
            self.map_loader = DFMMapLoader(*args)
        elif self.model_type == 'ANUGA': 
            self.map_loader = ANUGAMapLoader(*args)

        # To save for building animations later
        self.timestrings = []
        self.img_paths = []


    def run(self):
        ti_0, ti_f, dti = self.get_time_indices()
        self.check_time_inputs()
        for i in tqdm(range(ti_0, ti_f, dti)):
            self.create_timestrings(i)
            hydro_i, ets, eco_year = self.get_map_indices(i)
            map_vars, species_names = self.map_loader.load(hydro_i, ets, eco_year)
            z_grid, grid2plot = self.get_quantity_grids(map_vars)
            self.plot_quantity(z_grid, grid2plot, species_names)

        if self.animate: 
            self.create_gif()


    def compute_ecofac(self, ecofac_input):
        if ecofac_input is None:
            # Compute ideal (continuous) ecofac
            ecofac_est = (24 * 365.) / (self.veg_int_hr * self.n_ets)
            # Round to nearest whole number
            ecofac = int(round(ecofac_est))
            msg = f"Computed ecofac = {ecofac:d}, derived from veg_interval and n_ets (rounded from {ecofac_est:.3f})"
        else:
            ecofac = ecofac_input
            msg = f"Using provided ecofac = {ecofac:d}"

        # Compute associated days_per_year after ecofac rounding
        days_per_year = (ecofac * self.veg_int_hr * self.n_ets) / 24.

        r.report(f"Computing simulation times based on n_ets = {self.n_ets} and veg_interval = {self.veg_int_hr*3600}")
        r.report(msg)
        r.report(f"ecofac = {ecofac:d} corresponds to {days_per_year:.1f} days per year.")

        return ecofac, days_per_year
        

    def create_timestrings(self, i):
        time_parts = self.get_time_breakdown(i)
        fname_timestring = self.format_fname_time(time_parts)
        title_timestring = self.format_title_time(time_parts)
        self.timestrings.append((fname_timestring, title_timestring))


    def get_time_breakdown(self, i):
        """Returns a dict with time components for formatting titles and filenames."""
        sim_hours = float(i*self.plotHR_int if self.eco_plot else i*self.plotHR_div)
        
        veg_years = sim_hours*self.ecofac/24/self.days_per_year
        veg_days_rem = (veg_years % 1)*self.days_per_year  # days after full years
        veg_days_tot = sim_hours*self.ecofac/24.

        return {"sim_days": int(sim_hours // 24),
                "sim_hrs_rem": int(sim_hours % 24),
                "sim_mins_rem": int(sim_hours*60 % 60),  # for if we have a plot interval less than one hour, to avoid duplicate file names
                "veg_years": int(veg_years),
                "veg_days_rem": int(round(veg_days_rem)),
                "veg_days_tot": int(round(veg_days_tot))}
    

    def format_fname_time(self, time_parts):
        # Always use hydrodynamic time for the filenames
        if self.plotHR_int >= 1:
            return f"{time_parts['sim_days']}days_{time_parts['sim_hrs_rem']}hrs"
        else:
            return f"{time_parts['sim_days']}days_{time_parts['sim_hrs_rem']}hrs_{time_parts['sim_mins_rem']}mins"


    def format_title_time(self, time_parts):
        if self.plot_label_time == 'eco-morphodynamic':
            if time_parts["veg_years"] > 0:
                return f"{time_parts['veg_years']} years, {time_parts['veg_days_rem']} days"
            else:
                return f"{time_parts['veg_days_tot']} days"
        elif self.plotHR_int >= 1:
            return f"{time_parts['sim_days']} days, {time_parts['sim_hrs_rem']} hrs"
        else:
            return f"{time_parts['sim_days']} days, {time_parts['sim_hrs_rem']} hrs, {time_parts['sim_mins_rem']} mins"    
    

    def get_time_indices(self):
        if self.eco_plot:
            self.plotHR_0 = max(self.plot_times['plotHR_0'], self.veg_int_hr)
            self.plotHR_int = self.veg_int_hr
            self.plotHR_div = self.veg_int_hr
        else:
            self.plotHR_0 = self.plot_times['plotHR_0']
            self.plotHR_int = self.plot_times['plotHR_int']
            self.plotHR_div = self.plot_times['mapHR_int']

        ti_0 = int(self.plotHR_0/self.plotHR_div)
        # Adding 1 to include the final time step. TODO: check this
        ti_f = int(self.plot_times['plotHR_f']/self.plotHR_div) + 1
        dti = int(self.plotHR_int/self.plotHR_div)

        return ti_0, ti_f, dti


    def check_time_inputs(self):
        final_ind = int(math.ceil(self.plot_times['plotHR_f'] / self.plot_times['mapHR_int']))
        self.map_loader.check_final_index(final_ind)


    def get_map_indices(self, i):
        # How we loop over the model time slices will depend on whether we are plotting hydro or veg quantities
        # First need to determine the year, ets, and hydro_ind based on the value of i
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

    def get_hydro_eco_times(self, i):
        ets = ((i // self.veg_int_hr) % self.n_ets) + 1
        eco_year = (i // (self.veg_int_hr * self.n_ets)) + 1
        return ets, eco_year

    def get_model_type(self):
        # Determine from files in model directory whether ththe model is DFM or ANUGA
        return 'DFM' if (self.simdir/'dflowfm').exists() else 'ANUGA'

    def get_modeldir(self):
        return self.simdir/'dflowfm' if self.model_type=='DFM' else self.simdir

    def get_model_name(self):
        # Find file with the model file extension, this could be ".mdu" or "_orig.mdu" for DFM, ".mdf" for D3D4, ".sww" for ANUGA
        if self.model_type == 'DFM':
            model_files = list(self.modeldir.glob(f"*.mdu"))
        elif self.model_type == 'ANUGA':
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
        # Adjust "mortality" key if necessary, from full name to generic name
        if self.quantity == "Mortality":
            map_vars["Mortality"] = map_vars.pop(self.full_quantity_name)

        # Create interpolation function the first time -> quick interpolations for all other time steps
        if self.plot_method == 'interp' and self.interp_func is None:
            self.interp_func = self.create_interp_func(map_vars)

        if self.plot_method == 'exact' and self.empty_grid is None:
            self.empty_grid, self.x_idx, self.y_idx = self.get_exact_grid_idx(map_vars)
            
        # to be populated if we are doing vector plots
        self.vector_components = None

        # For DFM, z_grid should be recomputed every step in case morphology is turned on. For ANUGA, oh well it's fast enough
        z_grid = self.create_grid(map_vars['Bathymetry']) * self.quantity_units['Bathymetry'][1]

        if self.quantity == 'Bathymetry':
            return z_grid, None
        elif not self.eco_plot:
            depth_grid = self.create_grid(map_vars['Depth']) * self.quantity_units[self.quantity][1]            
            show_inds = depth_grid < self.quantity_units['Depth'][2]
            if self.quantity == 'Depth':
                return z_grid, np.ma.masked_where(show_inds, depth_grid)
            elif self.quantity == 'WSE':
                wse_grid = self.create_grid(map_vars[self.quantity]) * self.quantity_units[self.quantity][1]
                return z_grid, np.ma.masked_where(show_inds, wse_grid)
            elif self.quantity == 'Velocity':
                V_grid = self.create_grid(map_vars[self.quantity]) * self.quantity_units[self.quantity][1]
                if self.plot_vectors:
                    Vx_grid = self.create_grid(map_vars['Vel_x']) * self.quantity_units[self.quantity][1]
                    Vy_grid = self.create_grid(map_vars['Vel_y']) * self.quantity_units[self.quantity][1]
                    self.vector_components = (np.ma.masked_where(show_inds, Vx_grid),
                                              np.ma.masked_where(show_inds, Vy_grid),
                                              np.ma.masked_where(show_inds, V_grid))
                return z_grid, np.ma.masked_where(show_inds, V_grid)
            elif self.quantity == 'Max Shear Stress':
                tau_grid = self.create_grid(map_vars[self.quantity]) * self.quantity_units[self.quantity][1]
                return z_grid, np.ma.masked_where(tau_grid < self.quantity_units['Max Shear Stress'][2], tau_grid)
            else:
                raise ValueError(f"Unexpected quantity name: {self.quantity}")
        else:
            if len(map_vars[self.quantity]) > 0:
                if self.quantity in ['Fractions', 'Mortality']:
                    fraction_grid_list = []
                    for frac in map_vars[self.quantity]:
                        veg_grid = self.create_grid(frac) * self.quantity_units[self.quantity][1]
                        fraction_grid_list.append(np.ma.masked_where(veg_grid < self.quantity_units[self.quantity][2], veg_grid))
                    return z_grid, fraction_grid_list
                elif self.quantity == 'Stem Density':  # sum_product for 'Stem Density'
                    veg_data = sum_product(map_vars['Fractions'], map_vars[self.quantity])
                else:  # cell_averaging for 'Stem Diameter', 'Stem Face Factor', and 'Stem Height'
                    veg_data = cell_averaging(map_vars['Fractions'], map_vars[self.quantity])
                veg_grid = self.create_grid(veg_data) * self.quantity_units[self.quantity][1]
                return z_grid, np.ma.masked_where(veg_grid < self.quantity_units[self.quantity][2], veg_grid)
            else:
                # If no veg data at all in this time step, return an empty grid
                empty_grid = np.ma.masked_all(z_grid.shape)
                return z_grid, empty_grid


    def create_interp_func(self, map_vars):
        # create interpolation function the first time -> quick interpolations for all other time steps
        interp_func = create_nn_interpFunc(map_vars['X'], map_vars['Y'], 
                                           grid_size=self.plot_interp_vars['cell_size'], 
                                           k_nn=self.plot_interp_vars['n_neighbors'],
                                           polygon_csv=self.mask_bndy_file, extents=self.extents,
                                           )
        return interp_func

    def get_exact_grid_idx(self, map_vars):
        # Get number of unique x/y values 
        x_unique = np.unique(map_vars['X'])
        y_unique = np.unique(map_vars['Y'])
        # Get number of rows and columns
        nx = x_unique.size
        ny = y_unique.size
        # Create empty template grid
        empty_grid = np.full((ny, nx), np.nan)
        # Map 1-D points to grid indices
        x_idx = np.searchsorted(x_unique, map_vars['X'])
        y_idx = np.searchsorted(y_unique, map_vars['Y'])

        return empty_grid, x_idx, y_idx
    
    def map_to_grid(self, var_1d):
        grid = self.empty_grid.copy()
        grid[self.y_idx, self.x_idx] = var_1d
        return grid
    
    def create_grid(self, var):
        # if plot_method is 'interp'
        if self.interp_func is not None:
            grid = self.interp_func(var)
        # if plot_method is 'exact'
        elif self.empty_grid is not None:
            grid = self.map_to_grid(var)
        return grid


    def plot_quantity(self, base_grid, main_grid, species_names):
        if type(main_grid) is list:  # for veg fractions, we plot each fraction separately
            for grid, name in zip(main_grid, species_names):
                self.plot_single_quantity(base_grid, grid, 
                        title=f"{self.full_quantity_name} -- {name} -- {self.timestrings[-1][1]}",
                        fname=f"{self.full_quantity_name.replace(' ', '')}_{name}_{self.timestrings[-1][0]}")
        else:
            self.plot_single_quantity(base_grid, main_grid, 
                        title=f"{self.full_quantity_name} -- {self.timestrings[-1][1]}",
                        fname=f"{self.full_quantity_name.replace(' ', '')}_{self.timestrings[-1][0]}")
            
            
    def plot_single_quantity(self, base_grid, main_grid, title, fname):

        fig, ax = plt.subplots(figsize=self.plot_specs['figsize'])

        base = ax.imshow(base_grid,
                         cmap=self.cmaps['Bathymetry'],
                         vmin=self.cmap_lims['Bathymetry'][0],
                         vmax=self.cmap_lims['Bathymetry'][1])
        
        if main_grid is not None:
            main = ax.imshow(main_grid,
                             cmap=self.cmaps[self.quantity],
                             vmin=self.cmap_lims[self.quantity][0],
                             vmax=self.cmap_lims[self.quantity][1])            

        if self.plot_vectors:
            self.do_quiver(ax, base_grid)

        if self.show_title:
            ax.set_title(title, fontsize=self.plot_specs['fontsize'])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.12)
        cbar = plt.colorbar(main if self.quantity != 'Bathymetry' else base, cax=cax)
        cbar.set_label(f"{self.quantity} {self.quantity_units[self.quantity][0]}", rotation=270, labelpad=25, fontsize=self.plot_specs['fontsize'])
        cbar.ax.tick_params(labelsize=self.plot_specs['fontsize'])

        if self.show_topo_cbar:
            cax_z = divider.append_axes("left", size="3%", pad=0.12)
            cbar_z = plt.colorbar(base, cax=cax_z)
            cbar_z.set_label(f"Bathymetry {self.quantity_units['Bathymetry'][0]}", rotation=90, labelpad=5, fontsize=self.plot_specs['fontsize'])
            cbar_z.ax.yaxis.set_label_position("left")
            cbar_z.ax.yaxis.tick_left()
            cbar_z.ax.tick_params(labelsize=self.plot_specs['fontsize'], direction="out", labelrotation=0)          

        if self.scalebar:
            unit_conv = 1.
            bar_text = self.scalebar_props['distance']
            if bar_text.split()[1] == 'km': unit_conv = 1000.
            elif bar_text.split()[1] == 'mi': unit_conv = 1609.
            scalebar = AnchoredSizeBar(ax.transData,
                                    float(bar_text.split()[0])/self.plot_interp_vars['cell_size']*unit_conv, 
                                    bar_text, 
                                    self.scalebar_props['loc'],
                                    pad=0.5,
                                    color='k',
                                    frameon=self.scalebar_props['frameon'],
                                    sep=4,
                                    size_vertical=2,
                        # fontproperties=FontProperties(size=fontsizes['ax_tix'], #weight='bold'
                        #                                             )
                                                                    )
            ax.add_artist(scalebar)

        fig.tight_layout()

        self.output_plot_dir.mkdir(parents=True, exist_ok=True)
        fig_path = self.output_plot_dir / fname
        plt.savefig(f"{fig_path}.png", dpi=self.plot_specs['output_dpi'], bbox_inches='tight')
        self.img_paths.append(f"{fig_path}.png")
        if self.save_grids:
            np.savez_compressed(f"{fig_path}.npz", data=main_grid.data, mask=main_grid.mask)
        plt.close()


    def do_quiver(self, ax, base_grid):  # See pyplot.quiver manual for more info
        # get x and y components of velocity
        vx, vy, V = self.vector_components
        # normalize to get unit vectors
        with np.errstate(divide='ignore', invalid='ignore'):
            vx, vy = vx/V, vy/V
        # get dimenions of grid
        ny, nx = base_grid.shape
        gridX, gridY = np.meshgrid(np.linspace(0, nx, nx), np.linspace(0, ny, ny))
        # Identify mesh array indices where we want to show vectors (cut down number of arrows)
        dr = int(self.vector_props['vect_spacing']/self.plot_interp_vars['cell_size'])
        ax.quiver(gridX[::dr,::dr], gridY[::dr,::dr], vx[::dr,::dr], vy[::dr,::dr],
                        color=self.vector_props['color'],
                        scale=self.vector_props['scale'],
                        pivot=self.vector_props['pivot'],
                        width=self.vector_props['width'],
                        )
    

    def create_gif(self):
        import imageio.v2 as imageio
        r.report("Creating animation...")
        fps = 5 if self.eco_plot else 10
        #img_paths = [self.output_plot_dir / f"{self.full_quantity_name.replace(' ', '')}_{ts[0]}.png" for ts in self.timestrings]
        images = [imageio.imread(p) for p in self.img_paths]
        gif_path = self.output_plot_dir / 'animation.gif'
        imageio.mimsave(str(gif_path), images, fps=fps, loop=0)  # type: ignore[arg-type]
        if self.delete_static_imgs:
            for img in self.img_paths:
                Path.unlink(img)
