"""
Class for plotting various output quantites from DYCOVE-ANUGA/DFM
"""

from typing import Optional, Any, Union
import numpy as np
from pathlib import Path
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap, ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

from dycove.utils.array_math import cell_averaging, sum_product
from dycove.utils.gridding import create_nn_interpFunc
from dycove.utils.model_loader import DFMMapLoader, ANUGAMapLoader


class ModelPlotter:
    def __init__(self, 
                 simdir, 
                 quantity, 
                 plot_times,
                 hydro_plot_label_time='hydrodynamic',
                 veg_plot_label_time='eco-morphodynamic',
                 n_ets_year=14,
                 hr_ets=12,
                 vegfac=None,
                 plot_interp_vars: Optional[dict[str, int]] = None,
                 plot_specs: Optional[dict[str, Any]] = None, 
                 cmap_lims: Optional[dict[str, tuple[float, float]]] = None, 
                 cmaps: Optional[dict[str, Colormap]] = None,
                 quantity_units: Optional[dict[str, tuple[str, Union[int, float]]]] = None, 
                 mask_bndy_file=None,
                 extents=None,
                 animate=False, 
                 delete_static_imgs=False,
                 save_grids=False,
                 ):
        
        self.simdir = Path(simdir)
        self.quantity = quantity
        self.plot_times = plot_times
        self.mask_bndy_file = mask_bndy_file
        self.extents = extents
        self.animate = animate
        self.delete_static_imgs = delete_static_imgs
        self.save_grids = save_grids

        # simplify quantity name if "Mortality"
        self.full_quantity_name = self.quantity
        if self.quantity[:9] == "Mortality":
            self.quantity = "Mortality"

        # for restart simulations, some additional temporal inputs
        self.n_ets_year = n_ets_year
        self.hr_ets = hr_ets

        # compute vegfac if not provided (based on n_ets_year and hr_ets)
        self.vegfac, self.days_per_year = self.compute_vegfac(vegfac)

        # only relevant for ANUGA and DFM, must build an interpolation function for converting 1-D arrays to a grid
        self.interp_func = None
            
        defaults = {
            'plot_cell_size': 5,  # in meters
            'n_neighbors'   : 1,   # typically 1 or 3
        }
        self.plot_interp_vars = {**defaults, **(plot_interp_vars or {})}  # merge provided custom values with default values

        defaults = {
            'figsize'       : (6, 6),
            'fontsize'      : 12,
            'output_dpi'    : 100,
        }
        self.plot_specs = {**defaults, **(plot_specs or {})}  # merge provided custom values with default values

        defaults = {
            'Bathymetry'      : (0, 1),
            'WSE'             : (0, 1),
            'Depth'           : (0, 3),
            'Velocity'        : (0, 0.5),
            'Max Shear Stress': (0, 0.1),
            'Fractions'       : (0, 1),
            'Stem Height'     : (0, 2),
            'Stem Density'    : (0, 300),
            'Stem Diameter'   : (0, 50),
            'Mortality'       : (0, 1.0),
        }
        self.cmap_lims = {**defaults, **(cmap_lims or {})}  # merge provided custom values with default values

        defaults = {
            'Bathymetry'      : ListedColormap(plt.colormaps["Greys_r"](np.linspace(0, 1.0, 256))),
            #'Bathymetry'      : ListedColormap(plt.colormaps["bone"](np.linspace(0, 1.0, 256))),
            'WSE'             : ListedColormap(plt.colormaps["Blues_r"](np.linspace(0.1, 1.0, 256))),
            'Depth'           : ListedColormap(plt.colormaps["Blues_r"](np.linspace(0.1, 1.0, 256))),
            'Velocity'        : ListedColormap(plt.colormaps["viridis"](np.linspace(0, 1.0, 256))),
            'Max Shear Stress': ListedColormap(plt.colormaps["plasma"](np.linspace(0, 1.0, 256))),
            'Fractions'       : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Height'     : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Density'    : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Diameter'   : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Mortality' : ListedColormap(plt.colormaps["Oranges"](np.linspace(0.2, 1.0, 256))),
        }
        self.cmaps = {**defaults, **(cmaps or {})}  # merge provided custom values with default values

        # units to use for plot title/colorbar, second value is multiplier based on most standard units
        defaults = {
            'Bathymetry'      : ('[m]', 1),
            'WSE'             : ('[m]', 1),
            'Depth'           : ('[m]', 1),
            'Velocity'        : ('[m/s]', 1),
            'Max Shear Stress': ('[N/m$^2$]', 1),
            'Fractions'       : ('[--]', 1),
            'Stem Height'     : ('[m]', 1),
            'Stem Density'    : ('[stems/m$^2$]', 1),
            'Stem Diameter'   : ('[mm]', 1000),
            'Mortality'       : ('[%]', 100),
        }
        self.quantity_units = {**defaults, **(quantity_units or {})}  # merge provided custom values with default values

        # hardcoded, for masking out quantities when they are small enough
        self.min_depth = 0.01
        self.min_tau   = 0.001
        self.min_veg_lims = {
            'Fractions'    : 0.01,
            'Stem Height'  : 0.1,
            'Stem Density' : 10,
            'Stem Diameter': 0.01,
            'Mortality'    : 0.01,
        }

        self.output_plot_dir = self.simdir / 'figures' / self.quantity.replace(' ', '')
        self.model_type = self.get_model_type()
        self.modeldir = self.get_modeldir()
        self.model_name = self.get_model_name()
        self.veg_plot = self.quantity in ['Fractions', 'Stem Height', 'Stem Density', 'Stem Diameter', 'Mortality']

        # change how we label time in our plots depending on if we are plotting vegetation (uses vegfac) or hydrodynamics (no morfac/vegfac)
        self.plot_label_time = veg_plot_label_time if self.veg_plot else hydro_plot_label_time

        # load model output files using the appropriate loader class
        args = (self.modeldir, self.model_name, self.full_quantity_name, self.veg_plot, self.n_ets_year)
        if self.model_type == 'DFM': 
            self.map_loader = DFMMapLoader(*args)
        elif self.model_type == 'ANUGA': 
            self.map_loader = ANUGAMapLoader(*args)

        # to save for building animations later
        self.timestrings = []


    def run(self):
        ti_0, ti_f, dti = self.get_time_indices()
        for i in tqdm(range(ti_0, ti_f, dti)):
            self.create_timestrings(i)
            hydro_i, ets, eco_year = self.get_map_indices(i)
            map_vars = self.map_loader.load(hydro_i, ets, eco_year)
            z_grid, grid2plot = self.get_quantity_grids(map_vars)
            self.plot_quantity(z_grid, grid2plot)

        if self.animate: self.create_gif()


    def compute_vegfac(self, vegfac_input):
        if vegfac_input is None:
            # Compute ideal (continuous) vegfac
            vegfac_est = (24 * 365.) / (self.hr_ets * self.n_ets_year)
            # Round to nearest whole number
            vegfac = int(round(vegfac_est))
            msg = f"Computed vegfac = {vegfac:d}, derived from veg_interval and n_ets_year (rounded from {vegfac_est:.3f})"
        else:
            vegfac = vegfac_input
            msg = f"Using provided vegfac = {vegfac:d}"

        # compute associated days_per_year after vegfac rounding
        days_per_year = (vegfac * self.hr_ets * self.n_ets_year) / 24.

        print(msg)
        print(f"vegfac = {vegfac:d} corresponds to {days_per_year:.1f} days per year.")

        return vegfac, days_per_year
        

    def create_timestrings(self, i):
        time_parts = self.get_time_breakdown(i)
        fname_timestring = self.format_fname_time(time_parts)
        title_timestring = self.format_title_time(time_parts)
        self.timestrings.append((fname_timestring, title_timestring))

    def get_time_breakdown(self, i):
        """Returns a dict with time components for formatting titles and filenames."""
        sim_hours = float(i*self.plotHR_int if self.veg_plot else i*self.plotHR_div)
        
        veg_years = sim_hours*self.vegfac/24/self.days_per_year
        veg_days_rem = (veg_years % 1)*self.days_per_year  # days after full years
        veg_days_tot = sim_hours*self.vegfac/24.

        return {"sim_days": int(sim_hours // 24),
                "sim_hrs_rem": int(sim_hours % 24),
                "sim_mins_rem": int(sim_hours*60 % 60),  # for if we have a plot interval less than one hour, to avoid duplicate file names
                "veg_years": int(veg_years),
                "veg_days_rem": int(round(veg_days_rem)),
                "veg_days_tot": int(round(veg_days_tot))}
    
    def format_fname_time(self, time_parts):
        # always use hydrodynamic time for the filenames
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
        if self.veg_plot:
            self.plotHR_0 = max(self.plot_times['plotHR_0'], self.hr_ets)
            self.plotHR_int = self.hr_ets
            self.plotHR_div = self.hr_ets
        else:
            self.plotHR_0 = self.plot_times['plotHR_0']
            self.plotHR_int = self.plot_times['plotHR_int']
            self.plotHR_div = self.plot_times['mapHR_int']

        ti_0 = int(self.plotHR_0/self.plotHR_div)
        # adding 1 to include the final time step. TODO: check this
        ti_f = int(self.plot_times['plotHR_f']/self.plotHR_div) + 1
        dti = int(self.plotHR_int/self.plotHR_div)

        return ti_0, ti_f, dti


    def get_map_indices(self, i):
        # how we loop over the model time slices will depend on whether we are plotting hydro or veg quantities
        # first need to determine the year, ets, and hydro_ind based on the value of i
        # hydro_ind is simple counter for hydro quantities, but must be multiplied by hr_ets for veg quantities
        hydro_ind, ets, eco_year = i, None, None

        if self.veg_plot:
            ets, eco_year = self.get_veg_eco_times(i)
            hydro_ind = i*self.hr_ets
        
        return hydro_ind, ets, eco_year


    def get_veg_eco_times(self, i):
        ets = ((i-1) % self.n_ets_year) + 1
        eco_year = ((i-1) // self.n_ets_year) + 1
        return ets, eco_year


    def get_hydro_eco_times(self, i):
        ets = ((i // self.hr_ets) % self.n_ets_year) + 1
        eco_year = (i // (self.hr_ets * self.n_ets_year)) + 1
        return ets, eco_year

    def get_model_type(self):
        # determine from files in model directory whether ththe model is DFM or ANUGA
        return 'DFM' if (self.simdir/'dflowfm').exists() else 'ANUGA'

    def get_modeldir(self):
        return self.simdir/'dflowfm' if self.model_type=='DFM' else self.simdir

    def get_model_name(self):
        # Find file with the model file extension, this could be ".mdu" or "_orig.mdu" for DFM, ".mdf" for D3D4, ".sww" for ANUGA
        if self.model_type == 'DFM':
            model_file = next(self.modeldir.glob(f"*.mdu"))
        elif self.model_type == 'ANUGA':
            model_file = next(self.modeldir.glob(f"*.sww"))
        return model_file.stem


    def create_interp_func(self, map_vars):
        # create interpolation function the first time -> quick interpolations for all other time steps
        interp_func = create_nn_interpFunc(map_vars['X'], map_vars['Y'], 
                                           grid_size=self.plot_interp_vars['plot_cell_size'], 
                                           k_nn=self.plot_interp_vars['n_neighbors'],
                                           polygon_csv=self.mask_bndy_file, extents=self.extents,
                                           )
        return interp_func
        

    def get_quantity_grids(self, map_vars):
        # create interpolation function the first time -> quick interpolations for all other time steps
        if self.interp_func is None:
            self.interp_func = self.create_interp_func(map_vars)
            
        # for DFM, z_grid should be recomputed every step in case morphology is turned on. For ANUGA, oh well it's fast enough
        z_grid = self.interp_func(map_vars['Bathymetry'])
        if self.quantity == 'Bathymetry':
            return z_grid, None
        elif not self.veg_plot:
            depth_grid = self.interp_func(map_vars['Depth'])
            if self.quantity == 'Depth':
                return z_grid, np.ma.masked_where(depth_grid < self.min_depth, depth_grid)
            elif self.quantity == 'WSE':
                wse_grid = self.interp_func(map_vars[self.quantity])
                return z_grid, np.ma.masked_where(depth_grid < self.min_depth, wse_grid)
            elif self.quantity == 'Velocity':
                V_grid = self.interp_func(map_vars[self.quantity])
                return z_grid, np.ma.masked_where(depth_grid < self.min_depth, V_grid)
            elif self.quantity == 'Max Shear Stress':
                tau_grid = self.interp_func(map_vars[self.quantity])
                return z_grid, np.ma.masked_where(tau_grid < self.min_tau, tau_grid)
            else:
                raise ValueError(f"Unexpected quantity name: {self.quantity}")
        else:
            if len(map_vars[self.quantity]) > 0:
                if self.quantity in ['Fractions', 'Mortality']:
                    fraction_grid_list = []
                    for frac in map_vars[self.quantity]:
                        veg_grid = self.interp_func(frac)
                        fraction_grid_list.append(np.ma.masked_where(veg_grid < self.min_veg_lims[self.quantity], veg_grid))
                    return z_grid, fraction_grid_list
                elif self.quantity == 'Stem Density':  # sum_product for 'Stem Density'
                    veg_data = sum_product(map_vars['Fractions'], map_vars[self.quantity])
                else:  # cell_averaging for 'Stem Diameter', 'Stem Face Factor', and 'Stem Height'
                    veg_data = cell_averaging(map_vars['Fractions'], map_vars[self.quantity])
                veg_grid = self.interp_func(veg_data)
                return z_grid, np.ma.masked_where(veg_grid < self.min_veg_lims[self.quantity], veg_grid)
            else:
                # may not be any veg data at all in this time step, so return an empty grid
                empty_grid = np.ma.masked_all(z_grid.shape)
                return z_grid, empty_grid


    def plot_quantity(self, base_grid, main_grid):
        if type(main_grid) is list:  # for veg fractions, we plot each fraction separately
            for i, grid in enumerate(main_grid):
                self._plot_single_quantity(base_grid, grid, 
                        title=f"{self.full_quantity_name} -- {self.timestrings[-1][1]} -- i={i}",
                        fname=f"{self.full_quantity_name.replace(' ', '')}_{self.timestrings[-1][0]}_i={i}")
        else:
            self._plot_single_quantity(base_grid, main_grid, 
                        title=f"{self.full_quantity_name} -- {self.timestrings[-1][1]}",
                        fname=f"{self.full_quantity_name.replace(' ', '')}_{self.timestrings[-1][0]}")
            
            
    def _plot_single_quantity(self, base_grid, main_grid, title, fname):

        fig, ax = plt.subplots(figsize=self.plot_specs['figsize'])

        base_grid *= self.quantity_units['Bathymetry'][1]
        base = ax.imshow(base_grid,
                         cmap=self.cmaps['Bathymetry'],
                         vmin=self.cmap_lims['Bathymetry'][0],
                         vmax=self.cmap_lims['Bathymetry'][1])
        
        if main_grid is not None:
            main_grid *= self.quantity_units[self.quantity][1]
            main = ax.imshow(main_grid,
                             cmap=self.cmaps[self.quantity],
                             vmin=self.cmap_lims[self.quantity][0],
                             vmax=self.cmap_lims[self.quantity][1])            

        ax.set_title(title, fontsize=self.plot_specs['fontsize'])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.12)
        colorbar = plt.colorbar(main if self.quantity != 'Bathymetry' else base, cax=cax)
        colorbar.set_label(f"{self.quantity} {self.quantity_units[self.quantity][0]}", rotation=270, labelpad=25, fontsize=self.plot_specs['fontsize'])
        colorbar.ax.tick_params(labelsize=self.plot_specs['fontsize'])

        fig.tight_layout()

        self.output_plot_dir.mkdir(parents=True, exist_ok=True)
        fig_path = self.output_plot_dir / fname
        plt.savefig(f"{fig_path}.png", dpi=self.plot_specs['output_dpi'], bbox_inches='tight')
        if self.save_grids:
            np.savez_compressed(f"{fig_path}.npz", data=main_grid.data, mask=main_grid.mask)
        plt.close()


    def create_gif(self):
        import imageio.v2 as imageio
        print("Creating animation...")
        fps = 5 if self.veg_plot else 10
        img_paths = [self.output_plot_dir / f"{self.full_quantity_name.replace(' ', '')}_{ts[0]}.png" for ts in self.timestrings]
        images = [imageio.imread(p) for p in img_paths]
        gif_path = self.output_plot_dir / 'animation.gif'
        imageio.mimsave(str(gif_path), images, fps=fps, loop=0)  # type: ignore[arg-type]
        if self.delete_static_imgs:
            for img in img_paths:
                Path.unlink(img)



### ----- EXECUTING THE ABOVE ----- ###

if __name__ == "__main__":
    plotter = ModelPlotter(
                #simdir = Path('C:/Users/neltull/Documents/ANUGA/WLD_local/sim01'),
                simdir = Path('C:/Users/neltull/Documents/Delft3D/BijDeVaateModel/DFM/sim50'),
                # --- 'Bathymetry [m]', 'WSE [m]', 'Depth [m]', 'Velocity [m/s]', 
                # --- 'Fractions [-]', 'Stem Height [m]', 'Mortality [-]', 'Stem Density [stems/m$^2$]' (DFM only)
                #quantity = 'Fractions [-]',
                #quantity = 'Mortality -- Flooding [-]',
                #quantity = 'Mortality -- Uprooting [-]',
                #quantity = 'Mortality -- Total [-]',
                #quantity = 'Stem Density [stems/m$^2$]',
                #quantity = 'Stem Face Factor [m$^(-1)$]',
                quantity = 'Stem Height [m]',
                #quantity = 'Velocity [m/s]',
                #quantity = 'WSE [m]',
                #quantity = 'Depth [m]',
                #quantity = 'Bathymetry [m]',
                #quantity = 'Max Shear Stress [N/m$^2$]',
                plot_times = {'plotHR_0': 0.,          # sim hr to start plotting
                              'plotHR_f': 21*24.,      # sim hr to stop plotting, not to exceed total sim length
                              'mapHR_int': 1,   # sim hrs between map outputs
                              'plotHR_int': 2,  # hrs between consecutive plots, cannot be less than map_output, unused if plotting vegetation
                              },
                mask_bndy_file = None,
                extents = None,
                animate = False,
                delete_static_imgs = False,
    )

    plotter.run()