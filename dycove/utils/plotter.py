### ---------- Use this class to plot ANUGA/DFM results ---------- ###

import numpy as np
from pathlib import Path
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable

from array_math import cell_averaging, sum_product
from gridding import create_nn_interpFunc
from model_loader import DFMMapLoader, ANUGAMapLoader


### ----- PLOTTING METHODS ----- ###

class ModelPlotter:
    def __init__(self, 
                 simdir, 
                 quantity, 
                 plot_times,
                 mask_bndy_file=None,
                 extents=None,
                 hydro_plot_label_time='hydrodynamic',
                 veg_plot_label_time='eco-morphodynamic',
                 n_ets_year=14,
                 hr_ets=12,
                 vegfac=None,
                 plot_interp_vars=None,
                 plot_specs=None, cmap_lims=None, cmaps=None,
                 animate=False, delete_static_imgs=False,
                 save_grids=False,
                 ):
        
        self.simdir = Path(simdir)
        self.quantity = quantity
        self.quantity_name = ' '.join(quantity.split()[:-1])  # removing the unit label from the name
        self.plot_times = plot_times
        self.mask_bndy_file = mask_bndy_file
        self.extents = extents
        #self.vegfac = vegfac
        self.animate = animate
        self.delete_static_imgs = delete_static_imgs
        self.save_grids = save_grids

        # for restart simulations, some additional temporal inputs
        self.n_ets_year = n_ets_year
        self.hr_ets = hr_ets

        # compute vegfac if not provided (based on n_ets_year and hr_ets)
        self.vegfac, self.days_per_year = self.compute_vegfac(vegfac)

        # only relevant for ANUGA and DFM, must build an interpolation function for converting 1-D arrays to a grid
        self.interp_func = None
            
        self.plot_interp_vars = plot_interp_vars or {
            'plot_cell_size': 10,  # in meters
            'n_neighbors'   : 1,   # typically 1 or 3
        }

        self.plot_specs = plot_specs or {
            'figsize'       : (6, 6),
            'fontsize'      : 12,
            'output_dpi'    : 150,
        }

        self.cmap_lims = cmap_lims or {
            'Bathymetry'      : (0, 1.5),
            'WSE'             : (-0.1, 0.1),
            'Depth'           : (0, 3),
            'Velocity'        : (0, 0.5),
            'Max Shear Stress': (0, 0.1),
            'Fractions'       : (0, 1.0),
            'Stem Height'     : (0, 2.0),
            'Stem Density'    : (0, 300),
            'Stem Diameter'   : (0, 0.1),
            'Stem Face Factor': (0, 15),
            'Mortality -- Flooding' : (0, 1.0),
            'Mortality -- Uprooting': (0, 1.0),
            'Mortality -- Total': (0, 1.0),
        }

        self.cmaps = cmaps or {
            'Bathymetry'      : ListedColormap((np.linspace(0, 1.0, 256))),
            #'Bathymetry'      : ListedColormap(plt.colormaps["bone"](np.linspace(0, 1.0, 256))),
            'WSE'             : ListedColormap(plt.colormaps["Blues_r"](np.linspace(0.1, 1.0, 256))),
            'Depth'           : ListedColormap(plt.colormaps["Blues_r"](np.linspace(0.1, 1.0, 256))),
            'Velocity'        : ListedColormap(plt.colormaps["viridis"](np.linspace(0, 1.0, 256))),
            'Max Shear Stress': ListedColormap(plt.colormaps["plasma"](np.linspace(0, 1.0, 256))),
            'Fractions'       : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Height'     : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Density'    : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Diameter'   : ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Stem Face Factor': ListedColormap(plt.colormaps["Greens"](np.linspace(0.2, 1.0, 256))),
            'Mortality -- Flooding' : ListedColormap(plt.colormaps["Oranges"](np.linspace(0.2, 1.0, 256))),
            'Mortality -- Uprooting': ListedColormap(plt.colormaps["Oranges"](np.linspace(0.2, 1.0, 256))),
            'Mortality -- Total'    : ListedColormap(plt.colormaps["Oranges"](np.linspace(0.2, 1.0, 256))),
            }

        # hardcoded, for masking out quantities when they are small enough
        self.min_depth  = 0.1
        self.min_tau    = 0.001
        self.min_veg_lims = {
            'Fractions': 0.01,
            'Stem Height'  : 0.1,
            'Stem Density' : 10,
            'Stem Diameter': 0.01,
            'Stem Face Factor': 0.025,
            'Mortality -- Flooding': 0.01,
            'Mortality -- Uprooting': 0.01,
            'Mortality -- Total': 0.01, 
        }

        self.output_plot_dir = self.simdir / 'figures' / self.quantity_name.replace(' ', '')
        self.model_type, self.many_dirs = self.get_model_type()
        self.modeldir = self.get_modeldir()
        self.model_name = self.get_model_name()
        self.veg_plot = self.quantity_name in ['Fractions', 'Stem Height', 'Stem Density', 'Stem Diameter'] or \
                            'Mortality' in self.quantity_name

        # change how we label time in our plots depending on if we are plotting vegetation (uses vegfac) or hydrodynamics (no morfac/vegfac)
        self.plot_label_time = veg_plot_label_time if self.veg_plot else hydro_plot_label_time

        # load model output files using the appropriate loader class
        if self.model_type == 'DFM': 
            self.map_loader = DFMMapLoader(self.modeldir, self.model_name, 
                                           self.quantity_name, self.veg_plot, self.n_ets_year)
        elif self.model_type == 'ANUGA': 
            self.map_loader = ANUGAMapLoader(self.modeldir, self.model_name, 
                                             self.quantity_name, self.veg_plot, self.n_ets_year)
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
        if self.veg_plot:
            sim_hours = float(i*self.plotHR_int)
        elif self.many_dirs:
            sim_hours = float(self.plotHR_0 + i*self.plotHR_div)
        else:
            sim_hours = float(i*self.plotHR_div)
        
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
        # how we loop over the model time slices will depend on the directory structure and whether we are plotting hydro or veg quantities
        # first need to determine the year, ets, and hydro_ind based on the value of i and if this run has many directories and veg plotting
        # veg_ind will always be i, and will always be a simple sequence of integers starting at 1 (if not veg_plot, then it doesn't matter)
        hydro_ind, ets, eco_year = i, None, None
        # if self.veg_plot and i == 0:
        #     print("No vegetation at time step 0")
        # else:
        if self.veg_plot:
            ets, eco_year = self.get_veg_eco_times(i)
            hydro_ind = i*self.hr_ets
        
        return hydro_ind, ets, eco_year


    def get_veg_eco_times(self, i):
        ets = ((i-1) % self.n_ets_year) + 1
        eco_year = ((i-1) // self.n_ets_year) + 1
        return ets, eco_year
    
    # def get_ETS_indices(self, n_ets_year, hr_ets):
    #     # returning: total number of ETS to plot, total number of eco years to plot (in case it's not all of them)
    #     ets2plot   = int(math.ceil(self.plot_times['plotHR_f']/hr_ets))
    #     years2plot = int(math.ceil(ets2plot/n_ets_year))
    #     return ets2plot, years2plot


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
        if self.quantity_name == 'Bathymetry':
            return z_grid, None
        elif not self.veg_plot:
            depth_grid = self.interp_func(map_vars['Depth'])
            if self.quantity_name == 'Depth':
                return z_grid, np.ma.masked_where(depth_grid < self.min_depth, depth_grid)
            elif self.quantity_name == 'WSE':
                wse_grid = self.interp_func(map_vars[self.quantity_name])
                return z_grid, np.ma.masked_where(depth_grid < self.min_depth, wse_grid)
            elif self.quantity_name == 'Velocity':
                V_grid = self.interp_func(map_vars[self.quantity_name])
                return z_grid, np.ma.masked_where(depth_grid < self.min_depth, V_grid)
            elif self.quantity_name == 'Max Shear Stress':
                tau_grid = self.interp_func(map_vars[self.quantity_name])
                return z_grid, np.ma.masked_where(tau_grid < self.min_tau, tau_grid)
            else:
                raise ValueError(f"Unexpected quantity name: {self.quantity_name}")
        else:
            if len(map_vars[self.quantity_name]) > 0:
                if self.quantity_name[:9] in ['Fractions', 'Mortality']:
                    fraction_grid_list = []
                    for frac in map_vars[self.quantity_name]:
                        veg_grid = self.interp_func(frac)
                        fraction_grid_list.append(np.ma.masked_where(veg_grid < self.min_veg_lims[self.quantity_name], veg_grid))
                    return z_grid, fraction_grid_list
                elif self.quantity_name == 'Stem Density':  # sum_product for 'Stem Density'
                    veg_data = sum_product(map_vars['Fractions'], map_vars[self.quantity_name])
                else:  # cell_averaging for 'Stem Diameter', 'Stem Face Factor', and 'Stem Height'
                    veg_data = cell_averaging(map_vars['Fractions'], map_vars[self.quantity_name])
                veg_grid = self.interp_func(veg_data)
                return z_grid, np.ma.masked_where(veg_grid < self.min_veg_lims[self.quantity_name], veg_grid)
            else:
                # may not be any veg data at all in this time step, so return an empty grid
                empty_grid = np.ma.masked_all(z_grid.shape)
                return z_grid, empty_grid


    def plot_quantity(self, base_grid, main_grid):
        if type(main_grid) is list:  # for veg fractions, we plot each fraction separately
            for i, grid in enumerate(main_grid):
                self._plot_single_quantity(base_grid, grid, 
                        title=f"{self.quantity_name} -- {self.timestrings[-1][1]} -- i={i}",
                        fname=f"{self.quantity_name.replace(' ', '')}_{self.timestrings[-1][0]}_i={i}")
        else:
            self._plot_single_quantity(base_grid, main_grid, 
                        title=f"{self.quantity_name} -- {self.timestrings[-1][1]}",
                        fname=f"{self.quantity_name.replace(' ', '')}_{self.timestrings[-1][0]}")
            
            
    def _plot_single_quantity(self, base_grid, main_grid, title, fname):
        fig, ax = plt.subplots(figsize=self.plot_specs['figsize'])

        base = ax.imshow(base_grid,
                         cmap=self.cmaps['Bathymetry'],
                         vmin=self.cmap_lims['Bathymetry'][0],
                         vmax=self.cmap_lims['Bathymetry'][1])
        if main_grid is not None:
            main = ax.imshow(main_grid,
                            cmap=self.cmaps[self.quantity_name],
                            vmin=self.cmap_lims[self.quantity_name][0],
                            vmax=self.cmap_lims[self.quantity_name][1])            

        ax.set_title(title, fontsize=self.plot_specs['fontsize'])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect('equal')

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.12)
        colorbar = plt.colorbar(main if self.quantity_name != 'Bathymetry' else base, cax=cax)
        colorbar.set_label(self.quantity, rotation=270, labelpad=25, fontsize=self.plot_specs['fontsize'])
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
        img_paths = [self.output_plot_dir / f"{self.quantity_name.replace(' ', '')}_{ts[0]}.png" for ts in self.timestrings]
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