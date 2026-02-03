
from abc import ABC, abstractmethod
import numpy as np
from pathlib import Path
import xarray as xr


class BaseMapLoader(ABC):
    """
    Abstract base class for loading map quantities from DYCOVE model outputs.

    Provides shared logic for vegetation quantities and file management between
    model-specific subclasses :class:`~dycove.utils.model_loader.DFMMapLoader` 
    and :class:`~dycove.utils.model_loader.ANUGAMapLoader`, and potentially 
    others in the future.

    Parameters
    ----------
    modeldir : str or Path
        Path to the model’s working directory.
    model_name : str
        Root name of the model file (without extension).
    quantity : str
        Quantity to load (e.g., ``'Velocity'``, ``'Mortality -- Flooding'``).
    eco_plot : bool
        Whether vegetation quantities are being plotted.
    n_ets_year : int
        Number of eco-time-steps per simulated year.

    Notes
    -----
    - Subclasses must implement the :meth:`load` method to return a dictionary of
      arrays containing the required fields (e.g. ``X``, ``Y``, ``Bathymetry``, and the 
      target ``quantity``).
    - As expected, loading vegetation data/DYCOVE outputs relies on logic independent
      of the numerical model being used. But even in the case of plotting vegetation,
      we would also want to plot bathymetry as a layer underneath, which is dependent
      on numerical model load methods.

    Example
    -------
    >>> loader = ANUGAMapLoader(modeldir, model_name, 'Depth', False, 14)
    >>> data = loader.load(hydro_i=0, ets=None, eco_year=None)
    >>> print(data.keys())
    dict_keys(['X', 'Y', 'Bathymetry', 'WSE', 'Depth'])
    """

    def __init__(self, modeldir, model_name, quantity_name, eco_plot, n_ets_year):
        self.modeldir = Path(modeldir)
        self.model_name = model_name
        self.quantity = quantity_name
        self.eco_plot = eco_plot
        self.n_ets_year = n_ets_year
        self.ecodir = self.modeldir / 'veg_output'

        # single location for storing names of vegetation variables stored in class `VegCohort`
        # since this is DYCOVE related, this isn't dependent on the numerical model being used
        self.veg_varnames = {'Fractions': 'fraction', 
                             'Stem Density': 'density', 
                             'Stem Diameter': 'diameter', 
                             'Stem Height': 'height',
                             'Potential Mortality -- Flooding': 'potential_mort_flood',
                             'Potential Mortality -- Desiccation': 'potential_mort_desic',
                             'Potential Mortality -- Uprooting': 'potential_mort_uproot',
                             'Potential Mortality -- Burial': 'potential_mort_burial',
                             'Potential Mortality -- Scour': 'potential_mort_scour',
                             'Mortality -- Flooding': 'applied_mort_flood',
                             'Mortality -- Desiccation': 'applied_mort_desic',
                             'Mortality -- Uprooting': 'applied_mort_uproot',
                             'Mortality -- Burial': 'applied_mort_burial',
                             'Mortality -- Scour': 'applied_mort_scour',
                             'Mortality -- Total': 'applied_mort_total',
                             }   
        
    @abstractmethod
    def load(self, hydro_i, ets, eco_year):
        pass

    @abstractmethod
    def _load_outputs(self, modeldir):
        pass

    def check_final_index(self, index):
        """ Raises error if final plot time exceeds simulation length """
        try:
            self._load_outputs(self._get_model_subdir())
            wse = self.cached_map_vars[self.hydro_varnames['WSE']][index]
        except IndexError as exc:
            raise ValueError("Final plot time in 'plot_times' exceeds the length of the simulation!") from exc
        
    def _load_veg(self, ets, eco_year):
        """
        Read vegetation files, this function abstracts similar logic among DFM and ANUGA
        No. of vegetation cohorts present is equal to the eco year we are plotting
        No cached variable here because these files are written every ETS
        """

        veg_fractions, veg_quantity, veg_names = [], [], []     

        # Loop through all cohort files saved for this ETS, load file to running list
        for file in self.ecodir.glob(f'cohort*_year{eco_year}_ets{ets}.npz'):
            c = dict(np.load(file, allow_pickle=True))
            veg_names.append(c["name"])
            veg_fractions.append(c["fraction"])
            veg_quantity.append(c[self.veg_varnames[self.quantity]])  # not used if quantity is 'Fractions'

        return veg_fractions, veg_quantity, veg_names

    def _pass_veg(self, veg_data, data):
        # For distributing veg data to correct keys in data dict
        data['Fractions'] = veg_data[0]
        if self.quantity != "Fractions":
            # For all other quantities, we need Fractions in order to do weighted averaging of cohorts in grid cells
            data[self.quantity] = veg_data[1]
        return data, veg_data[2]


class ANUGAMapLoader(BaseMapLoader):
    """
    Loader for ANUGA hydrodynamic and DYCOVE vegetation output files.

    Loads data from ANUGA `.sww` NetCDF files and DYCOVE ``.npz`` vegetation cohort files.

    Notes
    -----
    - ANUGA stores variables at cell centroids, but the ``.sww`` format provides only vertex 
      coordinates; centroids are computed on first load and cached.
    - Subsequent calls reuse saved x and y centroid coordinate arrays because it is an 
      intensive computation for large grids.
    - Velocity is derived from stored quantities Depth and Depth-averaged momentum.

    Returns
    -------
    dict
        Dictionary of NumPy arrays with keys:
        ``'X'``, ``'Y'``, ``'Bathymetry'``, and (if applicable)
        ``'WSE'``, ``'Depth'``, ``'Velocity'``, or vegetation fields.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.current_subdir = None
        self.cached_map_vars = None
        self.cached_veg_mortality = None

        self.mesh_varnames = {'X': 'x', 
                              'Y': 'y', 
                              'Z': 'elevation_c',
                              'triangles': 'volumes',
                              }
        self.hydro_varnames = {'WSE': 'stage_c', 
                               'x-momentum': 'xmomentum_c',
                               'y-momentum': 'ymomentum_c',
                               }
            
    def load(self, hydro_i, ets, eco_year):
        self._load_outputs(self._get_model_subdir())
        assert self.cached_map_vars is not None  # for Pylance...

        data = {'X': np.asarray(self.xx_c),
                'Y': np.asarray(self.yy_c),
                # change from DFM: removed index because no time dimension for elevation in ANUGA
                'Bathymetry': np.asarray(self.zz_c)}
        
        if self.quantity == 'Bathymetry':
            pass

        elif not self.eco_plot:
            data['WSE'] = np.asarray(self.cached_map_vars[self.hydro_varnames['WSE']][hydro_i])
            data['Depth'] = data['WSE'] - data['Bathymetry']
            if self.quantity not in ['WSE', 'Depth']:
                if self.quantity == 'Velocity':
                    xmom = np.asarray(self.cached_map_vars[self.hydro_varnames['x-momentum']][hydro_i])
                    ymom = np.asarray(self.cached_map_vars[self.hydro_varnames['y-momentum']][hydro_i])
                    with np.errstate(divide='ignore', invalid='ignore'):
                        data['Vel_x'] = xmom/data['Depth']
                        data['Vel_y'] = ymom/data['Depth']
                    data['Velocity'] = np.sqrt(data['Vel_x']**2 + data['Vel_y']**2)
        else:
            veg_data = self._load_veg(ets, eco_year)
            data = self._pass_veg(veg_data, data)

        return data

    def _load_outputs(self, subdir):
        # Only load hydro outputs if the subdirectory has changed 
        # (for ANUGA, as of now, there will only be one file to load)
        if self.cached_map_vars is None:

            self.cached_map_vars = xr.open_dataset(subdir / f'{self.model_name}.sww')
            assert self.cached_map_vars is not None  # for Pylance...

            # only create the mesh centroid variables if they don't exist for this mesh
            # -- this is a time consuming step, so save these files for future plotting with the same mesh
            files_exist = False
            x_fname = "xCentroidsSavedForFastRecomputation.npy"
            y_fname = "yCentroidsSavedForFastRecomputation.npy"    
            if (subdir / x_fname).exists() and (subdir / y_fname).exists():
                self.xx_c = np.load(subdir / x_fname)
                self.yy_c = np.load(subdir / y_fname)
                # make sure existing files are not leftover from previous run
                if len(self.xx_c) == len(self.cached_map_vars[self.mesh_varnames['Z']]):
                    files_exist = True
            if not files_exist:
                # load mesh VERTEX coordinates (centroid coordinates not available in ANUGA netCDF file)
                xx = self.cached_map_vars[self.mesh_varnames['X']]
                yy = self.cached_map_vars[self.mesh_varnames['Y']]
                # convert vertex coordinates to centroid coordinates using 'volumes' variable
                self.xx_c = [(xx[i]+xx[j]+xx[k])/3. for i, j, k in self.cached_map_vars[self.mesh_varnames['triangles']]]
                self.yy_c = [(yy[i]+yy[j]+yy[k])/3. for i, j, k in self.cached_map_vars[self.mesh_varnames['triangles']]]
                
                np.save(subdir / x_fname, self.xx_c)
                np.save(subdir / y_fname, self.yy_c)
            
            self.zz_c = self.cached_map_vars[self.mesh_varnames['Z']]

    def _get_model_subdir(self):
        # Needed because DYCOVE-DFM creates a separate folder containing the output
        return self.modeldir
    

class DFMMapLoader(BaseMapLoader):
    """
    Loader for DFM hydrodynamic and DYCOVE vegetation output files.

    Loads data from DFM ``*_map.nc`` NetCDF file and DYCOVE ``.npz`` vegetation cohort
    files.

    Notes
    -----
    - Automatically detects dynamic morphology (``'mesh2d_mor_bl'``).
    - Caches variables for efficient access.

    Returns
    -------
    dict
        Dictionary of NumPy arrays with fields among:
        ``'X'``, ``'Y'``, ``'Bathymetry'``, ``'WSE'``, ``'Depth'``, ``'Velocity'``,
        and vegetation fields.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.current_subdir = None
        self.cached_map_vars = None
        self.cached_veg_mortality = None
        self.cached_species_index = None

        # names of variables in DFM netCDF output files
        self.mesh_varnames = {'X': 'mesh2d_face_x', 
                              'Y': 'mesh2d_face_y', 
                              'Z': 'mesh2d_flowelem_bl',          # static bathymetry
                              'Bathymetry mor': 'mesh2d_mor_bl'}  # dynamic bathymetry
        self.hydro_varnames = {'WSE': 'mesh2d_s1', 
                               'Velocity': ('mesh2d_ucmag', 'mesh2d_ucx', 'mesh2d_ucy'),
                               'Max Shear Stress': 'mesh2d_tausmax'}

    def load(self, hydro_i, ets, eco_year):
        self._load_outputs(self._get_model_subdir())
        assert self.cached_map_vars is not None  # for Pylance...

        # load mesh2d data
        data = {'X': np.asarray(self.cached_map_vars[self.mesh_varnames['X']]),
                'Y': np.asarray(self.cached_map_vars[self.mesh_varnames['Y']])}
        
        # need to handle case where morphology is off, then this variable won't be present in output file
        if self.mesh_varnames['Bathymetry mor'] in self.cached_map_vars:
            data['Bathymetry'] = np.asarray(self.cached_map_vars[self.mesh_varnames['Bathymetry mor']][hydro_i])
        else:
            data['Bathymetry'] = np.asarray(self.cached_map_vars[self.mesh_varnames['Z']])

        if self.quantity == 'Bathymetry':
            pass

        elif not self.eco_plot:
            data['WSE'] = np.asarray(self.cached_map_vars[self.hydro_varnames['WSE']][hydro_i])
            data['Depth'] = data['WSE'] - data['Bathymetry']
            if self.quantity not in ['WSE', 'Depth']:
                if self.quantity == 'Velocity':
                    data['Velocity'] = np.asarray(self.cached_map_vars[self.hydro_varnames[self.quantity][0]][hydro_i])
                    data['Vel_x'] = np.asarray(self.cached_map_vars[self.hydro_varnames[self.quantity][1]][hydro_i])
                    data['Vel_y'] = np.asarray(self.cached_map_vars[self.hydro_varnames[self.quantity][2]][hydro_i])
                else:
                    data[self.quantity] = np.asarray(self.cached_map_vars[self.hydro_varnames[self.quantity]][hydro_i])
        else:
            veg_data = self._load_veg(ets, eco_year)
            data = self._pass_veg(veg_data, data)

        return data
    
    def _load_outputs(self, subdir):
        # Only load hydro outputs if the subdirectory has changed 
        # (for DFM, it will only change in restart mode)
        if subdir != self.current_subdir:
            self.cached_map_vars = xr.open_dataset(subdir / f'{self.model_name}_map.nc')
            self.current_subdir = subdir

    def _get_model_subdir(self):
        # Needed because DYCOVE-DFM creates a separate folder containing the output
        return self.modeldir / 'output'