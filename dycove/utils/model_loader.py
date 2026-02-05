
from abc import ABC, abstractmethod
from typing import Mapping, Any
import numpy as np
from pathlib import Path
import xarray as xr


# ------------------------------------------------------------ #
# Some utility functions first, for internal and external use
# ------------------------------------------------------------ #

def get_veg_file_index(eco_dir: Path | str) -> dict:
    """ 
    Lightweight first pass to read eco year & ETS from output and generate file index.

    Parameters
    ----------
    eco_dir : Path or str
        Path to DYCOVE's ``veg_output`` directory.
     
    """
    eco_dir = Path(eco_dir)
    index = {}
    old_file_format = list(eco_dir.glob(".npz"))
    if not old_file_format:
        for fpath in eco_dir.iterdir():
            with xr.open_dataset(fpath) as ds:
                year = ds.attrs["eco_year"]
                ets = ds.attrs["ets"]
            if (year, ets) not in index:
                index[(year, ets)] = [fpath]
            else:
                index[(year, ets)].append(fpath)
    # return empty index if using old format, won't be used
    return index


def get_anuga_centroid_coords(anuga_vars: xr.Dataset) -> tuple[list, list]:
    """ 
    Convert ANUGA vertex coordinates to centroid coordinates (centroid coordinates 
    not available in ANUGA sww file).

    Parameters
    ----------
    anuga_vars : xr.Dataset
        Dataset containing netCDF variables.
    """
    # Load mesh vertex coordinates
    xx = anuga_vars['x']
    yy = anuga_vars['y']
    # Convert vertex coordinates to centroid coordinates using 'volumes' variable
    xx_c = [(xx[i]+xx[j]+xx[k])/3. for i, j, k in anuga_vars['volumes']]
    yy_c = [(yy[i]+yy[j]+yy[k])/3. for i, j, k in anuga_vars['volumes']]

    return xx_c, yy_c


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

        # Single location for storing names of vegetation variables stored in class `VegCohort`
        # DYCOVE related and independent on the numerical model being used
        self.veg_varnames = {'Fractions': 'fraction', 
                             'Stem Density': 'density', 
                             'Stem Diameter': 'diameter', 
                             'Stem Height': 'height',
                             'Species': 'name',
                             'Cohort': 'cohort',
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
        
        if self.eco_plot:
            self.veg_file_index = get_veg_file_index(self.ecodir)

        
    @abstractmethod
    def load(self, hydro_i, ets, eco_year):
        pass

    @abstractmethod
    def _load_outputs(self, modeldir):
        pass
    
    def check_final_index(self, index):
        """ Raises error if final plot time exceeds simulation length """
        self._load_outputs(self._get_model_subdir())
        wse = self.cached_map_vars[self.hydro_varnames['WSE']]
        if index >= len(wse):
            raise ValueError("Final plot time in 'plot_times' exceeds the length of the simulation!")
        
    def _load_veg(self, ets, eco_year):
        """
        Read vegetation files, this function abstracts similar logic among DFM and ANUGA
        No. of vegetation cohorts present is equal to the eco year we are plotting
        No cached variable here because these files are written every ETS
        """

        veg_data = {"fractions": [],
                    "quantities": [],
                    "cohorts": [],
                    }

        old_file_format = list(self.ecodir.glob(".npz"))
        if old_file_format:
            return self._load_veg_old_format(veg_data, ets, eco_year)
        else:
            return self._load_veg_new_format(veg_data, ets, eco_year)
        
        
    def _load_veg_new_format(self, veg_data, ets, eco_year):
        if (eco_year, ets) in self.veg_file_index:
            for fpath in self.veg_file_index[(eco_year, ets)]:
                c = xr.load_dataset(fpath)

                veg_data["cohorts"].append((c.attrs[self.veg_varnames["Species"]], 
                                            c.attrs[self.veg_varnames["Cohort"]]))
                veg_data["fractions"].append(c.data_vars[self.veg_varnames["Fractions"]])
                if self.veg_varnames[self.quantity] in c.attrs:
                    veg_data["quantities"].append(c.attrs[self.veg_varnames[self.quantity]])
                else:  # redundant if quantity is 'Fractions'
                    veg_data["quantities"].append(c.data_vars[self.veg_varnames[self.quantity]])

        return veg_data


    def _load_veg_old_format(self, veg_data, ets, eco_year):
        # Loop through all cohort files saved for this ETS, load attributes to running list
        for file in self.ecodir.glob(f"cohort*_year{eco_year}_ets{ets}.npz"):
            c = dict(np.load(file, allow_pickle=True))
            c_id = file.stem.split("_")[0][6:]  # old format didn't have species name, give generic name
            veg_data["cohorts"].append(f"Cohort {c_id}") 
            #veg_data["species"].append(c[self.veg_varnames["Species"]])
            veg_data["fractions"].append(c[self.veg_varnames["Fractions"]])
            veg_data["quantities"].append(c[self.veg_varnames[self.quantity]])  # not used if quantity is 'Fractions'

        return veg_data


    def _pass_veg(self, veg_data, data):
        # Distribute veg_data to correct keys in data dict
        data["Cohort Names"] = veg_data["cohorts"]
        data["Fractions"] = veg_data["fractions"]
        # For all other quantities, we still need Fractions in order to do weighted averaging of cohorts in grid cells
        data[self.quantity] = veg_data["quantities"]  # if quantity is Fractions, this line is redundant
        return data


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
                'Bathymetry': np.asarray(self.zz_c),
                'Cohort Names': None  # need this for output consistency
                }
        
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

            self.cached_map_vars = xr.load_dataset(subdir / f'{self.model_name}.sww')
            assert self.cached_map_vars is not None  # for Pylance...

            # Only create the mesh centroid variables if they don't exist for this mesh
            # -- this is a time consuming step, so save these files for future plotting with the same mesh
            files_exist = False
            x_fname = "xCentroidsSavedForFastRecomputation.npy"
            y_fname = "yCentroidsSavedForFastRecomputation.npy"    
            if (subdir / x_fname).exists() and (subdir / y_fname).exists():
                self.xx_c = np.load(subdir / x_fname)
                self.yy_c = np.load(subdir / y_fname)
                # Make sure existing files are not leftover from previous run
                if len(self.xx_c) == len(self.cached_map_vars[self.mesh_varnames['Z']]):
                    files_exist = True
            if not files_exist:
                self.xx_c, self.yy_c = get_anuga_centroid_coords(self.cached_map_vars)
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
                'Y': np.asarray(self.cached_map_vars[self.mesh_varnames['Y']]),
                'Cohort Names': None  # need this for output consistency
                }
        
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
            self.cached_map_vars = xr.load_dataset(subdir / f'{self.model_name}_map.nc')
            self.current_subdir = subdir

    def _get_model_subdir(self):
        # Needed because DYCOVE-DFM creates a separate folder containing the output
        return self.modeldir / 'output'