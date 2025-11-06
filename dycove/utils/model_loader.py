from abc import ABC, abstractmethod
import numpy as np
from pathlib import Path
from netCDF4 import Dataset  # type: ignore


class BaseMapLoader(ABC):
    def __init__(self, modeldir, model_name, quantity, veg_plot, n_ets_year):
        self.modeldir = Path(modeldir)
        self.model_name = model_name
        self.quantity = quantity
        self.veg_plot = veg_plot
        self.n_ets_year = n_ets_year
        self.vegdir = self.modeldir / 'veg_output'

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

    def _load_DFM_ANUGA_veg(self, ets, eco_year):
        # Read vegetation files, this function abstracts similar logic among DFM and ANUGA
        # No. of vegetation cohorts present is equal to the eco year we are plotting
        # No cached variable here because these files are written every ETS
        veg_fractions, veg_quantity = [], []     

        # loop through all cohort files saved for this ETS, load file to running list
        for file in self.vegdir.glob(f'cohort*_year{eco_year}_ets{ets}.npz'):
            c = dict(np.load(file, allow_pickle=True))
            veg_fractions.append(c["fraction"])
            veg_quantity.append(c[self.veg_varnames[self.quantity]])  # not used if quantity is 'Fractions'

        return veg_fractions, veg_quantity

    # def _load_DFM_ANUGA_mortality(self, veg_i):
    #     # Now that mortality is stored with cohort files, this function is parallel to the _load_DFM_ANUGA_veg
    #     veg_fractions, veg_densities, veg_diameters, veg_heights = [], [], [], []

    #     # loop through all cohort files saved for this ETS, load file to running list
    #     for file in self.vegdir.glob(f'cohort*_year{eco_year}_ets{ets}.npz'):
    #         c = dict(np.load(file, allow_pickle=True))
    #         veg_fractions.append(c["fraction"])
    #         veg_densities.append(c["density"])
    #         veg_diameters.append(c["diameter"])
    #         veg_heights.append(c["height"])

    #     return veg_fractions, veg_densities, veg_diameters, veg_heights

    #     # if self.cached_veg_mortality is None:
    #     #     with open(self.vegdir / f'{self.mort_filenames[self.quantity]}.json', 'r') as f:
    #     #         self.cached_veg_mortality = json.load(f)
    #     # return self.cached_veg_mortality[veg_i]

    def _pass_DFM_ANUGA_veg(self, veg_data, data):
        # just for distributing veg data to correct keys in data dict
        data['Fractions'] = veg_data[0]
        if self.quantity != "Fractions":
            # for all other quantities, we need Fractions in order to do weighted averaging of cohorts in grid cells
            data[self.quantity] = veg_data[1]
        return data

class ANUGAMapLoader(BaseMapLoader):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.current_subdir = None
        self.cached_map_vars = None
        self.cached_veg_mortality = None

    def _load_anuga_outputs(self, subdir):
        # only load hydro outputs if the subdirectory has changed 
        # (for ANUGA, as of now, there will only be one file to load)
        if self.cached_map_vars is None:

            self.cached_map_vars = Dataset(subdir / f'{self.model_name}.sww').variables
            assert self.cached_map_vars is not None  # for Pylance...

            # only create the mesh centroid variables if they don't exist for this mesh
            # -- this is a time consuming step, so save these files for future plotting with the same mesh
            files_exist = False
            if (subdir / "xx_c.npy").exists() and (subdir / "yy_c.npy").exists():
                self.xx_c = np.load(subdir / "xx_c.npy")
                self.yy_c = np.load(subdir / "yy_c.npy")
                # make sure existing files are not leftover from previous run
                if len(self.xx_c) == len(self.cached_map_vars['elevation_c']):
                    files_exist = True
            if not files_exist:
                # load mesh VERTEX coordinates (centroid coordinates not available in ANUGA netCDF file)
                xx = self.cached_map_vars['x']
                yy = self.cached_map_vars['y']
                # convert vertex coordinates to centroid coordinates using 'volumes' variable
                self.xx_c = [(xx[i]+xx[j]+xx[k])/3. for i, j, k in self.cached_map_vars['volumes']]
                self.yy_c = [(yy[i]+yy[j]+yy[k])/3. for i, j, k in self.cached_map_vars['volumes']]
                
                np.save(subdir / "xx_c.npy", self.xx_c)
                np.save(subdir / "yy_c.npy", self.yy_c)
            
            self.zz_c = self.cached_map_vars['elevation_c']

            
    def load(self, hydro_i, ets, eco_year):
        self._load_anuga_outputs(self.modeldir)
        assert self.cached_map_vars is not None  # for Pylance...

        data = {'X': np.asarray(self.xx_c),
                'Y': np.asarray(self.yy_c),
                # change from DFM: removed index because no time dimension for elevation in ANUGA
                'Bathymetry': np.asarray(self.zz_c)}
        
        if self.quantity == 'Bathymetry':
            pass

        elif not self.veg_plot:
            data['WSE'] = np.asarray(self.cached_map_vars['stage_c'][hydro_i])
            data['Depth'] = data['WSE'] - data['Bathymetry']
            if self.quantity not in ['WSE', 'Depth']:
                if self.quantity == 'Velocity':
                    xmom = np.asarray(self.cached_map_vars['xmomentum_c'][hydro_i])
                    ymom = np.asarray(self.cached_map_vars['ymomentum_c'][hydro_i])
                    with np.errstate(divide='ignore', invalid='ignore'):
                        xvel = xmom/data['Depth']
                        yvel = ymom/data['Depth']
                    data[self.quantity] = np.sqrt(xvel**2 + yvel**2)
        else:
            # if 'Mortality' in self.quantity:
            #     data[self.quantity] = self.load_DFM_ANUGA_mortality(veg_i)
            # else:
            veg_data = self._load_DFM_ANUGA_veg(ets, eco_year)
            data = self._pass_DFM_ANUGA_veg(veg_data, data)

        return data


class DFMMapLoader(BaseMapLoader):
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
                               'Velocity': 'mesh2d_ucmag', 
                               'Max Shear Stress': 'mesh2d_tausmax'}

    def _load_DFM_outputs(self, subdir):
        # only load hydro outputs if the subdirectory has changed (for DFM, it will only change in restart mode)
        if subdir != self.current_subdir:
            ncfile = Dataset(subdir / f'{self.model_name}_map.nc')
            self.cached_map_vars = ncfile.variables
            self.current_subdir = subdir


    def load(self, hydro_i, ets, eco_year):
        # get model output subdirectory and read data
        singledir = self.modeldir / 'output'
        multidir = self.modeldir / f'output_year{eco_year}_ets{ets}'
        self._load_DFM_outputs(multidir if multidir.exists() else singledir)
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

        elif not self.veg_plot:
            data['WSE'] = np.asarray(self.cached_map_vars[self.hydro_varnames['WSE']][hydro_i])
            data['Depth'] = data['WSE'] - data['Bathymetry']
            if self.quantity not in ['WSE', 'Depth']:
                data[self.quantity] = np.asarray(self.cached_map_vars[self.hydro_varnames[self.quantity]][hydro_i])
        else:
            # if 'Mortality' in self.quantity:
            #     data[self.quantity] = self._load_DFM_ANUGA_mortality(veg_i)
            # else:
            veg_data = self._load_DFM_ANUGA_veg(ets, eco_year)
            data = self._pass_DFM_ANUGA_veg(veg_data, data)

        return data