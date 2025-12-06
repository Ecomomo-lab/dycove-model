###############################################################
#  DFM_hydro.py
###############################################################

import re
import os
from datetime import datetime
from pathlib import Path
import numpy as np

from dycove.sim.base import HydroSimulationBase, HydroEngineBase
from dycove.utils.simulation_reporting import Reporter

r = Reporter()

def _import_bmi():
    try:
        from bmi.wrapper import BMIWrapper
        return BMIWrapper
    except ImportError:
        msg = ("The `bmi` package is not installed. "
               "Refer to the documentation for installation instructions.")
        r.report(msg, level="ERROR")
        raise ImportError(msg)
    

class DFM(HydroSimulationBase):
    """
    Hydrodynamic simulation wrapper for the Delft3D FM model.

    This class connects the generic :class:`~dycove.sim.base.HydroSimulationBase`  
    to :class:`~dycove.sim.base.engines.DFM_hydro.DFMEngine`, providing a
    consistent Python interface for running D-Flow FM through its BMI and DIMR 
    interfaces.

    Notes
    -----
    - All higher-level logic that can be abstracted from the engine classes is
      handled in :class:`~dycove.sim.base.HydroSimulationBase`; all low-level model 
      interactions are delegated to 
      :class:`~dycove.sim.base.engines.DFM_hydro.DFMEngine`.
    """

    def __init__(self, dfm_path, config_path, mdu_path, vegetation=None):

        # build DFM engine
        engine = DFMEngine(dfm_path, config_path, mdu_path, vegetation)
        # pass DFM engine to the base class
        super().__init__(engine)
    

class DFMEngine(HydroEngineBase):
    """
    Engine interface for DFM hydro-morphodynamic model.

    This engine:
    
    - Loads and initializes DFM executables (DIMR + D-Flow FM BMI).
    - Manages exchange of hydrodynamic and vegetation state variables though
      DFM-specific ``BMI-python`` wrapper.
    - Ensures that required input files are present and are consistent 
      with simulation settings.

    Parameters
    ----------
    dfm_path : Path or str
        Path to the root Delft3D-FM installation directory. Might look like this:
        'C:/Program Files (x86)/Deltares/Delft3D Flexible Mesh Suite HM (2021.03)/plugins/DeltaShell.Dimr/kernels/x64'
    config_path : Path or str
        Path to DIMR configuration file ``dimr_config.xml``
    mdu_path : Path or str
        Path to DFM MDU file.
    vegetation : VegetationSpecies or MultipleVegetationSpecies, optional
        Vegetation object passed down from the base simulation.

    Notes
    -----
    - Vegetation files (.xyz) required by DFM vegetation module are auto-created 
      if missing.
    - Parallel mode is not currently implemented.

    """

    def __init__(self, dfm_path, config_path, mdu_path, vegetation=None):
        # lazy bmi loading
        self.BMIWrapper = _import_bmi()

        # paths to Delft3D FM dll files
        self.dflowfm_path = Path(dfm_path) / ("dflowfm/bin/dflowfm.dll")
        self.dimr_path    = Path(dfm_path) / ("dimr/bin/dimr_dll.dll")
        self.mdu_path     = Path(mdu_path)  # location of MDU file that contains model directions/inputs
        self.model_dir    = mdu_path.parent  # model directory containing MDU and other model files
        self.config_path  = config_path  # location of config file used for running DFM using dimr

        # currently only need this passed here for the file checks that happens under initialize()
        # but we are keeping it and using it with all possible engines to make output writing more consistent
        #   (accessing vegetation via engine rather than passing both the the outputs class)
        self.veg = vegetation

        # add DLL paths to env before calling BMI
        # NOTE: for versions of python 3.7 and earlier, you will need to set the env variables differently:
        # os.environ['PATH'] = os.path.join(dfm_path, 'share', 'bin') + ";" + os.path.join(dfm_path, 'dflowfm', 'bin') + ...
        os.add_dll_directory(dfm_path / Path("dflowfm/bin"))
        os.add_dll_directory(dfm_path / Path("dimr/bin"))
        os.add_dll_directory(dfm_path / Path("share/bin"))
        os.add_dll_directory(dfm_path / Path("dwaves/bin"))
        os.add_dll_directory(dfm_path / Path("esmf/scripts"))
        os.add_dll_directory(dfm_path / Path("swan/scripts"))

        # BMI wrapper object that interacts with the actual numerical model (e.g., getting and setting variables)
        self.dflowfm = self.BMIWrapper(engine=str(self.dflowfm_path), configfile=str(self.mdu_path))
        # BMI wrapper object that handles the deployment of the model executables
        self.dimr = self.BMIWrapper(engine=str(self.dimr_path), configfile=str(self.config_path))

    def initialize(self):
        ### ----- Required for vegetation module ----- ###
        # get model inputs and assign to dictionaries
        self.mdu_vars = self.get_model_inputs()
        self.morphology, self.morph_vars = self.get_morphodynamic_inputs()
        # create (empty) vegetation .xyz files for each of the three vegetation variables, if they don't exist
        self.vegetation_file_check()
        ### ----- Required for numerical model ----- ###
        # initialize DFM dimr wrapper
        self.dimr.initialize()

    def step(self, seconds):
        # advance model by specified number of seconds
        self.dimr.update(seconds)

    def cleanup(self):
        # finalize/cleanup DFM dimr wrapper
        self.dimr.finalize()

    def get_rank(self):
        # TODO: implement parallel processing for DFM
        raise NotImplementedError("Parallel mode not currently implemented for Delft3D FM")

    def get_cell_count(self):
        # get number of cells in numerical model grid
        return int(self.dflowfm.get_var("ndxi"))  # number of non-boundary boxes, i.e. within-domain boxes

    def get_refdate(self):
        # dependent on input file string format
        # this input has a line in the MDU file, but most other models probably don't care what the date is and we can just hardcode the date
        refdatestr = self.mdu_vars["RefDate"]
        return datetime(int(refdatestr[:4]), int(refdatestr[4:6]), int(refdatestr[6:]))
    
    def get_elevation(self):
        # return array of elevation values at cell centers
        # DFM returns arrays with boundary values included, slice those out first 
        n_cells = self.get_cell_count()
        return np.array(self.dflowfm.get_var("bl"))[:n_cells]

    def get_velocity_and_depth(self):
        # get "accumulators" for time-varying quantities
        # the BMI/DIMR uses this roundabout way of getting time-varying quantities (pull sums of values, divide by time interval)
        is_dtint     = self.dflowfm.get_var("is_dtint")
        is_sumvalsnd = self.dflowfm.get_var("is_sumvalsnd")        
        vel_full     = np.array(is_sumvalsnd[:, 1]/is_dtint)
        depth_full   = np.array(is_sumvalsnd[:, 2]/is_dtint)
        # remove extra (boundary) cells from arrays
        n_cells = self.get_cell_count()
        velocity, depth = vel_full[:n_cells], depth_full[:n_cells]   
        # reset "accumulators" back to zero for next step
        is_sumvalsnd.fill(0.0)
        is_dtint.fill(0.0)
        return velocity, depth

    def get_vegetation(self):
        # convert to numpy arrays because DFM returns pointers and we don't want to accidentally modify them
        stemdensity = np.array(self.dflowfm.get_var("rnveg"))
        stemdiameter = np.array(self.dflowfm.get_var("diaveg"))
        stemheight = np.array(self.dflowfm.get_var("stemheight"))
        return stemdensity, stemdiameter, stemheight

    def set_vegetation(self, stemdensity, stemdiameter, stemheight):
        self.dflowfm.set_var("rnveg", stemdensity)
        self.dflowfm.set_var("diaveg", stemdiameter)
        self.dflowfm.set_var("stemheight", stemheight)


    # --------------------------------------------------------
    # Extra, required DFM methods related to input processing
    # --------------------------------------------------------

    def get_model_inputs(self):
        # Read model file
        mdu_lines = self.mdu_path.read_text().splitlines()
        mdu_vars = {}
        # Parse lines into mdu_vars
        for line in mdu_lines:
            if "=" in line:
                slist = re.split("=|#", line)
                mdu_vars[slist[0].strip()] = slist[1].strip()
        return mdu_vars

    def get_morphodynamic_inputs(self):
        # same filename as in above function, different extension
        # morphology could be turned off, so we check if there are morph files in the model directory
        mor_filepath = self.model_dir / (self.mdu_path.stem + ".mor")
        morph_vars = {}
        morphology = True
        if mor_filepath.exists():
            with open(mor_filepath) as f:
                for line in f:
                    if "=" in line:
                        slist = re.split("=|#", line)
                        morph_vars[slist[0].strip()] = slist[1].strip()
        else:
            msg = "Morphology file NOT FOUND; proceeding with morphology off."
            r.report(msg, level="WARNING")
            morphology = False
            # if morphology is off, ensure that vegetation mor variable is set to 0 (no burial/scour)
            if self.veg is not None:
                self.veg.mor = 0
        return morphology, morph_vars
    

    def vegetation_file_check(self):
        """
        MODIFIES .mdu file if certain vegetation-related lines are not present:

        - Adds a filename next to 'ExtForceFile' if blank, creates the file too
        - Adds drag coefficient from VegetationAttributes if [veg] block is present
        - Adds [veg] block if it is not present, including drag coefficient from VegetationAttributes

        Creates empty text files for stem density, stem diameter, and stem height, if they don't 
        already exist. 
        
        The filenames are those specified in the vegetation .ext file in the model directory 
        (e.g., "FlowFM_veg.ext").

        These files can be created beforehand if prior vegetation establishment is desired.
        
        Otherwise, blank files are required so that DFM knows to store these variables through time.
        """

        # Read .mdu file lines
        self.mdu_lines = self.mdu_path.read_text().splitlines()

        # Track for .mdu file modification
        self.mdu_modified = False
        
        self.add_mdu_extforcefile()
        self.add_mdu_veg_module()
        self.create_ext_file()
        self.create_veg_xyz_files()

        if self.mdu_modified:
            self.write_modified_mdu()


    def add_mdu_extforcefile(self):

        # Add ExtForceFile to .mdu if not there (and if vegetation is active)
        if self.mdu_vars["ExtForceFile"] == "" and self.veg is not None:
            self.mdu_modified = True

            drag = self.veg.attrs.drag  # get drag value from VegetationAttributes 

            # Get name of model/file based on name of "new" .ext file
            try:
                replacement = self.mdu_vars["ExtForceFileNew"].replace("_bnd", "")
            except:
                msg = ("Either the 'ExtForceFileNew' file name in the .mdu file does not end in the expected "
                       "'_bnd.ext', or there is no 'ExtForceFileNew' file defined in the .mdu file")
                r.report(msg, level="ERROR")
                raise NameError(msg)

            self.mdu_vars["ExtForceFile"] = replacement

            for i, line in enumerate(self.mdu_lines):
                # Replace blank space with name of required .ext file
                if line.startswith("ExtForceFile "):
                    slist = re.split("=|#", line)
                    n_spaces = len(slist[1])
                    self.mdu_lines[i] = f"{slist[0]}= {replacement}{' '*max(n_spaces - len(replacement) - 1, 1)}#{slist[2]}"
                # Replace drag coefficient value with the one provided in input .json file (if [veg] block is present)
                if line.strip().startswith("Cdveg"):
                    slist = re.split("=|#", line)
                    self.mdu_lines[i] = f"{slist[0]}= {drag:.1f}{' '*13}#{slist[2]}"


    def add_mdu_veg_module(self):

        # Add [veg] section to .mdu if not there (and if vegetation is active)
        veg_block_present = any(line.strip().startswith("[veg]") for line in self.mdu_lines)
        if not veg_block_present and self.veg is not None:
            self.mdu_modified = True

            drag = self.veg.attrs.drag  # get drag value from VegetationAttributes
            
            self.mdu_lines.append("")
            self.mdu_lines.extend([
                "[veg]",
                "Vegetationmodelnr                 = 1               # 1: Baptist et al. (2007) equation for calculation of vegetation roughness",
                "Clveg                             = 0.8             # Stem distance factor, default=0.8",
                f"Cdveg                             = {drag:.1f}             # Stem Cd coefficient, default=0.7",
                "Cbveg                             = 0.7             # Stem stiffness coefficient, default=0.7",
            ])
                

    def write_modified_mdu(self):
        self.mdu_path.write_text("\n".join(self.mdu_lines) + "\n")
        msg = "DFM MDU file updated and rewritten to include required inputs for vegetation module."
        r.report(msg)



    def create_ext_file(self):
        # Create file in the model directory if it doesn't exist
        ext_force_file = self.model_dir / self.mdu_vars["ExtForceFile"]
        if not ext_force_file.exists():
            content = """QUANTITY=stemdensity
FILENAME=stemdensity.xyz
FILETYPE=7
METHOD=5
OPERAND=O

QUANTITY=stemdiameter
FILENAME=stemdiameter.xyz
FILETYPE=7
METHOD=5
OPERAND=O

QUANTITY=stemheight
FILENAME=stemheight.xyz
FILETYPE=7
METHOD=5
OPERAND=O
"""
            with open(ext_force_file, "w") as f:
                f.write(content)


    def create_veg_xyz_files(self):
        # get file names from vegetation boundary file, write blank files to model directory if they don't exist
        #ext_force_file = self.model_dir / self.mdu_vars["ExtForceFile"]
        req_veg_files = ["stemdensity.xyz", "stemdiameter.xyz", "stemheight.xyz"]
        # if not ext_force_file.exists() or self.mdu_vars["ExtForceFile"] == "":
        #     msg = ("An *.ext file for external forcing (old format) must be specified next to 'ExtForceFile' in the MDU file, "
        #             "and that file must exist in the model folder (at the same level as the MDU file)")
        #     r.report(msg, level="ERROR")
        #     raise FileNotFoundError(msg)
        # with open(ext_force_file, "r") as f:
        #     for line in f:
        #         slist = line.strip().split("=")
        #         if slist[0] == "FILENAME" and slist[1] in req_veg_files:
        #             veg_file = self.model_dir / slist[1]

        for filename in req_veg_files:
            veg_file = self.model_dir / filename
            if not veg_file.exists():
                with open(veg_file, "w") as f:
                    f.write("")


    def check_simulation_inputs(self, simstate):
        """
        Compare DYCOVE simulation time to MDU simulation time.

        All MDU files (DFM models) will have a simulation time specifed, but DYCOVE will run
        DFM for a period of time based on how many veg years we want to simulate. Basically, 
        the time specified in the MDU needs to be arbitrarily large enough so that we never 
        run into the issue of the model stopping prematurely.
        """

        if simstate.hydro_sim_days*86400 > int(self.mdu_vars["TStop"]):
            msg = (f"Model simulation time specified in MDU file (TStop = {self.mdu_vars['TStop']}) not long enough based on "
                   "the inputs provided for sim_years, n_ets, and veg_interval, which give simulation length of "
                   f"{simstate.hydro_sim_days*86400}. Please provide an arbitrarily larger number for TStop in MDU file.")
            r.report(msg, level="ERROR")
            raise ValueError(msg)

    def is_parallel(self):
        # DFM not currently set up for parallel processing in this model, but setting this up for future use
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            size = comm.Get_size()
            return True if size > 1 else False
        except:
            return False
        
    def merge_parallel_veg(self, OutputManager):
        raise NotImplementedError("Parallel mode not currently implemented for Delft3D FM")
    