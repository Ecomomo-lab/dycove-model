###############################################################
#  DFM_hydro.py
###############################################################

import re
import os
import sys
from types import ModuleType
from datetime import datetime
from pathlib import Path
import numpy as np
import xarray as xr

from dycove.sim.base import HydroSimulationBase, HydroEngineBase
from dycove.utils.simulation_reporting import Reporter
from dycove.constants import H_LIM_VELOCITY


r = Reporter()

def _import_bmi():
    """ Lazy loading of bmi wrapper to avoid import errors when DFM will not be tested/used. """
    try:
         # need to mock this due to incompatibility with updated setuptools versions (python 3.7+)
        sys.modules['pkg_resources'] = ModuleType('pkg_resources')
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


    def finalize_simulation(self):
        """
        Finalize DFM output and model resources.

        Serial DFM uses the standard DYCOVE finalization path unchanged.
        Parallel DFM additionally reconstructs global DFM and vegetation
        outputs, with MPI-safe propagation of rank-0 reconstruction errors.
        """
        if not self.engine.is_parallel():
            return super().finalize_simulation()

        r.report("Merging outputs, cleaning up, and finalizing simulation...")

        # All ranks must finish writing partition-local outputs first.
        self.engine.parallel_barrier()

        map_error = None

        if self.engine.get_rank() == 0:
            try:
                self.engine.merge_parallel_dfm_map()
            except Exception as exc:
                map_error = (
                    "Parallel DFM map reconstruction failed: "
                    f"{type(exc).__name__}: {exc}"
                )

        self.engine.raise_parallel_root_error(map_error)

        veg_error = None

        try:
            self.outputs.reconcile_vegetation_output(self.simstate)
        except Exception as exc:
            if self.engine.get_rank() == 0:
                veg_error = (
                    "Parallel vegetation-output reconstruction failed: "
                    f"{type(exc).__name__}: {exc}"
                )

        self.engine.raise_parallel_root_error(veg_error)

        self.engine.cleanup()
        r.report("Simulation complete!")


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
    - Parallel DFM execution is supported through MPI domain partitioning.
    - Parallel vegetation output is reconstructed using DFM global cell IDs
      and the native merged DFM map ordering.

    """

    def __init__(self, dfm_path, config_path, mdu_path, vegetation=None):

        self.dfm_path = Path(dfm_path)

        if os.name == "nt":
            self.dll_dirs = self.add_dll_directories(self.dfm_path)
            self.dflowfm_path = self.dfm_path / "dflowfm/bin/dflowfm.dll"
            self.dimr_path = self.dfm_path / "dimr/bin/dimr_dll.dll"
        else:
            self.dll_dirs = []
            self.dflowfm_path = self.dfm_path / "lib/libdflowfm.so"
            self.dimr_path = self.dfm_path / "lib/libdimr.so"

        # ---- Verify libraries exist ----
        if not self.dflowfm_path.exists():
            raise FileNotFoundError(
                f"D-Flow FM library not found: {self.dflowfm_path}"
            )

        if not self.dimr_path.exists():
            raise FileNotFoundError(
                f"DIMR library not found: {self.dimr_path}"
            )

        self.BMIWrapper = _import_bmi()


        self.mdu_path     = mdu_path  # location of MDU file that contains model directions/inputs
        self.model_dir    = mdu_path.parent  # model directory containing MDU and other model files
        self.config_path  = config_path  # location of config file used for running DFM using dimr

        self.veg = vegetation

        self.open_bmi_wrappers()


    def add_dll_directories(self, dfm_path):
        """
        Add DLL paths to env before calling BMI.

        Return list to retain handle and avoid accidental garbage collecting.

        Note that for versions of python 3.7 and earlier, you would need to set the env
        variables differently:

        .. code-block:: python

           os.environ['PATH'] = os.path.join(dfm_path, 'share', 'bin') + ";" +
                                os.path.join(dfm_path, 'dflowfm', 'bin') + ";" + ... )
        """
        return [os.add_dll_directory(dfm_path / Path("dflowfm/bin")),
                os.add_dll_directory(dfm_path / Path("dimr/bin")),
                os.add_dll_directory(dfm_path / Path("share/bin")),
                os.add_dll_directory(dfm_path / Path("dwaves/bin")),
                os.add_dll_directory(dfm_path / Path("esmf/scripts")),
                os.add_dll_directory(dfm_path / Path("swan/scripts")),
                ]


    def open_bmi_wrappers(self):
        """ Create BMI wrapper objects for DFM and DIMR """
        # BMI wrapper object that interacts with the actual numerical model (e.g., getting and setting variables)
        self.dflowfm = self.BMIWrapper(engine=str(self.dflowfm_path), configfile=str(self.mdu_path))
        # BMI wrapper object that handles the deployment of the model executables
        self.dimr = self.BMIWrapper(engine=str(self.dimr_path), configfile=str(self.config_path))


    def initialize(self):
        ### ----- Required for vegetation module ----- ###
        self.mdu_vars = self.get_model_inputs()
        self.morphology, self.morph_vars = self.get_morphodynamic_inputs()
        self.vegetation_file_check()

        ### ----- Required for numerical model ----- ###
        if self.is_parallel():
            from mpi4py import MPI

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()

            self.dimr.set_var(
                "useMPI",
                np.array([1], dtype=np.int32)
            )
            self.dimr.set_var(
                "myRank",
                np.array([rank], dtype=np.int32)
            )
            self.dimr.set_var(
                "numRanks",
                np.array([size], dtype=np.int32)
            )

            print(
                f"DYCOVE DFM parallel setup: rank={rank}, size={size}",
                flush=True,
            )

            comm.Barrier()

        self.dimr.initialize()

        # Load MPI cell ownership/global numbering in DFM BMI order.
        if self.is_parallel():
            self.load_partition_mapping()


    def step(self, seconds):
        self.dimr.update(seconds)


    def parallel_barrier(self):
        """Synchronize DFM MPI ranks."""
        if self.is_parallel():
            from mpi4py import MPI
            MPI.COMM_WORLD.Barrier()


    def raise_parallel_root_error(self, error):
        """
        Broadcast a rank-0 error to all MPI ranks and raise it consistently.

        This prevents non-root ranks from waiting indefinitely at a later
        synchronization point if a rank-0 output-reconstruction step fails.
        """
        if not self.is_parallel():
            if error is not None:
                raise RuntimeError(error)
            return

        from mpi4py import MPI

        error = MPI.COMM_WORLD.bcast(
            error if MPI.COMM_WORLD.Get_rank() == 0 else None,
            root=0,
        )

        if error is not None:
            raise RuntimeError(error)


    def cleanup(self):
        self.dimr.finalize()


    def get_cell_count(self):
        return int(self.dflowfm.get_var("ndxi"))  # number of non-boundary boxes, i.e. within-domain boxes


    def get_refdate(self):
        # This input has a line in the MDU file, but most other models probably don't care what the date is and we can just hardcode the date
        refdatestr = self.mdu_vars["RefDate"]
        return datetime(int(refdatestr[:4]), int(refdatestr[4:6]), int(refdatestr[6:]))


    def get_elevation(self):
        # DFM returns arrays with boundary values included, slice those out first
        n_cells = self.get_cell_count()
        return np.array(self.dflowfm.get_var("bl"))[:n_cells]


    def get_velocity_and_depth(self):
        n_cells = self.get_cell_count()
        depth = np.array(self.dflowfm.get_var("hs"))[:n_cells]
        velocity = np.array(self.dflowfm.get_var("ucmag"))[:n_cells]
        # Ignore velocities where depth is insufficient
        velocity = np.where(depth < H_LIM_VELOCITY, 0., velocity)

        return velocity, depth


    def get_vegetation(self):
        # Convert to numpy arrays because DFM returns pointers and we don't want to accidentally modify them
        stemdensity = np.array(self.dflowfm.get_var("rnveg"))
        stemdiameter = np.array(self.dflowfm.get_var("diaveg"))
        stemheight = np.array(self.dflowfm.get_var("stemheight"))
        return stemdensity, stemdiameter, stemheight


    def set_vegetation(self, stemdensity, stemdiameter, stemheight):
        self.dflowfm.set_var("rnveg", stemdensity)
        self.dflowfm.set_var("diaveg", stemdiameter)
        self.dflowfm.set_var("stemheight", stemheight)


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


    # --------------------------------------------------------
    # Some additional required, DFM-specific methods
    # --------------------------------------------------------

    def get_model_inputs(self):
        """ Read lines from MDU file into a dictionary """
        mdu_lines = self.mdu_path.read_text().splitlines()
        mdu_vars = {}
        for line in mdu_lines:
            if "=" in line:
                slist = re.split("=|#", line)
                mdu_vars[slist[0].strip()] = slist[1].strip()
        return mdu_vars

    def get_morphodynamic_inputs(self):
        """ Read lines from morphology file into a dictionary, if mprph files exist """
        # Same filename as MDU, different extension
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
            # # If morphology is off, ensure that vegetation mor variable is set to 0 (no burial/scour)
            # # I don't think this is necessary, burial/scour will be zero (nelson)
            # if self.veg is not None:
            #     self.veg.mor = 0
        return morphology, morph_vars


    def vegetation_file_check(self):
        """
        MODIFIES .mdu file if certain vegetation-related lines are not present:

        - Adds a filename next to 'ExtForceFile' if blank, creates the file too
        - Adds drag coefficient from VegetationAttributes if [veg] block is present
        - Adds appropriate Baptist model number (1 if no morph, 2 if morph)
        - Adds [veg] block if it is not present, including drag coefficient from VegetationAttributes
          and appropriate Baptist model number

        Creates empty text files for stem density, stem diameter, and stem height, if they don't
        already exist.

        The filenames are those specified in the vegetation .ext file in the model directory
        (e.g., "FlowFM_veg.ext").

        These files can be created beforehand if prior vegetation establishment is desired.

        Otherwise, blank files are required so that DFM knows to store these variables through time.
        """

        if self.veg is None:
            return

        # Read .mdu file lines
        self.mdu_lines = self.mdu_path.read_text().splitlines()

        # Track for .mdu file modification
        self.mdu_modified = False

        # All only execute if self.veg is not None
        self.add_extforcefile_to_mdu()
        self.add_veg_module_to_mdu()
        self.create_extforcefile()
        self.create_veg_xyz_files()

        if self.mdu_modified:
            self.write_modified_mdu()


    def add_extforcefile_to_mdu(self):
        """
        Add ExtForceFile to .mdu line if it's not there (and if vegetation is active).

        [external forcing]
        ExtForceFile = FlowFM.ext  # Old format for external forcings file ...
        """
        assert self.veg is not None  # for Pylance...

        if self.mdu_vars["ExtForceFile"] == "":
            self.mdu_modified = True

            # Get name of model/file based on name of "new" .ext file
            try:
                replacement = self.mdu_vars["ExtForceFileNew"].replace("_bnd", "")
            except:
                msg = ("Either the 'ExtForceFileNew' file name in the .mdu file does not end in the expected "
                       "'_bnd.ext', or there is no 'ExtForceFileNew' file defined in the .mdu file. If it was "
                       "purposeful that no boundaries were specified for this model, then this method must be "
                       "updated: please get in touch with us or create an Issue on GitHub.")
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
                    drag = self.veg.get_drag()
                    self.mdu_lines[i] = f"{slist[0]}= {drag:.1f}{' '*13}#{slist[2]}"


    def add_veg_module_to_mdu(self):
        """
        Add [veg] section to .mdu if it's not there (and if vegetation is active).

        Format:
        [veg]
        Vegetationmodelnr = 2     # 1: Baptist, 2: Baptist with morphology correction factor (lambda)
        Clveg             = 0.8   # Stem distance factor, default=0.8
        Cdveg             = 1.1   # Drag coefficient, pulled from input veg.json file
        Cbveg             = 0.7   # Stem stiffness coefficient, default=0.7
        """
        assert self.veg is not None  # for Pylance...

        veg_block_present = any(line.strip().startswith("[veg]") for line in self.mdu_lines)
        if not veg_block_present:
            drag = self.veg.get_drag()
            veg_model_num = 2 if self.veg.mor == 1 else 1
            self.mdu_modified = True

            self.mdu_lines.append("")
            self.mdu_lines.extend([
                "[veg]",
                f"Vegetationmodelnr                 = {veg_model_num}               # 1: Baptist et al. (2007) equation for calculation of vegetation roughness",
                "Clveg                             = 0.8             # Stem distance factor, default=0.8",
                f"Cdveg                             = {drag:.1f}             # Stem Cd coefficient, default=0.7",
                "Cbveg                             = 0.7             # Stem stiffness coefficient, default=0.7",
            ])


    def create_extforcefile(self):
        """ Create .ext file in the model directory if it doesn't exist """
        assert self.veg is not None  # for Pylance...

        ext_force_file = self.model_dir / self.mdu_vars["ExtForceFile"]
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
        if not ext_force_file.exists():
            with open(ext_force_file, "w") as f:
                f.write(content)
        # It may already exist if other spatially varying parameters are in use; need to append our content
        else:
            with open(ext_force_file, "r") as f:
                lines = f.read()
            # If the vegetation blocks are already there, skip this step (could have been copied from previous run dir)
            quantities = ["stemdensity", "stemdiameter", "stemheight"]
            if not all(s in lines for s in quantities):
                with open(ext_force_file, "w") as f:
                    f.write(lines)
                    f.write("\n")
                    f.write(content)


    def create_veg_xyz_files(self):
        """ Create required files for [veg] module to run, even if they are blank """
        assert self.veg is not None  # for Pylance...

        req_veg_files = ["stemdensity.xyz", "stemdiameter.xyz", "stemheight.xyz"]
        for filename in req_veg_files:
            veg_file = self.model_dir / filename
            if not veg_file.exists():
                with open(veg_file, "w") as f:
                    f.write("")


    def write_modified_mdu(self):
        """ Write modified .mdu lines back to file """
        self.mdu_path.write_text("\n".join(self.mdu_lines) + "\n")
        msg = "DFM MDU file updated and rewritten to include required inputs for vegetation module."
        r.report(msg)


    # --------------------------------------------------------
    # Parallel methods
    # --------------------------------------------------------

    def load_partition_mapping(self):
        """
        Build local DFM cell ownership/global numbering in BMI order.

        Global cell IDs are obtained directly from the D-Flow FM BMI
        variable ``iglobal_s``. Ownership information is read from this
        rank's partitioned NetCDF mesh and reordered into BMI-local order.
        """
        rank = self.get_rank()
        n_cells = self.get_cell_count()

        # Global cell IDs in the exact ordering used by DFM BMI arrays.
        bmi_globalnr = np.asarray(
            self.dflowfm.get_var("iglobal_s")
        ).reshape(-1)

        if len(bmi_globalnr) < n_cells:
            raise RuntimeError(
                f"Rank {rank}: iglobal_s length "
                f"{len(bmi_globalnr)} < ndxi {n_cells}"
            )

        self.partition_globalnr = np.asarray(
            bmi_globalnr[:n_cells],
            dtype=np.int64,
        )

        # Read this rank's partition MDU and use its NetFile entry.
        #
        # The partition mesh name is not necessarily derived from the
        # MDU filename. Use the partition MDU as the authoritative source
        # for the rank-local mesh filename.
        partition_mdu = (
            self.model_dir
            / f"{self.mdu_path.stem}_{rank:04d}.mdu"
        )

        if not partition_mdu.exists():
            raise FileNotFoundError(
                f"Partition MDU file not found: {partition_mdu}"
            )

        partition_mdu_vars = {}

        for line in partition_mdu.read_text().splitlines():
            if "=" in line:
                slist = re.split("=|#", line)

                if len(slist) >= 2:
                    partition_mdu_vars[
                        slist[0].strip().lower()
                    ] = slist[1].strip()

        if "netfile" not in partition_mdu_vars:
            raise RuntimeError(
                f"Partition MDU {partition_mdu.name} "
                "does not contain a NetFile entry."
            )

        net_file = (
            self.model_dir
            / partition_mdu_vars["netfile"]
        )

        if not net_file.exists():
            raise FileNotFoundError(
                f"Partition mesh file not found: {net_file}"
            )

        ds = xr.open_dataset(
            net_file,
            decode_cf=False,
        )

        try:
            if (
                "mesh2d_netelem_globalnr" in ds
                and "mesh2d_netelem_domain" in ds
            ):
                net_globalnr = np.asarray(
                    ds["mesh2d_netelem_globalnr"].values,
                    dtype=np.int64,
                )

                net_domain = np.asarray(
                    ds["mesh2d_netelem_domain"].values,
                    dtype=np.int32,
                )

            elif (
                "iglobal_s" in ds
                and "idomain" in ds
            ):
                net_globalnr = np.asarray(
                    ds["iglobal_s"].values,
                    dtype=np.int64,
                )

                net_domain = np.asarray(
                    ds["idomain"].values,
                    dtype=np.int32,
                )

            else:
                raise RuntimeError(
                    f"Partition mesh {net_file.name} does not contain "
                    "a recognized global-cell/ownership metadata pair. "
                    "Expected either "
                    "'mesh2d_netelem_globalnr' + "
                    "'mesh2d_netelem_domain' or "
                    "'iglobal_s' + 'idomain'."
                )
        finally:
            ds.close()

        if len(net_globalnr) != n_cells:
            raise RuntimeError(
                f"Rank {rank}: partition NetCDF contains "
                f"{len(net_globalnr)} cells but ndxi={n_cells}"
            )

        owner_lookup = {
            int(gid): int(owner)
            for gid, owner
            in zip(net_globalnr, net_domain)
        }

        try:
            self.partition_domain = np.asarray(
                [
                    owner_lookup[int(gid)]
                    for gid in self.partition_globalnr
                ],
                dtype=np.int32,
            )
        except KeyError as exc:
            raise RuntimeError(
                f"Rank {rank}: BMI global cell {exc.args[0]} "
                f"not found in local partition mesh {net_file.name}"
            ) from exc

        self.partition_owned = (
            self.partition_domain == rank
        )

        print(
            f"DYCOVE DFM partition mapping rank={rank}: "
            f"local={n_cells}, "
            f"owned={np.count_nonzero(self.partition_owned)}, "
            f"ghost={np.count_nonzero(~self.partition_owned)}",
            flush=True,
        )

        # Gather each rank's BMI-order global-cell mapping and ownership
        # metadata on rank 0 for reconstruction of global vegetation output.
        from mpi4py import MPI

        comm = MPI.COMM_WORLD

        local_mapping = (
            self.partition_globalnr.copy(),
            self.partition_owned.copy(),
        )

        self.parallel_partition_maps = comm.gather(
            local_mapping,
            root=0,
        )

        if rank == 0:
            r.report(
                f"DYCOVE DFM: collected partition mappings "
                f"for {len(self.parallel_partition_maps)} ranks"
            )


    def get_rank(self):
        if self.is_parallel():
            from mpi4py import MPI
            return MPI.COMM_WORLD.Get_rank()
        return 0


    def is_parallel(self):
        try:
            from mpi4py import MPI
        except ImportError:
            return False

        return MPI.COMM_WORLD.Get_size() > 1


    def merge_parallel_dfm_map(self):
        """
        Merge rank-local DFM map files into one global map file using
        Deltares' native ``dfmoutput mapmerge`` utility.

        This runs only on rank 0. The merged map provides the authoritative
        global DFM face ordering for reconstruction of DYCOVE vegetation output.
        """
        if self.get_rank() != 0:
            return

        import subprocess

        output_dir = self.model_dir / "output"
        dfmoutput = self.dfm_path / "bin/dfmoutput"

        if not dfmoutput.exists():
            raise FileNotFoundError(
                f"DFM output utility not found: {dfmoutput}"
            )

        nprocs = len(self.parallel_partition_maps)

        map_files = [
            output_dir / f"{self.mdu_path.stem}_{rank:04d}_map.nc"
            for rank in range(nprocs)
        ]

        missing = [
            str(f)
            for f in map_files
            if not f.exists()
        ]

        if missing:
            raise FileNotFoundError(
                "Missing parallel DFM map file(s): "
                + ", ".join(missing)
            )

        merged_map = (
            output_dir
            / f"{self.mdu_path.stem}_map.nc"
        )

        cmd = [
            str(dfmoutput),
            "mapmerge",
            "--infile",
            *[str(f) for f in map_files],
            "--outfile",
            str(merged_map),
            "--force",
        ]

        r.report(
            "DYCOVE DFM: merging parallel DFM map files"
        )

        subprocess.run(
            cmd,
            check=True,
            cwd=output_dir,
        )

        if not merged_map.exists():
            raise RuntimeError(
                f"DFM map merge did not create {merged_map}"
            )

        r.report(
            f"DYCOVE DFM: merged map created: {merged_map.name}"
        )

    def merge_parallel_veg(self, OutputManager):
        """
        Merge rank-local DYCOVE vegetation files into global files.

        Uses the DFM BMI-order global-cell mappings collected during
        initialization. Only cells owned by each MPI rank are inserted
        into the reconstructed global arrays; ghost-cell copies are ignored.
        """
        if self.get_rank() != 0:
            return

        outputdir = OutputManager.veg_dir
        file_index = OutputManager.file_index

        if not hasattr(self, "parallel_partition_maps"):
            raise RuntimeError(
                "Parallel partition mappings are unavailable on rank 0."
            )

        partition_maps = self.parallel_partition_maps
        nprocs = len(partition_maps)

        # DFM global cell numbering is 1-based.
        all_owned_gids = []

        for globalnr, owned in partition_maps:
            gids = np.asarray(globalnr, dtype=np.int64)
            owned = np.asarray(owned, dtype=bool)

            all_owned_gids.append(gids[owned])

        owned_global_ids = np.concatenate(all_owned_gids)

        if len(np.unique(owned_global_ids)) != len(owned_global_ids):
            raise RuntimeError(
                "Duplicate owned global DFM cell IDs found during "
                "vegetation-output reconstruction."
            )

        n_global = len(owned_global_ids)

        # The native Deltares parallel map merge defines the authoritative
        # global DFM face ordering. DYCOVE vegetation output must follow
        # this same order so it aligns directly with FlowFM_map.nc.
        merged_map = (
            self.model_dir
            / "output"
            / f"{self.mdu_path.stem}_map.nc"
        )

        if not merged_map.exists():
            raise FileNotFoundError(
                f"Merged DFM map file not found: {merged_map}"
            )

        with xr.open_dataset(
            merged_map,
            decode_cf=False,
        ) as ds:
            if "mesh2d_flowelem_globalnr" not in ds:
                raise RuntimeError(
                    f"{merged_map.name} does not contain "
                    "mesh2d_flowelem_globalnr."
                )

            target_global_ids = np.asarray(
                ds["mesh2d_flowelem_globalnr"].values,
                dtype=np.int64,
            ).reshape(-1)

        if len(target_global_ids) != n_global:
            raise RuntimeError(
                f"Merged DFM map contains {len(target_global_ids)} cells, "
                f"but DYCOVE found {n_global} owned global cells."
            )

        if len(np.unique(target_global_ids)) != n_global:
            raise RuntimeError(
                "Merged DFM map contains duplicate global cell IDs."
            )

        if not np.array_equal(
            np.sort(target_global_ids),
            np.sort(owned_global_ids),
        ):
            raise RuntimeError(
                "Merged DFM map global IDs do not match the complete "
                "set of owned DFM cells."
            )

        global_to_output = {
            int(gid): i
            for i, gid in enumerate(target_global_ids)
        }

        r.report(
            f"DYCOVE DFM: merging vegetation output from "
            f"{nprocs} ranks onto {n_global} global cells "
            f"in merged DFM map order"
        )

        for year in file_index:
            for ets in file_index[year]:
                for fname in file_index[year][ets]:

                    merged = {}
                    merged_attrs = None
                    processor_files = []

                    for p_rank in range(nprocs):

                        c_file = (
                            outputdir
                            / f"{fname}_proc{p_rank}.nc"
                        )

                        if not c_file.exists():
                            raise FileNotFoundError(
                                f"Parallel vegetation file not found: "
                                f"{c_file}"
                            )

                        globalnr, owned = partition_maps[p_rank]

                        globalnr = np.asarray(
                            globalnr,
                            dtype=np.int64,
                        )

                        owned = np.asarray(
                            owned,
                            dtype=bool,
                        )

                        owned_local_ids = np.flatnonzero(owned)

                        # Map each owned physical DFM cell into the
                        # exact face ordering of the native merged DFM map.
                        owned_output_ids = np.asarray(
                            [
                                global_to_output[int(gid)]
                                for gid in globalnr[owned]
                            ],
                            dtype=np.int64,
                        )

                        with xr.open_dataset(c_file) as c_sub:

                            if merged_attrs is None:
                                merged_attrs = dict(c_sub.attrs)

                            for key, var in c_sub.data_vars.items():

                                values = np.asarray(
                                    var.values
                                ).reshape(-1)

                                if len(values) != len(globalnr):
                                    raise RuntimeError(
                                        f"{c_file.name}: variable "
                                        f"{key!r} has length "
                                        f"{len(values)}, expected "
                                        f"{len(globalnr)}."
                                    )

                                if key not in merged:
                                    merged[key] = np.zeros(
                                        n_global,
                                        dtype=values.dtype,
                                    )

                                merged[key][
                                    owned_output_ids
                                ] = values[
                                    owned_local_ids
                                ]

                        processor_files.append(c_file)

                    OutputManager.save_netcdf(
                        outputdir,
                        fname,
                        merged,
                        saved_attrs=merged_attrs,
                    )

                    for c_file in processor_files:
                        c_file.unlink()

                    r.report(
                        f"DYCOVE DFM: merged {fname}.nc"
                    )
