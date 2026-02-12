import time
import logging


class Reporter:
    """ 
    Class that handles all printing of info/warning/errors to the screen.
    """

    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)

            # bootstrap logger
            try:
                from mpi4py import MPI
                comm = MPI.COMM_WORLD
                rank = comm.Get_rank()
            except ImportError:
                comm, rank = None, 0  # no MPI, single-core fallback

            logger = logging.getLogger("DynamicVegetationModel")
            logger.setLevel(logging.INFO)

            if not logger.handlers:  # avoid duplicate handlers
                handler = logging.StreamHandler()
                formatter = logging.Formatter(
                    f"%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S"
                )
                handler.setFormatter(formatter)
                logger.addHandler(handler)

            if rank != 0:  # suppress non-root ranks
                logger.setLevel(logging.CRITICAL)

            logger.propagate = False
            cls._instance.logger = logger
            cls._instance.rank = rank

        return cls._instance

    def report(self, msg, level="INFO"):
        if level == "INFO":
            self.logger.info(msg)
        elif level == "WARNING":
            self.logger.warning(msg)
        elif level == "ERROR":
            self.logger.error(msg)
        else:
            self.logger.log(level, msg)


    def print_model_time_info(self, simstate, veg_active):
        """ Print duration of simulation in hydrodynamic days and eco-morphodynamic years. """
        msg1 = (f"Running DYCOVE with {simstate.n_ets} Eco Time Steps per year with coupling "
                f"interval of {int(simstate.veg_interval/3600.)} hours")
        msg2 = f"Hydrodynamic model duration: {int(simstate.hydro_sim_days)} days"
        msg3 = f"Eco-morphodynamic model duration: {simstate.veg_sim_years} years"
        if veg_active:
            self.report(msg1)
            self.report(msg2)
            self.report(msg3)
        else:
            self.report(msg2)


    def print_runtime_updates(self, simstate, i, veg_active):
        """ Print current simulation time, elapsed wall time, and projected total wall time. """
        elapsed_units = 1., "seconds", 0      
        simstate.times_elapsed.append(time.time() - simstate.time_0)
        if simstate.times_elapsed[-1] > 60:
            elapsed_units = 60., "minutes", 1
            if simstate.times_elapsed[-1] > 3600:
                elapsed_units = 3600., "hours", 2
        t_str = round(simstate.times_elapsed[-1]/elapsed_units[0], elapsed_units[2])

        self.report(f"Current hydrodynamic model time = {round(simstate.hydrotime_seconds/86400., 1)} days")
        if veg_active:
            self.report(f"Current eco-morphodynamic model time = {round(simstate.hydrotime_seconds/86400./simstate.days_per_year*simstate.ecofac, 1)} years")
            self.report(f"Current eco-morphodynamic model date = {simstate.vegtime_date}")
        self.report(f"Total time elapsed: {t_str} {elapsed_units[1]}")        
        
        proj_units_avg = 1., "seconds", 0
        time_per_step_avg = simstate.times_elapsed[-1] if i == 0 else simstate.times_elapsed[-1] / i
        n_steps = simstate.n_veg_steps if veg_active else simstate.n_hydro_steps
        proj_time_avg = time_per_step_avg * n_steps  # (avg time per step) * (n steps total)
        if proj_time_avg > 60:
            proj_units_avg = 60., "minutes", 1
            if proj_time_avg > 3600:
                proj_units_avg = 3600., "hours", 2

        t_str = round(proj_time_avg/proj_units_avg[0], proj_units_avg[2])
        self.report(f"Projected total run time based on average yield step: {t_str} {proj_units_avg[1]}"
                    "\n\n##### --------------------------------------------- #####\n")
