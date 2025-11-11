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
                #handler.flush = sys.stdout.flush
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

