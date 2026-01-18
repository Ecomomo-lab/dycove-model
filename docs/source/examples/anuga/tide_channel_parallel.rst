.. _tide-channel-anuga-pll:

Example - Simple beach and tide channel (ANUGA in parallel)
===========================================================

This example corresponds to the Python file `tide_channel_ANUGA_pll.py` located in the `examples/ANUGA/tide_channel_parallel/` directory, along with supporting input files.
This example uses the same bathymetry and vegetation inputs as the :doc:`serial ANUGA example <../anuga/tide_channel>`.

DYCOVE-ANUGA can be run in parallel as long as `mpi4py` is installed (see :doc:`Installation Instructions <../../installation/index>`) and the user has access to more than one processor.
DYCOVE commands remain the same, but the ANUGA domain must be generated in parallel using a slightly modified approach.
For this example, the parallel domain generation approach is provided in `gen_anuga_domain_pll.py` in the working directory. 
Users with existing ANUGA model codes that run in parallel should be able to seemlessly integrate DYCOVE in those codes by passing their parallelized domain in the same way as the serial examples:

.. code-block:: python

   HydroModel = ANUGA_hydro.ANUGA(domain, vegetation=veg_1)

This example uses a Windows batch file `tide_channel_ANUGA_pll.bat` in the working directory, which contains instructions on how to run it.
As with any batch file, view and edit it by alt-clicking and choosing the "Edit" option.
For users with HPC access, the DYCOVE-ANUGA Python script can be run via a Slurm scheduling script.

Note that for very large ANUGA domains, the parallel file merging that needs to occur at the end of the simulation can take a long time (we are working on a better solution).
We tested a simulation of 370,000 elements running on 180 cores that ran for about 2:15 hours, and the final merging took about one hour.
This time may depend on the frequency of simulation outputs, load balancing (of vegetation), and number of vegetation cohorts (colonization events).
If running on HPC, request a couple additional simulation hours just in case, until you get a sense of how long it takes at the end.