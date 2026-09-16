.. _tide-channel-dfm-parallel:

Example - Simple beach and tide channel (Delft3D FM in parallel)
================================================================

This example demonstrates running DYCOVE coupled with Delft3D Flexible Mesh
(DFM) using MPI domain decomposition.

The example is located in
``examples/DFM/tide_channel_parallel/`` and is based on the serial
:doc:`DFM tide-channel example <tide_channel>`. It uses the same model
geometry, boundary conditions, and vegetation configuration.

This example does not consider morphology.

Parallel DFM workflow
---------------------

Parallel DFM execution requires three settings to use the same number of
MPI processes:

1. The DFM model must be partitioned into that number of domains.
2. The DIMR ``<process>`` list must contain the corresponding process IDs.
3. DYCOVE must be launched with the same number of MPI processes.

For example, a four-process run corresponds to:

.. code-block:: text

   dflowfm --partition:ndomains=4 FlowFM.mdu
   <process>0 1 2 3</process>
   mpiexec -np 4 python run_tide_channel.py

The supplied ``run_parallel_linux.sh`` script keeps these settings
synchronized using a single ``NUM_PROCESSES`` value.

Requirements
------------

A Linux Delft3D FM installation with MPI support is required. The Python
environment must also include ``mpi4py`` and the normal DYCOVE
dependencies.

Set ``D3D_HOME`` to the root of the Linux Delft3D FM installation. The
directory is expected to contain ``bin/`` and ``lib/`` directories.

For example:

.. code-block:: bash

   export D3D_HOME=/path/to/delft3d/lnx64

Running the example
-------------------

The default number of MPI processes is two:

.. code-block:: bash

   ./run_parallel_linux.sh

A different number can be selected with ``NUM_PROCESSES``:

.. code-block:: bash

   NUM_PROCESSES=4 ./run_parallel_linux.sh

The script:

1. updates the DIMR ``<process>`` list;
2. removes old partition products;
3. partitions the DFM model using ``dflowfm --partition``; and
4. launches the DYCOVE Python driver with ``mpiexec``.

DYCOVE determines the active MPI rank and communicator size from
``MPI.COMM_WORLD``; the Python model setup itself does not require a
hard-coded rank count.

Parallel output
---------------

During a parallel DFM run, each MPI rank operates on its local DFM
partition, which contains both owned and ghost cells.

DYCOVE uses DFM global-cell numbering and ownership information to
identify the physical cells owned by each rank. Rank-local vegetation
outputs are written during the simulation.

At finalization, DYCOVE:

1. uses the native Deltares ``dfmoutput mapmerge`` utility to reconstruct
   the global DFM map;
2. uses the merged DFM map global-cell ordering as the authoritative
   physical-cell order;
3. reconstructs vegetation output from owned cells only, excluding ghost
   copies; and
4. writes merged vegetation files aligned with the global DFM map.

Notes
-----

The supplied MDU comments out the legacy ``Writebalancefile`` setting
because this keyword is unsupported by some DFM kernels.

MPI launch configuration varies between systems. On HPC systems, users
may need to place the same workflow inside the scheduler or container
configuration appropriate for their installation.
