Installation Instructions
=========================

DYCOVE is compatible with Python 3.10-3.13 (see Issues for some current incompatibilities). 
There are technically only two general dependencies, which are `numpy <https://numpy.org/install/>`_ and `xarray <https://pypi.org/project/xarray/>`_.
However, there are additional dependencies depending on the underlying numerical model to be used, and if the :class:`~dycove.utils.plotting.ModelPlotter` will be used.


Install ANUGA
-------------

For DYCOVE with ANUGA, `anuga <https://github.com/GeoscienceAustralia/anuga_core>`_ must be installed.
ANUGA may have additional dependencies depending on the application, such as `mpi4py <https://pypi.org/project/mpi4py/>`_ for running in parallel.
The easiest way to install ANUGA is to create a new ``conda`` environment and install via ``conda-forge`` (from the `ANUGA documentation page <https://anuga.readthedocs.io/en/latest/installation/install_anuga.html>`_):

.. code-block:: python

   conda install -c conda-forge anuga mpi4py

**Update**: ANUGA is not currently importing with Python 3.12 and is causing Python to crash without an error message.
For now, use Python 3.10, 3.11, or 3.13, and let us know if it starts working again with 3.12 by responding to the relevant Issue on GitHub.


Install Delft3D FM
------------------

For DYCOVE with Delft3D FM, the `software itself <https://oss.deltares.nl/web/delft3dfm>`_ will of course need to be licensed and installed locally.
From a Python standpoint, it is also recommended to create a new environment to install DYCOVE and its dependency `bmi <https://github.com/openearth/bmi-python>`_, which provides the interface between Python and the Delft3D FM software.
From a dedicated environment, the easiest way to install ``bmi`` is from source:

.. code-block:: python

   git clone https://github.com/openearth/bmi-python.git
   cd bmi-python
   pip install -e .

Windows users may also need to install `pypiwin32 <https://pypi.org/project/pypiwin32/>`_:

.. code-block:: python

   pip install pypiwin32

**Update**: The DYCOVE implementation of the DFM ``bmi`` library is not currently working with Python 3.11 and 3.12.


Install Plotting Libraries
--------------------------

DYCOVE's :class:`~dycove.utils.plotting.ModelPlotter` is a convenient way to plot outputs.
The additional dependencies are `matplotlib <https://matplotlib.org/3.2.2/users/installing.html>`_, `scipy <https://www.scipy.org/install.html>`_,  `tqdm <https://pypi.org/project/tqdm/>`_, and `imageio <https://pypi.org/project/ImageIO/>`_.
Of course, users will want to have some of these libraries installed whether or not they plan to use :class:`~dycove.utils.plotting.ModelPlotter`.


Installing DYCOVE
-----------------

DYCOVE itself can be installed from source like any other package:

.. code-block:: python

   git clone https://github.com/Ecomomo-lab/dycove-model.git
   cd dycove-model
   pip install -e .