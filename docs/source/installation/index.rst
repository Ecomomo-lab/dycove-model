Installation Instructions
=========================

DYCOVE is compatible with Python 3.10-3.13 (see Issues or notes below for some current incompatibilities). 
There are three general dependencies, which are `numpy <https://numpy.org/install/>`_, `xarray <https://pypi.org/project/xarray/>`_, and `scipy <https://www.scipy.org/install.html>`_.
However, there are additional dependencies depending on the underlying numerical model to be used, and if the :class:`~dycove.utils.plotting.ModelPlotter` will be used.


Install ANUGA
-------------

For DYCOVE with ANUGA, `anuga <https://github.com/GeoscienceAustralia/anuga_core>`_ [1]_ must be installed.
ANUGA may have additional dependencies depending on the application, such as `mpi4py <https://pypi.org/project/mpi4py/>`_ for running in parallel.
The easiest way to install ANUGA is to create a new ``conda`` environment and install via ``conda-forge`` (from the `ANUGA documentation page <https://anuga.readthedocs.io/en/latest/installation/install_anuga.html>`_):

.. code-block:: python

   conda install -c conda-forge anuga mpi4py

**Update**: ANUGA is not currently importing with Python 3.12 and is causing Python to crash without an error message.
For now, use Python 3.10, 3.11, or 3.13, and let us know if it starts working again with 3.12 by responding to the relevant Issue on GitHub.


Install Delft3D FM
------------------

For DYCOVE with Delft3D FM, the `software itself <https://oss.deltares.nl/web/delft3dfm>`_ [2]_ will of course need to be licensed and installed locally.
From a Python standpoint, it is also recommended to create a new environment to install DYCOVE and its dependency `bmi-python <https://github.com/openearth/bmi-python>`_, which provides a Basic Model Interface [3]_ between Python and the Delft3D FM software.
From a dedicated environment, the easiest way to install ``bmi-python`` is from source:

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
The additional dependencies are `matplotlib <https://matplotlib.org/3.2.2/users/installing.html>`_, `tqdm <https://pypi.org/project/tqdm/>`_, and `imageio <https://pypi.org/project/ImageIO/>`_.
Of course, users will want to have some of these libraries installed whether or not they plan to use :class:`~dycove.utils.plotting.ModelPlotter`.


Installing DYCOVE
-----------------

DYCOVE itself can be installed from source like any other package:

.. code-block:: python

   git clone https://github.com/Ecomomo-lab/dycove-model.git
   cd dycove-model
   pip install -e .


References
----------

.. [1] Lesser, G. R., Roelvink, J. A., van Kester, J. A. T. M., and Stelling, G. S., 2004. Development and validation of a three-dimensional morphological model. Coastal engineering, 51(8-9), 883–915. https://doi.org/10.1016/j.coastaleng.2004.07.014

.. [2] Davies, G., and Roberts, S., 2015. Open source flood simulation with a 2D discontinuous-elevation hydrodynamic model. In Proceedings of MODSIM 2015. https://doi.org/10.1016/j.coastaleng.2011.03.010

.. [3] Hutton, E. W. H., Piper, M. D., and Tucker, G. E., 2020. The Basic Model Interface 2.0: A standard interface for coupling numerical models in the geosciences. Journal of Open Source Software, 5(51), 2317. https://doi.org/10.21105/joss.02317.
