===============================
Aerosols scattering model
===============================
|Build| |Coverage| |PyPI| |Status| |Version| |Python| |License|

.. |Build| image:: https://travis-ci.org/seignovert/python-aerosols-scattering.svg?branch=master
        :target: https://travis-ci.org/seignovert/python-aerosols-scattering
.. |Coverage| image:: https://coveralls.io/repos/github/seignovert/python-aerosols-scattering/badge.svg?branch=master
        :target: https://coveralls.io/github/seignovert/python-aerosols-scattering?branch=master
.. |PyPI| image:: https://img.shields.io/badge/PyPI-aerosols--scattering-blue.svg
        :target: https://pypi.python.org/project/aerosols-scattering
.. |Status| image:: https://img.shields.io/pypi/status/aerosols-scattering.svg?label=Status
.. |Version| image:: https://img.shields.io/pypi/v/aerosols-scattering.svg?label=Version
.. |Python| image:: https://img.shields.io/pypi/pyversions/aerosols-scattering.svg?label=Python
.. |License| image:: https://img.shields.io/pypi/l/aerosols-scattering.svg?label=License

Python package for Titan's aerosols scattering

Install
-------
With ``pip``:

.. code:: bash

    $ pip install aerosols-scattering

Or the ``source files``:

.. code:: bash

    $ git clone https://github.com/seignovert/python-aerosols-scattering.git
    $ cd aerosols-scattering ; python setup.py install

Usage
------

.. code:: python

    >>> from aerosols import mie
    >>> from aerosols import fractals

Note
----
This package is an early attempt to model Titan's aerosols scattering based on Tomasko et al. 2008 paper (doi:`10.1016/j.pss.2007.11.019`_)

.. _`10.1016/j.pss.2007.11.019`: https://dx.doi.org/10.1016/j.pss.2007.11.019