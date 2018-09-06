===============================
Titan aerosols models
===============================
|Build| |Coverage| |PyPI| |Status| |Version| |Python| |License|

.. |Build| image:: https://travis-ci.org/seignovert/python-titan-aerosols.svg?branch=master
        :target: https://travis-ci.org/seignovert/python-titan-aerosols
.. |Coverage| image:: https://coveralls.io/repos/github/seignovert/python-titan-aerosols/badge.svg?branch=master
        :target: https://coveralls.io/github/seignovert/python-titan-aerosols?branch=master
.. |PyPI| image:: https://img.shields.io/badge/PyPI-aerosols--scattering-blue.svg
        :target: https://pypi.python.org/project/titan-aerosols
.. |Status| image:: https://img.shields.io/pypi/status/titan-aerosols.svg?label=Status
        :target: https://pypi.python.org/project/titan-aerosols
.. |Version| image:: https://img.shields.io/pypi/v/titan-aerosols.svg?label=Version
        :target: https://pypi.python.org/project/titan-aerosols
.. |Python| image:: https://img.shields.io/pypi/pyversions/titan-aerosols.svg?label=Python
        :target: https://pypi.python.org/project/titan-aerosols
.. |License| image:: https://img.shields.io/pypi/l/titan-aerosols.svg?label=License
        :target: https://pypi.python.org/project/titan-aerosols

Python package for Titan's aerosols models

Install
-------
With ``pip``:

.. code:: bash

    $ pip install titan-aerosols

Or directly from the ``source files``:

.. code:: bash

    $ git clone https://github.com/seignovert/python-titan-aerosols.git
    $ cd python-titan-aerosols ; python setup.py install

Usage
------

.. code:: python

    >>> from aerosols import mie
    >>> from aerosols import fractals

Note
----
This package is an early attempt to model Titan's aerosols scattering based on Tomasko et al. 2008 paper (doi:`10.1016/j.pss.2007.11.019`_)

.. _`10.1016/j.pss.2007.11.019`: https://dx.doi.org/10.1016/j.pss.2007.11.019