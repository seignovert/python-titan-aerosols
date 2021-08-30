Titan aerosols models
=====================

|Build| |Coverage| |PyPI| |Status| |Version| |Python| |License| |Citation|

.. |Build| image:: https://travis-ci.org/seignovert/python-titan-aerosols.svg?branch=master
        :target: https://travis-ci.org/seignovert/python-titan-aerosols
.. |Coverage| image:: https://coveralls.io/repos/github/seignovert/python-titan-aerosols/badge.svg?branch=master
        :target: https://coveralls.io/github/seignovert/python-titan-aerosols?branch=master
.. |PyPI| image:: https://img.shields.io/badge/PyPI-aerosols--scattering-blue.svg
        :target: https://pypi.org/project/titan-aerosols/
.. |Status| image:: https://img.shields.io/pypi/status/titan-aerosols.svg?label=Status
        :target: https://pypi.org/project/titan-aerosols/
.. |Version| image:: https://img.shields.io/pypi/v/titan-aerosols.svg?label=Version
        :target: https://pypi.org/project/titan-aerosols/
.. |Python| image:: https://img.shields.io/pypi/pyversions/titan-aerosols.svg?label=Python
        :target: https://pypi.org/project/titan-aerosols/
.. |License| image:: https://img.shields.io/pypi/l/titan-aerosols.svg?label=License
        :target: https://pypi.org/project/titan-aerosols/
.. |Citation| image:: https://zenodo.org/badge/147735627.svg
        :target: https://zenodo.org/badge/latestdoi/147735627

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

Python usage
-------------

.. code:: python

    >>> from aerosols import index_tholins

    >>> nr, ni = index_tholins(338e-9)
    (1.6489699384541059, 0.2392676321412895)


    >>> from aerosols import mie_tholins

    >>> wvln = 338e-9 # Wavelength (m)
    >>> rm = 50e-9 # Monomer radius (m)

    >>> qsct, qext, qabs, gg, theta, P = mie_tholins(wvln, rm)
    (2.150748326506086e-15,
     6.519732093912762e-15,
     4.368983767406676e-15,
     0.19301947916187234,
     array([0., 0.01745329, ..., 3.14159265]),
     array([2.23653193, ..., 0.88785229]))


    >>> from aerosols import fractals_tholins

    >>> Df = 2.0
    >>> N = 266

    >>> qsct, qext, qabs, gg, theta, P = fractals_tholins(wvln, rm, Df, N)
    (1.5986535423863113e-12,
     2.5652821769307767e-12,
     9.666286345444654e-13,
     None,
     array([0, ..., 3.14159265]),
     array([135.83547352468324, ..., 0.16033083012643]))

    >>> N = 3000

    >>> qsct, qext, qabs, gg, theta, P = fractals_tholins(wvln, rm, Df, N)
    ValueError: Model tested only for N = 2 - 1024 (received N=3000)

    >>> qsct, qext, qabs, gg, theta, P = fractals_tholins(wvln, rm, Df, N, force=True)
    (1.877008401099561e-11,
     2.829777018602765e-11,
     9.527686175032043e-12,
     None,
     array([0, ..., 3.14159265]),
     array([1.20358413e+03, ..., 1.27914327e-01]))

    >>> from aerosols import mie

    >>> qsct, qext, qabs, gg, theta, P = mie(wvln, nr, ni, rm)
    (...)

    >>> from aerosols import fractals

    >>> qsct, qext, qabs, gg, theta, P = fractals(wvln, nr, ni, rm, Df, N)
    (...)

A static notebook is also available `here <https://nbviewer.jupyter.org/github/seignovert/python-titan-aerosols/blob/main/examples/Tholins_examples.ipynb>`_.


CLI usage
----------

.. code:: bash

    $ fractal_tholins --help
    usage: fractal_tholins [-h] [--phase-function] [--nang NANG]
                       [--fractal-dimension FRACTAL_DIMENSION] [--force]
                       wvln rm N
    Fractals cross-sections and phase function for tholin aggregate. Use default
    tholins indexes (CVD) and Tomasko et al. 2008.

    positional arguments:
    wvln                  Wavelength (m)
    rm                    Monomer radius (m)
    N                     Number of monomers

    optional arguments:
    -h, --help            show this help message and exit
    --phase-function, -p  Display the phase function
    --nang NANG           Number of angles for the phase function (0 -> pi/2)
    --fractal-dimension FRACTAL_DIMENSION, -df FRACTAL_DIMENSION
                            Fractal dimension
    --force, -f           Bypass validity checks


    $ fractal_tholins 338e-9 60e-9 266
    # Cross sections:
    Scattering: 2.715e-12 m^-2
    Absorption: 1.558e-12 m^-2
    Extinction: 4.273e-12 m^-2


    $ fractal_tholins -p 338e-9 60e-9 266
    # Phase function
    0.0     1.86e+02
    1.0     1.78e+02
    ...
    179.0   1.15e-01
    180.0   1.15e-01


    $ fractal_tholins -p --nang 10 338e-9 60e-9 266
    # Phase function
    0.0     2.52e+02
    10.0    2.37e+01
    ...
    170.0   1.55e-01
    180.0   1.56e-01


    $ fractal_tholins 338e-9 60e-9 266 -df 2.3
    Model tested only for Df = 2 (received Df=2.30)


    $ fractal_tholins 338e-9 60e-9 266 -df 2.3 --force
    # Cross sections:
    Scattering: 2.657e-12 m^-2
    Absorption: 1.351e-12 m^-2
    Extinction: 4.008e-12 m^-2


Note
----
This package is an early attempt to model Titan's aerosols scattering based on Tomasko et al. 2008 paper (doi:`10.1016/j.pss.2007.11.019`_)

.. _`10.1016/j.pss.2007.11.019`: https://dx.doi.org/10.1016/j.pss.2007.11.019
