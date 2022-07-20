"""Tholins database module."""

import sqlite3 as sqlite
from pathlib import Path

import numpy as np

from .fractals import fractals
from .mie import mie


DEFAULT_DB = Path(__file__).parent / 'data' / 'optical_index.db'
DEFAULT_TABLE = 'Tholins_Doose'


class Database:
    """Optical indexes constant database.

    Parameters
    ----------
    fname: str or pathlib.Path, optional
        Database location.
    table: str, optional
        Tholins indexes table name.

    Raises
    ------
    FileNotFoundError
        If the database file was not found.
    ValueError
        If the tholins indexes table is not found.

    """
    def __init__(self, fname=DEFAULT_DB, table=DEFAULT_TABLE):
        self.fname = fname
        self.table = table

    def __str__(self):
        return self.fname.name

    def __repr__(self):
        return f'<{self.__class__.__name__} {self} | Table: {self.table}>'

    @property
    def fname(self):
        """Tholin database file name."""
        return self.__fname

    @fname.setter
    def fname(self, fname):
        """Database file name setter."""
        self.__fname = Path(fname)

        if not self.__fname.exists():
            raise FileNotFoundError(f"Database not found: {self}")

        self.con = sqlite.connect(self.fname)
        self.db = self.con.cursor()

    @property
    def table(self):
        """Tholins indexes table name."""
        return self.__table

    @table.setter
    def table(self, table):
        """Tholins indexes table name setter."""
        self.db.execute(
            f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table}'")

        if not self.db.fetchone():
            raise ValueError(f"Table `{table}` not found in the database")

        self.__table = table

    def fetchone(self):
        """Fetch from the database."""
        out = self.db.fetchone()

        if out is None:
            raise ValueError('Data not found.')

        return out

    def execute(self, cmd):
        """Execute SQL string in the database"""
        return self.db.execute(cmd)


def index_tholins(wvln, db=Database()):
    """Get laboratory produced optical indexes for tholins aerosols.

    Parameters
    ----------
    wvln: float
        Wavelength (m).
    db: Database, optional
        Optical index database.

    Returns
    -------
    nr: float
        Real part of the optical index.
    ni: float
        Imaginary part of the optical index.

    Note
    ----
    The values are extrapolated as constant for wavelengths above and
    below the maximum and minimum values in the database.
    In between,
        - the real part is linear interpolated
        - the imaginary part is LOG interpolated

    For `Tholins_CVD` table between 935 nm and 1.5 µm, the bump of the
    imaginary part is removed and fixed at 7.19e-3.

    """

    # Convert wavelength meters in micrometers
    wvln *= 1e6

    # Values SUP
    db.execute(
        f"SELECT wvln, nr, ni FROM {db.table} WHERE wvln >= {wvln:.4f}"
        " ORDER BY wvln ASC LIMIT 1"
    )

    try:
        wvln_sup, nr_sup, ni_sup = db.fetchone()
    except ValueError:
        db.execute(f"SELECT MAX(wvln), nr, ni FROM {db.table}")
        wvln_max, nr, ni = db.fetchone()
        print(
            f">>WARNING: wvln = {wvln:.3e} m > wvln_max_db = {wvln_max * 1.e-6:.3e} m"
            " => extrapolation cst"
        )
        return nr, ni

    # Values INF
    db.execute(
        f"SELECT wvln, nr, ni FROM {db.table} WHERE wvln <= {wvln:.4f}"
        " ORDER BY wvln DESC LIMIT 1"
    )

    try:
        wvln_inf, nr_inf, ni_inf = db.fetchone()
    except ValueError:
        db.execute(f"SELECT MIN(wvln), nr, ni FROM {db.table}")
        wvln_min, nr, ni = db.fetchone()
        print(
            f">>WARNING: wvln = {wvln:.3e} m < wvln_min_db = {wvln_min*1.e-6:.3e} m"
            " => extrapolation cst"
        )
        return nr, ni

    # Remove the bump @ 1 um (only for Tholin_CVD)
    if db.table == 'Tholins_CVD' and .935 <= wvln <= 1.5:
        ni_sup = 7.19e-3

    # Known value (no interpolation)
    if wvln_sup == wvln_inf:
        return nr_sup, ni_sup

    # Interpolation factor (in LOG wvln)
    factor = (np.log(wvln) - np.log(wvln_inf)) / (np.log(wvln_sup) - np.log(wvln_inf))
    # Real part is linear interpolated
    nr = nr_inf + factor * (nr_sup - nr_inf)
    # Imaginary part is LOG interpolated
    ni = np.exp(np.log(ni_inf) + factor * (np.log(ni_sup) - np.log(ni_inf)))

    return nr, ni


def mie_tholins(wvln, r, db=Database(), **kwargs):
    """Mie cross-sections and phase function for tholin particle.

    Use default tholins indexes and Bohren and Huffman theory.

    Parameters
    ----------
    wvln: float
        Wavelength (m).
    r: float
        Particle radius (m) (range from 0 to π/2).
    db: Database, optional
        Optical index database.
    nang: int, optional
        Number of angles for the phase function.

    Returns
    -------
    qsct: float
        Scattering cross section (m^-2).
    qext: float
        Extinction cross section (m^-2).
    qabs: float
        Absorption cross section (m^-2).
    gg: float
        Asymmetry parameter.
    theta: numpy.ndarray
        Phase function angles (radians).
    P: numpy.ndarray
        Phase function.

    """
    nr, ni = index_tholins(wvln, db)
    return mie(wvln, nr, ni, r, **kwargs)


def fractals_tholins(wvln, rm, Df, N, db=Database(), **kwargs):
    """Fractals cross-sections and phase function for tholin aggregate.

    Use default tholins indexes and Tomasko et al. 2008.

    Parameters
    ----------
    wvln: float
        Wavelength (m).
    rm: float
        Monomer radius (m).
    Df: float
        Fractal dimension.
    N: int
        Number of monomers.
    db: Database, optional
        Optical index database.
    nang: int, optional
        Number of angles for the phase function (range from 0 to π/2).
    force: bool, optional
        Bypass validity checks.

    Returns
    -------
    qsct: float
        Scattering cross section (m^-2).
    qext: float
        Extinction cross section (m^-2).
    qabs: float
        Absorption cross section (m^-2).
    gg: float
        Asymmetry parameter.
    theta: numpy.ndarray
        Phase function angles (radians).
    P: numpy.ndarray
        Phase function.

    """
    nr, ni = index_tholins(wvln, db)
    return fractals(wvln, nr, ni, rm, Df, N, **kwargs)
