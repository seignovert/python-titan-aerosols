# -*- coding: utf-8 -*-
import os
import sys
import sqlite3 as sqlite
import numpy as np

from .mie import NANG, mie
from .fractals import FORCE, fractals

DB_NAME = 'optical_index.db'
TABLE = 'Tholins_Doose'

class Database(object):
    """Optical indexes constant database"""
    def __init__(self, name=DB_NAME, table=TABLE):
        self.set_db(name)
        self.set_table(table)

    def set_db(self, name):
        """Check if the database exists and init the connection"""
        self.path = os.path.join(
            os.path.abspath(os.path.dirname(__file__)),
            'data',
            name
        )

        if not os.path.isfile(self.path):
            raise IOError("Database not found: {}".format(self.path))

        self.con = sqlite.connect(self.path)
        self.db = self.con.cursor()
    
    def set_table(self, table):
        """Check table exists in the database"""
        self.db.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='{}'".format(table))
        if not self.db.fetchone():
            raise AttributeError("Table {} not found in the database".format(table))
        self.table = table

    def fetchone(self):
        out = self.db.fetchone()
        if out is None:
            raise TypeError
        else:
            return out

    def execute(self, cmd):
        return self.db.execute(cmd)


def index_tholins(wvln, db=Database()):
        """
        Get laboratory produced optical index constant for tholins aerosols at
        a specific wavelength.
        
        Input
        ------
            wvln (float)    Wavelength in meter
            [db] (object)   Optical index database

        Output
        -------
            nr (float) Real part of the optical index
            ni (float) Imaginary part of the optical index

        Note
        ----
            The values are extrapolated as constant for wavelengths above and
            below the maximum and minimum values in the database.
            In between,
                - the real part is linear interpolated
                - the imaginary part is LOG interpolated
            
            For 'Tholins_CVD' table between 935 nm and 1.5 um, the bump of the
                imaginary part is removed and fixed at 7.19e-3.
        """

        # Convert wavelength meters in micrometers
        wvln *= 1.e6

        # Values SUP
        db.execute("SELECT wvln, nr, ni FROM {} WHERE wvln >= {} ORDER BY wvln ASC LIMIT 1".format(db.table, wvln))
        try:
            wvln_sup, nr_sup, ni_sup = db.fetchone()
        except TypeError:
            db.execute("SELECT MAX(wvln), nr, ni FROM {}".format(db.table))
            wvln_max, nr, ni = db.fetchone()
            print(">>WARNING: wvln = {:.3e} m > wvln_max_db = {:.3e} m => extrapolation cst".format(wvln, wvln_max*1.e-6))
            return nr, ni

        # Values INF
        db.execute("SELECT wvln, nr, ni FROM {} WHERE wvln <= {} ORDER BY wvln DESC LIMIT 1".format(db.table, wvln))
        try:
            wvln_inf, nr_inf, ni_inf = db.fetchone()
        except TypeError:
            db.execute("SELECT MIN(wvln), nr, ni FROM {}".format(db.table))
            wvln_min, nr, ni = db.fetchone()
            print(">>WARNING: wvln = {:.3e} m < wvln_min_db = {:.3e} m => extrapolation cst".format(wvln, wvln_min*1.e-6))
            return nr, ni

        # Known value (no interpolation)
        if wvln_sup == wvln_inf:
            return nr_sup, ni_inf

        # Interpolation factor (in LOG wvln)
        factor = (np.log(wvln)-np.log(wvln_inf)) / (np.log(wvln_sup) - np.log(wvln_inf))
        # Real part is linear interpolated
        nr = nr_inf + factor*(nr_sup - nr_inf)
        # Imaginary part is LOG interpolated
        ni = np.exp(np.log(ni_inf) + factor*(np.log(ni_sup)-np.log(ni_inf)))

        # Remove the bump @ 1 um (only for Tholin_CVD)
        if db.table == 'Tholins_CVD' and wvln >= .935 and wvln <= 1.5:
            ni = 7.19e-3
        return nr, ni


def mie_tholins(wvln, r, nang=NANG, db=Database()):
    """
    Mie cross-sections and phase function for tholin particule.
    Use default tholins indexes and Bohren and Huffman theory.

    Input
    ------
        wvln (float)    Wavelength (m)
        r    (float)    Particule radius (m)
        [nang] (int)    Number of angles for the phase function 
                         (range from 0 to pi/2)
        [db] (object)   Optical index database

    Outputs
    --------
        qsct  (float)    Scattering cross section (m^-2)
        qext  (float)    Extinction cross section (m^-2)
        qabs  (float)    Absorption cross section (m^-2)
        gg    (float)    Asymmetry parameter
        theta (float[])  Phase function angles (radians)
        P     (float[])  Phase function ()
    """
    nr, ni = index_tholins(wvln, db)
    return mie(wvln, nr, ni, r, nang)


def fractals_tholins(wvln, rm, Df, N, nang=NANG, force=FORCE, db=Database()):
    """
    Fractals cross-sections and phase function for tholin aggregate.
    Use default tholins indexes and Tomasko et al. 2008.

    Input
    ------
        wvln    (float) Wavelength (m)
        rm      (float) Monomer radius (m)
        Df      (float) Fractal dimension
        N       (int)   Number of monomers
        [nang]  (int)   Number of angles for the phase function
                         (range from 0 to pi/2)
        [force] (bool)  Bypass validity checks
        [db] (object)   Optical index database

    Outputs
    --------
        qsct  (float)    Scattering cross section (m^-2)
        qext  (float)    Extinction cross section (m^-2)
        qabs  (float)    Absorption cross section (m^-2)
        gg    (float)    Asymmetry parameter
        theta (float[])  Phase function angles (radians)
        P     (float[])  Phase function ()
    """
    nr, ni = index_tholins(wvln, db)
    return fractals(wvln, nr, ni, rm, Df, N, nang, force)
