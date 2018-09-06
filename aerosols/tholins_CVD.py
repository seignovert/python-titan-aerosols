#!/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import sqlite3 as sqlite
import numpy as np

DB_NAME = 'optical_index.db'
TABLE = 'Tholins_CVD'

class Database(object):
    """Optical indexes constant database""""
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
            raise IOError(f"Database not found: {self.path}")

        self.con = sqlite.connect(self.path)
        self.db = self.con.cursor()
    
    def set_table(self, table):
        """Check table exists in the database"""
        self.db.execute(f"SELECT name FROM sqlite_master WHERE type='table' AND name='{table}'")
        if not self.db.fetchone():
            raise AttributeError(f"Table {table} not found in the database")
        self.table = table

    def fetchone(self):
        out = self.db.fetchone()
        if out is None:
            raise TypeError
        else:
            return out

    def execute(self, cmd):
        return self.db.execute(cmd)


def tholins_CVD(wvln, db=Database()):
        """
        Get laboratory produced optical index constant for tholins aerosols at
        a specific wavelength.
        
        Input
        ------
            wvln (float) Wavelength in meter

        Output
        -------
            nr (float) Real part of the optical index
            ni (float) Imaginary part of the optical index

        Note
        ----
            The values are extrapolated as constant for wavelengths above 314 um and
            below 20 nm.
            In between,
                - the real part is linear interpolated
                - the imaginary part is LOG interpolated
            Between 935 nm and 1.5 um, the bump of the imaginary part is removed and
            fixed at 7.19e-3.
        """

        # Convert wavelength meters in micrometers
        wvln *= 1.e6

        # Values SUP
        db.execute(f"SELECT wvln, nr, ni FROM {db.table} WHERE wvln >= {wvln} ORDER BY wvln ASC LIMIT 1")
        try:
            wvln_sup, nr_sup, ni_sup = db.fetchone()
        except TypeError:
            db.execute(f"SELECT MAX(wvln), nr, ni FROM {db.table}")
            wvln_max, nr, ni = db.fetchone()
            print(f">>WARNING: wvln = {wvln:.3e} m > wvln_max_db = {wvln_max*1.e-6:.3e} m => extrapolation cst")
            return nr, ni

        # Values INF
        db.execute(f"SELECT wvln, nr, ni FROM {db.table} WHERE wvln <= {wvln} ORDER BY wvln DESC LIMIT 1")
        try:
            wvln_inf, nr_inf, ni_inf = db.fetchone()
        except TypeError:
            db.execute(f"SELECT MIN(wvln), nr, ni FROM {db.table}")
            wvln_min, nr, ni = db.fetchone()
            print(f">>WARNING: wvln = {wvln:.3e} m < wvln_min_db = {wvln_min*1.e-6:.3e} m => extrapolation cst")
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

        # Remove the bump @ 1 um
        if wvln >= .935 and wvln <= 1.5:
            ni = 7.19e-3
        return nr, ni
