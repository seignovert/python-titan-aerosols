#!/bin/python
##> Compute optical index constants based on laboratory experience 
## @details   Contains only a subroutine that computes (from sql tables)
##            refractive indexes for Tholin, CH4 and C2H6.
## @author    B. Seignovert (univ-reims@seignovert.fr)
## @version   1.0
## @date      2016/07/15
# -*- coding: utf-8 -*-

# --------
# Imports
# --------
import sqlite3 as sqlite
import numpy as np
import os,sys

class OpticalIndex:
    def __init__(self,root,name):
        self.name      = name
        self.root      = root
        self.db        = root + name
        self.connexion = sqlite.connect(self.db)
        self.cursor    = self.connexion.cursor()
        return

    def __repr__(self):
        return "Tholins Database Object. Location: %s" % self.db

    def get(self,wvln,table='Tholins_CVD'): # wvln in [m]
        """ @details The method computes the optical constants of laboratory produced
                     for @bti{waveln}. Above (\>314um) and below (\<20um) the data, values 
                     are extrapolated as constant. """

        # Check if table exists
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='%s'" % table);
        if not self.cursor.fetchone(): raise AttributeError("Index table unknown in db (received table=%s)" % table)
        
        # Values SUP
        self.cursor.execute("SELECT wvln,nr,ni FROM %s WHERE wvln >= %f ORDER BY wvln ASC  LIMIT 1;" % (table,wvln*1.e6) )
        try: wvln_sup,nr_sup, ni_sup = self.cursor.fetchone()
        except TypeError:
            self.cursor.execute("SELECT MAX(wvln),nr,ni FROM %s" % table )
            wvln_max,nr,ni = self.cursor.fetchone()
            print (">>WARNING: wvln = %.3e m > wvln_max_db = %.3e m => extrapolation cst" % (wvln,wvln_max*1.e-6) )
            return nr,ni

        # Values INF
        self.cursor.execute("SELECT wvln,nr,ni FROM %s WHERE wvln <= %f ORDER BY wvln DESC LIMIT 1;" % (table,wvln*1.e6) )
        try: wvln_inf,nr_inf, ni_inf = self.cursor.fetchone()
        except TypeError:
            self.cursor.execute("SELECT MIN(wvln),nr,ni FROM %s" % table )
            wvln_min,nr,ni = self.cursor.fetchone()
            print (">>WARNING: wvln = %.3e m < wvln_min_db = %.3e m => extrapolation cst" % (wvln,wvln_min*1.e-6) )
            return nr,ni

        #Interpolation
        factor = (np.log(wvln*1.e6)-np.log(wvln_inf) ) / (np.log(wvln_sup)-np.log(wvln_inf))  # All interpolation in LOG wvln
        nr = nr_inf + factor*( nr_sup - nr_inf )                                              # Real part is linear interpolated
        ni = np.exp( np.log(ni_inf) + factor*(np.log(ni_sup)-np.log(ni_inf)) )                # Imaginary part is LOF interpolated

        if table=='THOLINS_CVD' and wvln >= .935 and wvln <= 1.5: ni = 7.19e-3                # Remove the bump @ 1 um for THOLINS_CVD
        return nr,ni
            
if __name__ == '__main__':
    index_db = OpticalIndex('','optical_index.db')

    if (len(sys.argv) == 2):
        print index_db.get( float(sys.argv[1]) )
        
    elif (len(sys.argv) == 3):
        print index_db.get( float(sys.argv[1]), sys.argv[2] )
 
