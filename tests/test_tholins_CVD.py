# -*- coding: utf-8 -*-
import pytest

from aerosols.tholins_CVD import Database, tholins_CVD

def test_database_not_exists():
    with pytest.raises(IOError):
        Database(name='wrong.db')

def test_table_not_exists():
    with pytest.raises(AttributeError):
        Database(table='wrong_table')

def test_known_wvln():
    nr, ni = tholins_CVD(40e-9)
    assert nr == 0.8083
    assert ni == 0.389

def test_interpolation_wvln():
    nr, ni = tholins_CVD(30e-9)
    assert nr == 0.8568122639904139
    assert ni == 0.17585834352774574

def test_ni_1um():
    nr, ni = tholins_CVD(0.9924e-6)
    assert nr == 1.6515
    assert ni != 6.94e-3
    assert ni == 7.19e-3

def test_wvln_inf():
    nr, ni = tholins_CVD(1e-9)
    assert nr == 0.92
    assert ni == 0.0681
        
def test_wvln_sup():
    nr, ni = tholins_CVD(315e-6)
    assert nr == 1.9168
    assert ni == 0.0139
    
    
    
