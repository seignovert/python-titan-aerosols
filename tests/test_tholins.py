# -*- coding: utf-8 -*-
import pytest

from aerosols.tholins import Database, index_tholins, mie_tholins

def test_database_not_exists():
    with pytest.raises(IOError):
        Database(name='wrong.db')

def test_table_not_exists():
    with pytest.raises(AttributeError):
        Database(table='wrong_table')

def test_known_wvln():
    nr, ni = index_tholins(40e-9)
    assert nr == 0.8083
    assert ni == 0.389

def test_interpolation_wvln():
    nr, ni = index_tholins(30e-9)
    assert nr == 0.8568122639904139
    assert ni == 0.17585834352774574

def test_ni_1um():
    nr, ni = index_tholins(0.9924e-6)
    assert nr == 1.6515
    assert ni != 6.94e-3
    assert ni == 7.19e-3

def test_wvln_inf():
    nr, ni = index_tholins(1e-9)
    assert nr == 0.92
    assert ni == 0.0681
        
def test_wvln_sup():
    nr, ni = index_tholins(315e-6)
    assert nr == 1.9168
    assert ni == 0.0139
        
def test_mie():
    qsct, qext, qabs, gg, theta, P = mie_tholins(300e-9, 50e-9)
    assert qsct == 3.1862449318141566e-15
    assert qabs == 4.480266744979424e-15
    assert gg == 0.255474286672091
    assert P[0] == 2.495924768745142
    assert P[-1] == 0.7154853574567827

    
