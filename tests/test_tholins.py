# -*- coding: utf-8 -*-
import pytest

from aerosols.tholins import Database, index_tholins, mie_tholins, fractals_tholins

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
    nr, ni = index_tholins(338e-9)
    assert round(nr, 2) == 1.64
    assert round(ni, 2) == 0.17

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

    
def test_fractals():
    qsct, qext, qabs, gg, theta, P = fractals_tholins(338e-9, 60e-9, 2.0, 266)
    assert abs(qsct - 2.9e-12) < 0.1e-12
    assert abs(qext - 4.2e-12) < 0.1e-12
    assert abs(qabs - 1.3e-12) < 0.1e-12
    assert abs(P[0] - 185.8) < 0.1
    assert abs(P[1] - 178.1) < 0.1
    assert abs(P[2] - 157.2) < 0.1
    assert abs(P[4] - 98.9) < 0.1
    assert abs(P[6] - 52.6) < 0.1
    assert abs(P[8] - 28.4) < 0.1
    assert abs(P[10] - 17.4) < 0.1
    assert abs(P[15] - 8.0) < 0.1
    assert abs(P[20] - 4.7) < 0.1
    assert abs(P[30] - 2.2) < 0.1
    assert abs(P[40] - 1.1) < 0.1
    assert abs(P[50] - 0.58) < 0.01
    assert abs(P[60] - 0.36) < 0.01
    assert abs(P[80] - 0.21) < 0.01
    assert abs(P[100] - 0.144) < 0.001
    assert abs(P[120] - 0.117) < 0.001
    assert abs(P[140] - 0.112) < 0.001
    assert abs(P[160] - 0.115) < 0.001
    assert abs(P[180] - 0.117) < 0.001
