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
    nr, ni = index_tholins(250e-9)
    assert nr == 1.68
    assert ni == 0.391

def test_interpolation_wvln():
    nr, ni = index_tholins(338e-9)
    assert round(nr, 2) == 1.65
    assert round(ni, 2) == 0.24

def test_ni_1um():
    nr, ni = index_tholins(0.9924e-6, db=Database(table='Tholins_CVD'))
    assert nr == 1.6515
    assert ni != 6.94e-3
    assert ni == 7.19e-3

def test_wvln_inf():
    nr, ni = index_tholins(1e-9)
    assert nr == 0.92
    assert ni == 0.098
        
def test_wvln_sup():
    nr, ni = index_tholins(801e-6)
    assert nr == 1.9168
    assert ni == 0.0001
        
def test_mie():
    qsct, qext, qabs, gg, theta, P = mie_tholins(300e-9, 50e-9)
    assert qsct == pytest.approx(3.256812556887943e-15, 1e-6)
    assert qabs == pytest.approx(6.0880464818229446e-15, 1e-6)
    assert gg == pytest.approx(0.25734589112794115, 1e-6)
    assert P[0] == pytest.approx(2.5061596742176055, 1e-6)
    assert P[-1] == pytest.approx(0.7083061257489704, 1e-6)


def test_fractals():
    """Test based values against Seignovert et al. (2017) - Tab. 3.

    https://doi.org/10.1016/j.icarus.2015.09.039

    The values were adjusted after d21613b fix.

    """
    qsct, qext, qabs, _, _, P = fractals_tholins(
        338e-9, 60e-9, 2.0, 266, db=Database(table='Tholins_CVD'))
    assert qsct == pytest.approx(2.9e-12, abs=0.1e-12)
    assert qext == pytest.approx(4.2e-12, abs=0.1e-12)
    assert qabs == pytest.approx(1.3e-12, abs=0.1e-12)
    assert P[0] == pytest.approx(185.4, abs=0.1)
    assert P[1] == pytest.approx(177.7, abs=0.1)
    assert P[2] == pytest.approx(156.8, abs=0.1)
    assert P[4] == pytest.approx(98.7, abs=0.1)
    assert P[6] == pytest.approx(52.5, abs=0.1)
    assert P[8] == pytest.approx(28.4, abs=0.1)
    assert P[10] == pytest.approx(17.4, abs=0.1)
    assert P[15] == pytest.approx(8.0, abs=0.1)
    assert P[20] == pytest.approx(4.7, abs=0.1)
    assert P[30] == pytest.approx(2.2, abs=0.1)
    assert P[40] == pytest.approx(1.1, abs=0.1)
    assert P[50] == pytest.approx(0.58, abs=0.01)
    assert P[60] == pytest.approx(0.36, abs=0.01)
    assert P[80] == pytest.approx(0.21, abs=0.01)
    assert P[100] == pytest.approx(0.146, abs=1e-3)
    assert P[120] == pytest.approx(0.119, abs=1e-3)
    assert P[140] == pytest.approx(0.114, abs=1e-3)
    assert P[160] == pytest.approx(0.117, abs=1e-3)
    assert P[180] == pytest.approx(0.119, abs=1e-3)
