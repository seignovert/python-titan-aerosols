"""Test tholins database module."""
# pylint: disable=missing-function-docstring

from pytest import approx, raises

from aerosols.tholins import (
    Database, fractals_tholins, index_tholins, mie_tholins
)


def test_database_not_exists():
    with raises(FileNotFoundError):
        Database('wrong.db')


def test_table_not_exists():
    with raises(ValueError):
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
    qsct, _, qabs, gg, _, P = mie_tholins(300e-9, 50e-9)
    assert qsct == approx(3.256812556887943e-15, 1e-6)
    assert qabs == approx(6.0880464818229446e-15, 1e-6)
    assert gg == approx(0.25734589112794115, 1e-6)
    assert P[0] == approx(2.5061596742176055, 1e-6)
    assert P[-1] == approx(0.7083061257489704, 1e-6)


def test_fractals():
    """Test based values against Seignovert et al. (2017) - Tab. 3.

    https://doi.org/10.1016/j.icarus.2015.09.039

    The values were adjusted after d21613b fix.

    """
    qsct, qext, qabs, _, _, P = fractals_tholins(
        338e-9, 60e-9, 2.0, 266, db=Database(table='Tholins_CVD')
    )

    assert qsct == approx(2.9e-12, abs=0.1e-12)
    assert qext == approx(4.2e-12, abs=0.1e-12)
    assert qabs == approx(1.3e-12, abs=0.1e-12)
    assert P[0] == approx(185.4, abs=0.1)
    assert P[1] == approx(177.7, abs=0.1)
    assert P[2] == approx(156.8, abs=0.1)
    assert P[4] == approx(98.7, abs=0.1)
    assert P[6] == approx(52.5, abs=0.1)
    assert P[8] == approx(28.4, abs=0.1)
    assert P[10] == approx(17.4, abs=0.1)
    assert P[15] == approx(8.0, abs=0.1)
    assert P[20] == approx(4.7, abs=0.1)
    assert P[30] == approx(2.2, abs=0.1)
    assert P[40] == approx(1.1, abs=0.1)
    assert P[50] == approx(0.58, abs=0.01)
    assert P[60] == approx(0.36, abs=0.01)
    assert P[80] == approx(0.21, abs=0.01)
    assert P[100] == approx(0.146, abs=1e-3)
    assert P[120] == approx(0.119, abs=1e-3)
    assert P[140] == approx(0.114, abs=1e-3)
    assert P[160] == approx(0.117, abs=1e-3)
    assert P[180] == approx(0.119, abs=1e-3)
