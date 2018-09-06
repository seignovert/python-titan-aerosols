# -*- coding: utf-8 -*-
import pytest
import numpy as np

from aerosols.fractals import fractals, fractals_tomasko_2008

wvln = 338e-9
nr = 1.64
ni = 0.17
rm = 60e-9
Df = 2.0
N = 266
Xm = 2. * np.pi * rm / wvln

def test_fracts():
    qsct, qext, qabs, gg, theta, P = fractals(wvln, nr, ni, rm, Df, N)
    assert qsct == 2.9318512910130787e-12
    assert qabs == 1.28328626290106e-12
    assert abs(qabs + qsct - qext) < 0.001e-12
    assert gg is None
    assert theta[0] == 0
    assert theta[-1] == np.pi
    assert len(theta) == 181
    assert abs(P[0] - 185.8) < 0.1
    assert abs(P[-1] - 0.117) < 0.001
    assert len(P) == 181

def test_small_agg():
    Qs, Qa, Qe, _, _, _, _, _, _ = fractals_tomasko_2008(Df, 2, Xm, nr, ni)
    assert Qs == 0.6914269156771599
    assert Qa == 0.7118033316751097
    assert Qe == 1.4032302473522695

def test_df_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(1.9, N, Xm, nr, ni)

def test_n_min_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(Df, 1, Xm, nr, ni)

def test_n_max_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(Df, 1025, Xm, nr, ni)

def test_xm_min_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(Df, N, .9e-4, nr, ni)

def test_xm_max_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(Df, N, 1.6, nr, ni)

def test_nr_min_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(Df, N, Xm, 1.2, ni)

def test_nr_max_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(Df, N, Xm, 2.1, ni)

def test_ni_min_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(Df, N, Xm, nr, -.1)

def test_ni_max_err():
    with pytest.raises(ValueError):
        fractals_tomasko_2008(Df, N, Xm, nr, .8)
