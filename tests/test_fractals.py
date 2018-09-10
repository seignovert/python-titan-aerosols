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
    assert qsct == pytest.approx(2.9318512910130787e-12, 1e-6)
    assert qabs == pytest.approx(1.28328626290106e-12, 1e-6)
    assert qabs + qsct == pytest.approx(qext, 1e-6)
    assert gg is None
    assert theta[0] == 0
    assert theta[-1] == np.pi
    assert len(theta) == 181
    assert P[0] == pytest.approx(185.8, abs=0.1)
    assert P[-1] == pytest.approx(0.117, abs=1e-3)
    assert len(P) == 181

def test_small_agg():
    Qs, Qa, Qe, _, _, _, _, _, _ = fractals_tomasko_2008(Df, 2, Xm, nr, ni)
    assert Qs == pytest.approx(0.6914269156771599, 1e-6)
    assert Qa == pytest.approx(0.7118033316751097, 1e-6)
    assert Qe == pytest.approx(1.4032302473522695, 1e-6)

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
