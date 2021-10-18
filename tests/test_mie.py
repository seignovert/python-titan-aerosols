"""Test Mie module."""
# pylint: disable=missing-function-docstring

import numpy as np

from pytest import approx, raises

from aerosols.mie import mie, mie_bohren_huffman


def test_mie():
    qsct, qext, qabs, gg, theta, P = mie(300e-9, 0.8, 0.3, 50e-9)
    assert qsct == approx(7.363164550772519e-16, 1e-6)
    assert qabs == approx(5.020715990826407e-15, 1e-6)
    assert qabs + qsct == approx(qext, 1e-6)
    assert gg == approx(0.19041709245035676, 1e-6)
    assert theta[0] == 0
    assert theta[-1] == np.pi
    assert len(theta) == 181
    assert P[0] == approx(2.33396491281134, 1e-6)
    assert P[-1] == approx(0.8524586063042852, 1e-6)
    assert len(P) == 181


def test_nang_sup():
    with raises(ValueError):
        mie_bohren_huffman(1, complex(0.8, 0.3), nang=1001)


def test_nang_inf():
    with raises(ValueError):
        mie_bohren_huffman(1, complex(0.8, 0.3), nang=1)


def test_nmx_sup():
    with raises(ValueError):
        mie_bohren_huffman(150e3, complex(0.8, 0.3))
