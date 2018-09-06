# -*- coding: utf-8 -*-
import pytest
import numpy as np

from aerosols.mie import mie, mie_bohren_huffman

def test_mie():
    qsct, qext, qabs, gg, theta, P = mie(300e-9, 0.8, 0.3, 50e-9)
    assert qsct == 7.363164550772519e-16
    assert qabs == 5.020715990826407e-15
    assert qext == qabs + qsct
    assert gg == 0.19041709245035676
    assert theta[0] == 0
    assert theta[-1] == np.pi
    assert len(theta) == 181
    assert P[0] == 2.33396491281134
    assert P[-1] == 0.8524586063042852
    assert len(P) == 181

def test_nang_sup():
    with pytest.raises(AttributeError):
        mie_bohren_huffman(1, complex(0.8, 0.3), nang=1001)

def test_nang_inf():
    with pytest.raises(AttributeError):
        mie_bohren_huffman(1, complex(0.8, 0.3), nang=1)

def test_nmx_sup():
    with pytest.raises(ValueError):
        mie_bohren_huffman(150e3, complex(0.8, 0.3))
