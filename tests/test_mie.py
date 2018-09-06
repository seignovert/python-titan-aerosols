# -*- coding: utf-8 -*-
import pytest

from aerosols.mie import mie_tholins, mie_bohren_huffman

def test_mie_tholins():
    qsct, qext, qabs, gg, theta, P = mie_tholins(300e-9, 50e-9)
    assert qsct == 3.1862449318141566e-15
    assert qext == 7.66651167679358e-15
    assert qabs == qext - qsct
    assert gg == 0.255474286672091

def test_nang_sup():
    with pytest.raises(AttributeError):
        mie_bohren_huffman(1, complex(0.8, 0.3), nang=1001)

def test_nang_inf():
    with pytest.raises(AttributeError):
        mie_bohren_huffman(1, complex(0.8, 0.3), nang=1)

def test_nmx_sup():
    with pytest.raises(ValueError):
        mie_bohren_huffman(150e3, complex(0.8, 0.3))
