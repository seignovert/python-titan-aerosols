# -*- coding: utf-8 -*-
import pytest

from aerosols.cli import cli_fractal_tholins

def test_cli_fractal_tholins(capsys):
    argv = '338e-9 60e-9 266'.split()
    cli_fractal_tholins(argv)
    out, err = capsys.readouterr()

    stdout = """# Cross sections:
Scattering: 2.715e-12 m^-2
Absorption: 1.558e-12 m^-2
Extinction: 4.273e-12 m^-2
"""

    assert out == stdout
    assert err == ''

def test_cli_fractal_tholins_phase_function(capsys):
    argv = '-p --nang 3 338e-9 60e-9 266'.split()
    cli_fractal_tholins(argv)
    out, err = capsys.readouterr()

    stdout = """# Phase function
0.0\t5.88e+02
45.0\t2.49e+00
90.0\t5.37e-01
135.0\t3.51e-01
180.0\t3.63e-01
"""

    assert out == stdout
    assert err == ''


def test_cli_fractal_tholins_value_err(capsys):
    argv = '-df 2.3 338e-9 60e-9 266'.split()
    cli_fractal_tholins(argv)
    out, err = capsys.readouterr()

    stdout = "Model tested only for Df = 2 (received Df=2.30)\n"

    assert out == stdout
    assert err == ''


