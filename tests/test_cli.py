"""Test CLI module."""

from aerosols.cli import cli_fractal_tholins

def test_cli_fractal_tholins(capsys):
    argv = '338e-9 60e-9 266'.split()
    cli_fractal_tholins(argv)
    out, err = capsys.readouterr()

    stdout = (
        '# Cross sections:\n'
        'Scattering: 2.704e-12 m^-2\n'
        'Absorption: 1.585e-12 m^-2\n'
        'Extinction: 4.289e-12 m^-2\n'
    )

    assert out == stdout
    assert err == ''

def test_cli_fractal_tholins_phase_function(capsys):
    argv = '-p --nang 3 338e-9 60e-9 266'.split()
    cli_fractal_tholins(argv)
    out, err = capsys.readouterr()

    stdout = (
        '# Phase function\n'
        '0.0\t5.85e+02\n'
        '45.0\t2.48e+00\n'
        '90.0\t5.40e-01\n'
        '135.0\t3.55e-01\n'
        '180.0\t3.67e-01\n'
    )

    assert out == stdout
    assert err == ''


def test_cli_fractal_tholins_value_err(capsys):
    argv = '-df 2.3 338e-9 60e-9 266'.split()
    cli_fractal_tholins(argv)
    out, err = capsys.readouterr()

    stdout = "Model tested only for Df = 2 (received Df=2.30)\n"

    assert out == stdout
    assert err == ''
