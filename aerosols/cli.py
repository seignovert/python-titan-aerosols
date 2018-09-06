# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np

from .mie import NANG
from .tholins import fractals_tholins

def cli_fractal_tholins(argv=None):
    parser = argparse.ArgumentParser(description="""Fractals cross-sections and phase function for tholin aggregate.
Use default tholins indexes (CVD) and Tomasko et al. 2008.""")

    parser.add_argument('wvln', type=float, help='Wavelength (m)')
    parser.add_argument('rm', type=float, help='Monomer radius (m)')
    parser.add_argument('N', type=int, help='Number of monomers')
    parser.add_argument('--phase-function', '-p', action='store_true', 
                        help='Display the phase function')
    parser.add_argument('--nang', type=int, default=NANG,
                        help='Number of angles for the phase function (0 -> pi/2)')
    parser.add_argument('--fractal-dimension', '-df', type=float, default=2.0,
                        help='Fractal dimension')
    parser.add_argument('--force', '-f', action='store_true',
                        help='Bypass validity checks')

    args, others = parser.parse_known_args(argv)

    try:
        qsct, qext, qabs, gg, theta, P = fractals_tholins(
            args.wvln, args.rm, args.fractal_dimension, args.N, args.nang, args.force)
    except ValueError as e:
        print(e)
        return
    
    if not args.phase_function:
        print(f"# Cross sections:")
        print(f"Scattering: {qsct:.3e} m^-2")
        print(f"Absorption: {qabs:.3e} m^-2")
        print(f"Extinction: {qext:.3e} m^-2")
    else:
        print(f"# Phase function")
        for t,p in zip(np.degrees(theta), P):
            print(f"{t:.1f}\t{p:.2e}")

