# -*- coding: utf-8 -*-
import os
import argparse

from .tholins import fractals_tholins

def fractal_tholins(argv=None):
    parser = argparse.ArgumentParser(description='')

    parser.add_argument('-a', '--all', action='store_true',
                        help='Get all the output at once (disable: --limit/--page)')
    parser.add_argument('-l', '--limit', default=10, type=int,
                        help='Output limits (default: 10)')
    parser.add_argument('-p', '--page', default=1, type=int,
                        help='Page offset (default: 1)')
    parser.add_argument('--field', help='Field key and value (sep. ` ` or `=`)', metavar='VALUE')
