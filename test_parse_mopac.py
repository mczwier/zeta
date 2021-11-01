#!/usr/bin/env python

import sys
import os

from pprint import pprint

from zeta.mopac import MopacParser

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('outfile', help='MOPAC primary output')
parser.add_argument('-L', '--logfile', help='MOPAC detailed logfile')
parser.add_argument('-A', '--auxfile', help='MOPAC AUX file')
args = parser.parse_args()

outfile = open(args.outfile, 'rt')
logfile = open(args.logfile, 'rt') if args.logfile else None
auxfile = open(args.auxfile, 'rt') if args.auxfile else None

p = MopacParser(outfile, logfile, auxfile)
p.parse()

pprint(p.calculation)
pprint(p.calculation.provenance)
pprint(p.initial_geometry)

try:
    pprint(p.geoms[:5])
    pprint(p.geoms[-5:])
except (NameError,AttributeError):
    pass
