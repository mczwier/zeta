#!/usr/bin/env python

import sys
import os

import yaml

from pprint import pprint

from zeta.mopac import MopacParser

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('outfile', help='MOPAC primary output')
parser.add_argument('-A', '--auxfile', help='MOPAC AUX file [AUX]')
parser.add_argument('-An', '--auxnfile', help='MOPAC detailed auxfile [AUX(n,...)], '
                                              'may be stdout/err for some nn on some platforms')
args = parser.parse_args()

outfile = open(args.outfile, 'rt')
auxfile = open(args.auxfile, 'rt') if args.auxfile else None
auxnfile = open(args.auxnfile, 'rt') if args.auxnfile else None

p = MopacParser(outfile, auxfile, auxnfile)
p.parse()

print('MOPAC parser processed {:d} geometries'.format(len(p.calculation.geometries)))
try:
    pprint(p.geoms[:5])
    pprint(p.geoms[-5:])
except (NameError,AttributeError):
    pass

print('Final geometry:')
print(yaml.dump(p.calculation.geometries[-1]))
