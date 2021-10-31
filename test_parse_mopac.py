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

pprint(p.provenance)
pprint(p._parse_cache)
print(p.atoms)
pprint(p.initial_coords)

try:
    pprint(p.opt_energies)
    pprint(p.opt_grads)
    #pprint(p.next_job_input)
    print(p.geoms.shape)
except (NameError,AttributeError):
    pass
