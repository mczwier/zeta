#!/usr/bin/env python


# Dirty hack for now
try:
    import zeta # @UnusedImport
except ModuleNotFoundError:
    import sys
    from pathlib import Path
    p = Path(__file__).parent.parent
    sys.path.insert(0,str(p))
    del p
finally:
    import zeta  # @UnusedImport @Reimport

from argparse import ArgumentParser
import timeit

from zeta.engines.mopac import MopacParser


parser = ArgumentParser()
parser.add_argument('outfile', help='MOPAC primary output')
parser.add_argument('-A', '--auxfile', help='MOPAC AUX file [AUX]')
parser.add_argument('-An', '--auxnfile', help='MOPAC detailed auxfile [AUX(n,...)], '
                                              'may be stdout/err for some nn on some platforms')
parser.add_argument('-o', '--output', help='Write data to OUTFILE in yaml format')
args = parser.parse_args()

outfile = open(args.outfile, 'rt')
auxfile = open(args.auxfile, 'rt') if args.auxfile else None
auxnfile = open(args.auxnfile, 'rt') if args.auxnfile else None

begin = timeit.default_timer()
p = MopacParser(outfile, auxfile, auxnfile)
p.parse()
end = timeit.default_timer()

print('MOPAC parser processed {:d} geometries in {:.6g} seconds'
      .format(len(p.calculation.geometries),
              (end-begin)))

if args.output:
    begin = timeit.default_timer()
    with open(args.output, 'wt') as outputfile:
        zeta.data.io.yaml.write_yaml(p.calculation, outputfile)
    end = timeit.default_timer()
    print('Wrote calculation graph to {} ({:.6g} s)'.format(args.output,
                                                            end-begin))
    
    
    