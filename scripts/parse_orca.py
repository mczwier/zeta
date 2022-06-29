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

from zeta.engines.orca import OrcaParser


parser = ArgumentParser()
parser.add_argument('outfile', help='ORCA primary output')
parser.add_argument('-o', '--output', help='Write data to OUTFILE in yaml format')
args = parser.parse_args()

outfile = open(args.outfile, 'rt')

begin = timeit.default_timer()
p = OrcaParser(outfile)
p.parse()
end = timeit.default_timer()

print('READING OUTPUT FILE {}'.format(p.provenance['output_file_name']))
print('LAST CHANGED: {}'.format(p.provenance['output_file_datetime']))
print('CALCULATION ENGINE: {package_name:s} {package_version:s}'.format(**p.provenance))
print(p.provenance)
try:
    input_text = p.provenance['input_file_text']
except KeyError:
    pass
else:
    print('INPUT FILE ({}):'.format(p.provenance.get('input_file_name', '')))
    print(input_text)


if args.output:
    begin = timeit.default_timer()
    with open(args.output, 'wt') as outputfile:
        zeta.data.io.yaml.write_yaml(p.calculation, outputfile)
    end = timeit.default_timer()
    print('Wrote calculation graph to {} ({:.6g} s)'.format(args.output,
                                                            end-begin))
    
    
    