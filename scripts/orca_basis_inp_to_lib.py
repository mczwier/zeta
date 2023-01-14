#!/usr/bin/env python3

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

from zeta.data.elements import element_symbol_to_name
from zeta.data.numbers import subshell_to_l

import argparse, sys, re

re_end = re.compile('end', re.IGNORECASE)
re_gto = re.compile('\s*NewGTO\s+(\S+)')
re_ecp = re.compile('\s*NewECP\s+(\S+)')
re_ncore = re.compile('\s*N_core\s+(\S+)', re.IGNORECASE)
re_lmax = re.compile('\s*lmax\s+(\S+)', re.IGNORECASE)
re_l_nterms = re.compile('\s*(\S+)\s+(\S+)')

parser = argparse.ArgumentParser(description='Convert a basis set in ORCA input format to library (GAMESS-US) format')
parser.add_argument('input_format', help='Basis and ECPs in ORCA input format')
parser.add_argument('--output', '-o', help='Output file (default: stdout)')
args = parser.parse_args()

outfile = open(args.output, 'wt') if args.output else sys.stdout

data_open = True
ecp_open = False

outfile.write('$DATA\n')
with open(args.input_format, 'rt') as infile:
    line = infile.readline()
    while line:
        m = re_gto.match(line)
        if m:
            if ecp_open:
                outfile.write('$END\n')
                ecp_open = False
                outfile.write('$DATA\n')
                data_open = True
            symbol = m.group(1)
            element = element_symbol_to_name(symbol).upper()
            #outfile.write('! Recognized GTO for element {} ({})\n'.format(symbol, element))
            outfile.write(element + '\n')
            nprim = 0
            line = infile.readline().strip()
            while not re_end.search(line):
                nprim += 1
                outfile.write(line + '\n')
                line = infile.readline().strip()
            #outfile.write('! {} primitives\n'.format(nprim))
            outfile.write('\n') # tailing empty card
        m = re_ecp.match(line)
        if m:
            if data_open:
                outfile.write('$END\n')
                data_open = False
                outfile.write('$ECP\n')
                ecp_open = True
            symbol = m.group(1).upper()
            line = infile.readline()
            ncore = int(re_ncore.match(line).group(1))
            line = infile.readline()
            lmax = subshell_to_l(re_lmax.match(line).group(1))
            outfile.write('{}-ECP GEN    {}    {}\n'.format(symbol, ncore, lmax))
            line = infile.readline()
            while not re_end.search(line):
                m = re_l_nterms.match(line)
                subshell = m.group(1).lower()
                nterms = int(m.group(2))
                if nterms > 0:
                    outfile.write('{:<7d}----- {}-ul potential -----\n'.format(nterms, subshell))
                    for _i in range(nterms):
                        line = infile.readline()
                        (idx,b,c,d) = line.split()
                        outfile.write('{:>19.9f} {:>2d}  {:>19.9f}\n'.format(float(c), int(d), float(b)))
                line = infile.readline()
            
            outfile.write('\n') #tailing empty card
        
        line = infile.readline() # end while

outfile.write('$END\n')    
if args.output:
    outfile.close()
