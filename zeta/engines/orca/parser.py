'''
Created on Apr 21, 2022

@author: mzwier
'''

from zeta.parser import ParseError
from zeta.parser.textparser import TextFileParser, RegexpMatch, ContainsText, whitespace_only
from zeta.data import Provenance
from zeta.data.geom import Geometry, normalize_atoms
from zeta.data.helpers import normalize_multiplicity
import datetime, pathlib

import numpy as np

class OrcaParser:
    m_first_header = ContainsText('O   R   C   A')
    re_version = RegexpMatch('Program Version (\S+)')
    m_input_file = ContainsText('INPUT FILE')
    re_input_filename = RegexpMatch(r'^\s*NAME = (.+)')
    m_input_end = ContainsText('END OF INPUT')
    re_input_line = RegexpMatch(r'^\|\s*\d+> (.*)$')
    
    def __init__(self, outfile):
        self.outfile = TextFileParser(outfile)
        
        outpath = pathlib.Path(outfile.name)
        
        self.provenance = Provenance({'package_name': 'orca',
                                      'output_file_name': outpath.absolute(),
                                      'output_file_datetime': datetime.datetime.fromtimestamp(outpath.stat().st_mtime)
                                      })
        
    def parse(self):
        self.parse_prelude()
        self.parse_input_file()
    
    def parse_prelude(self):
        self.outfile.discard_until_match(self.m_first_header)
        self.outfile.discard_until_match(self.re_version)
        self.provenance['package_version'] = self.outfile.presult[1]
        
    def parse_input_file(self):
        self.outfile.discard_until_match(self.m_input_file)
        self.outfile.discard_until_match(self.re_input_filename)
        self.provenance['input_file_name'] = self.outfile.presult[1].strip()
        
        input_lines = []
        while self.outfile.read_until_match(self.m_input_end):
            self.outfile.assertp(self.re_input_line)
            input_lines.append(self.outfile.presult[1].rstrip())
        
        self.provenance['input_file_text'] = '\n'.join(input_lines)
        