'''
Created on Oct 30, 2021

@author: mzwier
'''

import re
from ..parser import TextFileParser, ffloat, FailedMatchError

import numpy as np

# Much faster than numpy.loadtxt
def _read_positions_py(source, natoms):
    coords = np.empty((natoms, 3), dtype=np.float64)
    
    for iatom in range(0, natoms):
        line = source.readline()

        coords[iatom,0] = float(line[0:18])
        coords[iatom,1] = float(line[18:36])
        coords[iatom,2] = float(line[36:54])
        
    return coords
read_positions = _read_positions_py

# cython is a factor of three improvement
try:
    import pyximport
    pyximport.install(inplace=True, language_level=3)
    from ._parsext import _assign_positions_cy
    
    def _read_positions_cy(source, natoms):
        coords = np.empty((natoms,3), dtype=np.float64)
        lines = [source.readline() for _i in range(natoms) ]
        _assign_positions_cy(lines, natoms, coords)
        return coords
    read_positions = _read_positions_cy
except:
    pass

# Numba is about 100% the cost of pure python, and 5-6x cyton
#try:
# import numba
#
# @numba.jit(nopython=False)
# def _assign_positions_nb(lines, natoms, coords):
#     for iatom in range(natoms):
#         line = lines[iatom]
#         coords[iatom,0] = float(line[0:18])
#         coords[iatom,1] = float(line[18:36])
#         coords[iatom,2] = float(line[36:54])        
#
# def _read_positions_nb(source, natoms):
#     coords = np.empty((natoms,3), dtype=np.float64)
#     lines = [source.readline() for _i in range(natoms) ]
#     _assign_positions_nb(lines, natoms, coords)
#     return coords
#
# read_positions = _read_positions_nb
#except:
#    pass

class MopacParser:
    re_section_break = re.compile(r'\*\*\*')
    re_mopac_version = re.compile(r'Version: (\S+)')
    re_prelude_header = re.compile(r'CALCULATION DONE:\s+([^\*]+)')
    re_prelude_keyword = re.compile(r'\*\s+(\S+)')
    re_system_charge = re.compile(r'CHARGE ON SYSTEM = (\S+)')

    re_begin_opt = re.compile(r'Geometry optimization')    
    
    # HEAT_OF_FORM_UPDATED:KCAL/MOL=-0.57545763423410D+04
    re_energy_updated = re.compile(r'HEAT_OF_FORM_UPDATED:([^=]+)=(\S+)') 
    re_grad_updated = re.compile(r'GRADIENT_UPDATED:([^=]+)=(\S+)')
    re_geom_updated = re.compile(r'ATOM_X_UPDATED')

    re_job_time = re.compile(r'JOB TIME:\s+(\S+)')
    
    def __init__(self, outfile, logfile=None, auxfile=None):
        self.provenance = {}
        self._parse_cache = {} # intermediate information
        self.outfile = TextFileParser(outfile)
        self.logfile = TextFileParser(logfile) if logfile else None
        self.auxfile = TextFileParser(auxfile) if auxfile else None
        self.geoms = []
        self.opt_energies = []
        self.opt_grads = []
        
    def parse(self):
        self.parse_prelude()
        self.parse_initial_geometry()
    
        # Skip extra Cartesian block if present
        if self.outfile.test_lookahead(re.compile('CARTESIAN COORDINATES'), 5):
            self.outfile.discard_to_match(re.compile('1'))
            self.outfile.skip_to_blank()
    
        while self.outfile.peek():
            self.outfile.scan_and_dispatch([(self.re_begin_opt, self.parse_opt),
                                            (self.re_job_time, self.parse_time)])
            
    def parse_prelude(self):
        m = self.outfile.discard_to_match(self.re_mopac_version)
        self.provenance['version'] = m.group(1).strip()
        
        m = self.outfile.discard_to_match(self.re_prelude_header)
        self.provenance['date'] = m.group(1).strip()
        
        m = self.outfile.discard_to_match(self.re_prelude_keyword)
        self._parse_cache['method'] = m.group(1)
        
        m = self.outfile.discard_to_match(self.re_prelude_keyword)
        self._parse_cache['multiplicity'] = m.group(1)
        
        m = self.outfile.discard_to_match(self.re_system_charge)
        self._parse_cache['charge'] = int(m.group(1))
        
        # Process any remaining echoed keywords
        # Need to do this for PRTXYZ
        m = self.outfile.discard_to_match(self.re_prelude_keyword)
        self._parse_cache['addtl_keywords'] = [m.group(1)]
        while not self.outfile.matches(self.re_section_break):
            line = self.outfile.readline()
            m = self.re_prelude_keyword.search(line)
            if m:
                self._parse_cache['addtl_keywords'].append(m.group(1))
        self.outfile.discard_to_match(self.re_section_break)
        
        self.provenance['keyword_line'] = self.outfile.readline().strip()
        self.provenance['title'] = self.outfile.readline().strip()
        
        self.outfile.skip_blanks()
        # Ready for geometry
        
    def parse_initial_geometry(self):
        '''Parse initial geometry'''
                
             
        '''   ATOM   CHEMICAL          X               Y               Z
  NUMBER   SYMBOL      (ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)
 
     1       C         79.07400000  *  34.56300000  *  48.37700000  *
     2       C         78.84000000  *  33.08600000  *  47.96200000  *
     3       C         78.64700000  *  34.71300000  *  49.84600000  *
     4       C         80.50100000  *  34.90900000  *  48.04300000  *
     5       C         79.90000000  *  32.14800000  *  47.73400000  *
'''
        
        assert 'ATOM' in self.outfile.peek()
        assert 'Y' in self.outfile.peek()
        
        self.outfile.readline()
        self.outfile.readline()
        self.outfile.skip_blanks()
        
        atoms = []
        coords = []
        
        line = self.outfile.readline()
        while line.strip() != '':
            fields = line.split()
            atoms.append(fields[1])
            coords.append((float(fields[2]), float(fields[4]), float(fields[6])))
            
            line = self.outfile.readline()
            
        coords = np.array(coords, dtype=np.float_)
        self.atoms = atoms
        self.initial_coords = coords
        
    def parse_opt(self, m, line):
        '''Parse an optimization sequence'''
        
        '''          Geometry optimization using L-BFGS
 CYCLE:     1 TIME:  13.977 TIME LEFT:  2.00D  GRAD.:   849.028 HEAT: -5754.576
 CYCLE:     2 TIME:   5.523 TIME LEFT:  2.00D  GRAD.:   829.477 HEAT: -5761.702
 CYCLE:     3 TIME:   7.398 TIME LEFT:  2.00D  GRAD.:   625.118 HEAT: -6042.151
 CYCLE:     4 TIME:   7.332 TIME LEFT:  2.00D  GRAD.:   830.612 HEAT: -5761.701
 CYCLE:     5 TIME:   7.359 TIME LEFT:  2.00D  GRAD.:   625.857 HEAT: -6042.301'''
        
        energies = []
        grads = []
        
        while self.outfile.readline(strip=True):
            fields = self.outfile.last_line.split()
            if fields[0] == 'CYCLE:':
                grads.append(float(fields[8]))
                energies.append(float(fields[10]))
                self.get_log_geom()
            elif fields[0] == 'RESTART':
                grads.append(float(fields[7]))
                energies.append(float(fields[9]))
                self.get_log_geom()
            else:
                break            
                          
                
        if self.logfile is not None:
            self.geoms = np.array(self.geoms)
            self.opt_energies = np.array(self.opt_energies)
            self.opt_grads = np.array(self.opt_grads)
        else:
            self.opt_energies = np.array(energies, dtype=np.float_)
            self.opt_grads = np.array(grads, dtype=np.float_)

        self.outfile.discard_to_match(re.compile('HEAT'))        
        if 'PRTXYZ' in self._parse_cache['addtl_keywords']:
            self.parse_prtxyz()
    
    def get_log_geom(self):
        if not self.logfile:
            return
       
        self.logfile.discard_to_match(self.re_energy_updated)
        self.opt_energies.append(ffloat(self.logfile.last_match.group(2)))
        self.logfile.read_and_match(self.re_grad_updated)
        self.opt_grads.append(ffloat(self.logfile.last_match.group(2)))
        self.logfile.read_and_match(self.re_geom_updated)
                
        # Should know an atom count by now
        natoms = len(self.atoms)
        coords = read_positions(self.logfile.textfile, natoms)     
        self.geoms.append(coords)        
        
    def parse_prtxyz(self):
        '''Parse the final geometry in input format'''
        
        lines = []
        lines.append(self.outfile.readline()) # keywords
        lines.append(self.outfile.readline()) # title
        lines.append(self.outfile.readline()) # geom sep
        while self.outfile.readline() != ' \n':
            lines.append(self.outfile.last_line)
        lines.append('\n')
        
        self.next_job_input=''.join(lines)
    
    def parse_time(self, m, line):
        self.provenance['jobtime'] = float(m.group(1))
        