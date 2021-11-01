'''
Created on Oct 30, 2021

@author: mzwier
'''

import re
from ..parser import TextFileParser, ffloat
from ..data.calc import Calculation, CalcType
from ..data.method import QMMethod
from ..data.geom import Geometry, GeometrySequence, normalize_atoms
from ..data.helpers import normalize_multiplicity


import numpy as np, pandas as pd

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
        self._parse_cache = {} # intermediate information
        self.outfile = TextFileParser(outfile)
        self.logfile = TextFileParser(logfile) if logfile else None
        self.auxfile = TextFileParser(auxfile) if auxfile else None
        
        # Calculation(s)
        self.calculation = Calculation(method=QMMethod())
        self.calculation.provenance['engine'] = 'mopac'
        # A copy stored here so that we can replace it with higher-precision data later if
        # possible
        self.initial_geometry = None
        self.final_geometry = None 
        
        self.geoms = []
        
    def parse(self):
        self.parse_prelude()
        self.parse_initial_geometry()
        
        self.calculation.calc_type = CalcType.ENERGY
    
        # Skip extra Cartesian block if present
        if self.outfile.test_lookahead(re.compile('CARTESIAN COORDINATES'), 5):
            self.outfile.discard_to_match(re.compile('1'))
            self.outfile.skip_to_blank()
    
        while self.outfile.peek():
            self.outfile.scan_and_dispatch([(self.re_begin_opt, self.parse_opt),
                                            (self.re_job_time, self.parse_time)])
            
        if self.calculation.calc_type == CalcType.OPTIMIZATION:
            self.calculation.geometries = self.geoms
        else:
            self.calculation.geometries = [self.initial_geometry]
             
    def parse_prelude(self):
        m = self.outfile.discard_to_match(self.re_mopac_version)
        self.calculation.provenance['version'] = m.group(1).strip()
        
        m = self.outfile.discard_to_match(self.re_prelude_header)
        self.calculation.provenance['date'] = pd.to_datetime(m.group(1).strip())
        
        m = self.outfile.discard_to_match(self.re_prelude_keyword)
        self.calculation.method.name = m.group(1)
        
        m = self.outfile.discard_to_match(self.re_prelude_keyword)
        self.calculation.method.multiplicity = normalize_multiplicity(m.group(1))
        
        m = self.outfile.discard_to_match(self.re_system_charge)
        self.calculation.method.charge = int(m.group(1))
        
        # Process any remaining echoed keywords
        # Need to do this for PRTXYZ
        m = self.outfile.discard_to_match(self.re_prelude_keyword)
        self.calculation.method.additional_keywords.append(m.group(1))
        while not self.outfile.matches(self.re_section_break):
            line = self.outfile.readline()
            m = self.re_prelude_keyword.search(line)
            if m:
                self.calculation.method.additional_keywords.append(m.group(1))
        self.outfile.discard_to_match(self.re_section_break)
        
        self.calculation.provenance['keyword_line'] = self.outfile.readline().strip()
        self.calculation.provenance['title'] = self.outfile.readline().strip()
        
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
    
        while self.outfile.readline(strip=True):
            fields = self.outfile.last_line.strip().split()
            atoms.append(fields[1])
            coords.append((float(fields[2]), float(fields[4]), float(fields[6])))
            
        coords = np.array(coords, dtype=np.float_)
        atoms = normalize_atoms(atoms)
            
        self.initial_geometry = Geometry(atoms, coords)

        
    def parse_opt(self, m, line):
        '''Parse an optimization sequence'''
        
        '''          Geometry optimization using L-BFGS
 CYCLE:     1 TIME:  13.977 TIME LEFT:  2.00D  GRAD.:   849.028 HEAT: -5754.576
 CYCLE:     2 TIME:   5.523 TIME LEFT:  2.00D  GRAD.:   829.477 HEAT: -5761.702
 CYCLE:     3 TIME:   7.398 TIME LEFT:  2.00D  GRAD.:   625.118 HEAT: -6042.151
 CYCLE:     4 TIME:   7.332 TIME LEFT:  2.00D  GRAD.:   830.612 HEAT: -5761.701
 CYCLE:     5 TIME:   7.359 TIME LEFT:  2.00D  GRAD.:   625.857 HEAT: -6042.301'''
        
        self.calculation.calc_type = CalcType.OPTIMIZATION
        
        energies = []
        grads = []
        geoms = []
        
        while self.outfile.readline(strip=True):
            fields = self.outfile.last_line.split()
            if fields[0] == 'CYCLE:':
                grads.append(float(fields[8]))
                energies.append(float(fields[10]))
                geoms.append(self.get_log_geom())
            elif fields[0] == 'RESTART':
                grads.append(float(fields[7]))
                energies.append(float(fields[9]))
                geoms.append(self.get_log_geom())
            else:
                break            
                          
                
        if self.logfile is not None:
            # Logfile contains initial but not final geometry
            self.geoms = geoms
            self.initial_geometry = geoms[0]
        else:
            self.geoms = [self.initial_geometry]
        
        if 'PRTXYZ' in map(str.upper, self.calculation.method.additional_keywords):
            self.parse_prtxyz()
    
    def get_log_geom(self):
        if not self.logfile:
            return
       
        self.logfile.discard_to_match(self.re_energy_updated)
        energy = ffloat(self.logfile.last_match.group(2))
        self.logfile.read_and_match(self.re_grad_updated)
        grad = ffloat(self.logfile.last_match.group(2))
        self.logfile.read_and_match(self.re_geom_updated)
                
        # Should know an atom count by now
        natoms = len(self.initial_geometry.atoms)
        coords = read_positions(self.logfile.textfile, natoms)     
        return Geometry(self.initial_geometry.atoms, coords, {'energy': energy, 'grad': grad})        
        
    def parse_prtxyz(self):
        '''Parse the final geometry in input format'''
        
        
        '''\
          CURRENT BEST VALUE OF HEAT OF FORMATION =  -7309.204109
 PM6-D3H4 Charge=0 Singlet XYZ PRNT=2 PRTXYZ Threads=1 AUX(6,COMP,PRECISION=8,XP,XS,XW)
 T+ target product xyz coordinates with 8 angstroms of solvent molecules

  C    79.35223568 +1  34.18884577 +1  48.50904523 +1
  C    79.16653308 +1  32.68104606 +1  48.36058949 +1
  C    78.82421764 +1  34.68922348 +1  49.85322545 +1
  C    80.75316617 +1  34.66567827 +1  48.12519736 +1
  C    80.24985037 +1  31.80129956 +1  48.29670655 +1
  C    80.02752032 +1  30.44542088 +1  48.02791793 +1
'''    
        
        self.outfile.discard_to_match(re.compile('HEAT OF FORMATION\s+=\s+(\S+)'))
        energy = ffloat(self.outfile.last_match.group(1))
        lines = []
        lines.append(self.outfile.readline()) # keywords
        lines.append(self.outfile.readline()) # title
        lines.append(self.outfile.readline()) # geom sep
        while self.outfile.readline() != ' \n':
            lines.append(self.outfile.last_line)
        lines.append('\n')
        
        # Save input that MOPAC has so nicely prepared for us
        self.next_job_input=''.join(lines)
        natoms = len(self.initial_geometry.atoms)
        coords = np.empty((natoms, 3), dtype=np.float_)
        
        # Parse final geometry
        for i, line in enumerate(lines[3:3+natoms]):
            line = line.strip()
            assert line != ''
            
            fields = line.split()
            coords[i, 0] = float(fields[1])
            coords[i, 1] = float(fields[3])
            coords[i, 2] = float(fields[5])
    
        self.final_geometry = Geometry(self.initial_geometry.atoms, coords, 
                                       properties={'energy': energy,
                                                   'grad': None})
    
    def parse_time(self, m, line):
        self.calculation.provenance['jobtime'] = float(m.group(1))
        