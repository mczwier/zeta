'''
Created on Oct 30, 2021

@author: mzwier
'''

from zeta.parser import ParseError
from zeta.parser.textparser import TextFileParser, ffloat, RegexpMatch, ContainsText, whitespace_only
from zeta.data.provenance import Provenance
from zeta.data.propertyset import PropertySet
from zeta.data.calc import Calculation, CalcType
from zeta.data.method import Method
from zeta.data.geom import Geometry, normalize_atoms
from zeta.data.helpers import normalize_multiplicity


import numpy as np, pandas as pd

# Much faster than numpy.loadtxt
def _read_positions_py(source, natoms):
    coords = np.empty((natoms, 3), dtype=np.float64)
    
    # field length and boundary positions
    fieldlen = 0
    b1 = 0
    b2 = 0
    b3 = 0
    
    for iatom in range(0, natoms):
        line = source.readline()
        if not fieldlen:
            fieldlen = (len(line)-1)//3
            b1 = fieldlen
            b2 = 2*fieldlen
            b3 = 3*fieldlen

        coords[iatom,0] = float(line[0:b1])
        coords[iatom,1] = float(line[b1:b2])
        coords[iatom,2] = float(line[b2:b3])
        
    return coords
read_positions = _read_positions_py

#cython is a factor of three improvement
try:
    import pyximport
    pyximport.install(inplace=True, language_level=3)
    from ._parsext import _assign_positions_cy

    def _read_positions_cy(source, natoms):
        coords = np.empty((natoms,3), dtype=np.float64)
        lines = [source.readline() for _i in range(natoms) ]
        fieldwidth = (len(lines[0]) - 1)//3
        _assign_positions_cy(lines, natoms, fieldwidth, coords)
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
    m_section_break = ContainsText('***')
    re_mopac_version = RegexpMatch(r'Version: (\S+)')
    re_prelude_header = RegexpMatch(r'CALCULATION DONE:\s+([^\*]+)')
    re_prelude_keyword = RegexpMatch(r'\*\s+(\S+)')
    re_system_charge = RegexpMatch(r'CHARGE ON SYSTEM = (\S+)')


    # Optimization
    re_begin_opt = RegexpMatch(r'Geometry optimization using (\S+)')
    
    # FLEPO
    re_opt_cycle_begin = RegexpMatch(r'AT THE BEGINNING OF CYCLE\s*(\S+)\s*THE FUNCTION VALUE IS\s*(\S+)')    
    re_opt_cycle_gradnorm = RegexpMatch(r'GRADIENT NORM\s*=\s*(\S+)')
    re_flepo_coord_line = RegexpMatch('XPARAM')
    re_flepo_grad_line = RegexpMatch('GRAD')
    
    
    # TIME LEFT:  1.92D  GRAD.:    29.945 HEAT: -7294.364
    re_opt_cycle_end = RegexpMatch(r'GRAD\.:\s*(\S+)\s+HEAT:\s*(\S+)')
    
    # Aux file
    # HEAT_OF_FORM_UPDATED:KCAL/MOL=-0.57545763423410D+04
    re_aux_geom = RegexpMatch('ATOM_X')
    re_energy_updated = RegexpMatch(r'HEAT_OF_FORM_UPDATED:([^=]+)=(\S+)') 
    re_grad_updated = RegexpMatch(r'GRADIENT_UPDATED:([^=]+)=(\S+)')
    re_geom_updated = RegexpMatch(r'ATOM_X_UPDATED')
    
    re_final_energy = RegexpMatch(r'FINAL HEAT OF FORMATION =\s*(\S+)')

    re_job_time = RegexpMatch(r'JOB TIME:\s+(\S+)')
    
    def __init__(self, outfile, auxfile=None, auxnfile=None):
        self.outfile = TextFileParser(outfile)
        self.auxfile = TextFileParser(auxfile) if auxfile else None
        
        # Confusingly, an AUX-format geometry trace can be put in a DIFFERENT
        # aux file. So for some jobs, auxfile and auxnfile are the same (only AUX
        # is specified, without a file number), and in others, they may be different
        # [AUX(n,...) is specified]
        self.auxnfile = TextFileParser(auxnfile) if auxnfile else self.auxfile
        
        # Calculation(s)
        self.provenance = Provenance()
        self.method = Method()
        self.calculation = Calculation()
        
        self.provenance['engine'] = 'mopac'

        self.initial_geometry = None
        self.final_geometry = None 
        self.geoms = []
        
    def parse(self):
        self.parse_prelude()
        self.calculation.calc_type = CalcType.ENERGY
        
        self.outfile.discard_while_match(whitespace_only)
        
        self.initial_geometry = self.parse_geometry()
        
        
        # Skip extra Cartesian block if present, as it likely has lower
        # precision than the geometry just read
        if self.outfile.testp_within_next(ContainsText('CARTESIAN COORDINATES'), 5):
            self.outfile.discard_until_match(ContainsText('1'))
            self.outfile.discard_until_match(whitespace_only)
    
        while self.outfile.line:
            self.outfile.scan_and_dispatch([(self.re_begin_opt, self.parse_opt),
                                            (self.re_final_energy, 
                                             self.parse_scf_result_auxfile if self.auxfile 
                                             else self.parse_scf_result_outfile),
                                            (self.re_job_time, self.parse_time)])
            
        self.finalize_calculation()    
        
    def finalize_calculation(self):
        if self.calculation.calc_type == CalcType.OPTIMIZATION:
            self.calculation.geometries = [self.initial_geometry] + self.geoms 
            if self.final_geometry:
                self.calculation.geometries.append(self.final_geometry)
        else:
            if self.final_geometry:
                self.calculation.geometries = [self.final_geometry]
            else:
                self.calculation.geometries = [self.initial_geometry]
            
        # Same run gets the same provenance and method info
        for geom in self.calculation.geometries:
            for ps in geom.property_sets:
                ps.method = self.method
                ps.provenance = self.provenance
        
    def parse_prelude(self):
        self.outfile.discard_until_match(self.re_mopac_version)
        self.provenance['version'] = self.outfile.presult[1].strip()
        
        self.outfile.discard_until_match(self.re_prelude_header)
        self.provenance['date'] = pd.to_datetime(self.outfile.presult[1].strip())
        
        self.outfile.discard_until_match(self.re_prelude_keyword)
        self.method['name'] = self.outfile.presult[1]
        
        self.outfile.discard_until_match(self.re_prelude_keyword)
        self.method['multiplicity'] = normalize_multiplicity(self.outfile.presult[1])
        
        self.outfile.discard_until_match(self.re_system_charge)
        self.method['charge'] = int(self.outfile.presult[1].strip())
        
        # Process any remaining echoed keywords
        additional_keywords = []
        self.outfile.discard_until_match(self.re_prelude_keyword)
        additional_keywords.append(self.outfile.presult[1])
        
        while self.outfile.read_until_match(self.m_section_break):
            if self.outfile.testp(self.re_prelude_keyword):
                additional_keywords.append(self.outfile.presult[1])
        self.method['additional_keywords'] = additional_keywords
        
        
        self.outfile.nextline()
        self.provenance['keyword_line'] = self.outfile.line.strip()
        
        self.outfile.nextline()
        self.provenance['title'] = self.outfile.line.strip()
        
    def parse_geometry(self):
        '''Parse initial geometry'''
                
             
        # Cartesian coordinates
        '''\
   ATOM   CHEMICAL          X               Y               Z
  NUMBER   SYMBOL      (ANGSTROMS)     (ANGSTROMS)     (ANGSTROMS)
 
     1       C         79.07400000  *  34.56300000  *  48.37700000  *
     2       C         78.84000000  *  33.08600000  *  47.96200000  *
     3       C         78.64700000  *  34.71300000  *  49.84600000  *
     4       C         80.50100000  *  34.90900000  *  48.04300000  *
     5       C         79.90000000  *  32.14800000  *  47.73400000  *
'''
   
        # Internal coordinates
        '''\
  ATOM    CHEMICAL      BOND LENGTH      BOND ANGLE     TWIST ANGLE 
 NUMBER    SYMBOL       (ANGSTROMS)      (DEGREES)       (DEGREES) 
   (I)                     NA:I           NB:NA:I       NC:NB:NA:I       NA    NB    NC 
     1       C          0.00000000       0.0000000       0.0000000        0     0     0
     2       C          1.52266930  *    0.0000000       0.0000000        1     0     0
     3       C          1.39620180  *  120.6151100  *    0.0000000        2     1     0
     4       C          1.39240400  *  120.1550500  *  177.5468800  *     3     2     1
     5       C          1.39378300  *  120.0936200  *   -0.5905398  *     4     3     2
     6       C          1.39331810  *  119.8763900  *   -0.4477492  *     5     4     3
     7       C          1.39283760  *  120.0707000  *    0.4505725  *     6     5     4
     8       H          1.09211460  *  119.7804000  * -179.0000900  *     7     6     5
'''     
                
        if 'ATOM' not in self.outfile.line:
            self.outfile.discard_until_match(ContainsText('ATOM'))
                
        if 'Y' in self.outfile.line:
            return self.parse_geom_cartesian()
        elif 'ANGLE' in self.outfile.line:
            return self.parse_geom_internal()
            
        # We've scanned through the output file to the appropriate point, but
        # if the aux file is present, we can use it to get higher precision
        # if the aux file was written with sufficient precision.
        # Maybe add a "prefer aux file" option
        
    def parse_geom_cartesian(self):
        self.outfile.discardn(2)
        self.outfile.discard_while_match(whitespace_only)
        
        atoms = []
        coords = []

        while True:
            fields = self.outfile.line.strip().split()
            atoms.append(fields[1])
            coords.append((float(fields[2]), float(fields[4]), float(fields[6])))
            
            if self.outfile.read_until_match(whitespace_only):
                continue
            else:
                break
            
        coords = np.array(coords, dtype=np.float_)
        atoms = normalize_atoms(atoms)
            
        return Geometry(atoms, coords)

    def parse_geom_internal(self):
        # for now, skip to Cartesians
        
        self.outfile.discard_until_match(ContainsText('CARTESIAN'))
        self.outfile.discardn(3)
        
        atoms = []
        coords = []
        while self.outfile.read_until_match(whitespace_only):
            fields = self.outfile.line.strip().split()
            atoms.append(fields[1])
            coords.append((float(fields[2]), float(fields[3]), float(fields[4])))

        coords = np.array(coords, dtype=np.float_)
        atoms = normalize_atoms(atoms)
            
        return Geometry(atoms, coords)
            
        
    def parse_opt(self, _parser):
        '''Parse an optimization sequence'''
        self.calculation.calc_type = CalcType.OPTIMIZATION
        self.method['opt_algorithm'] = self.outfile.presult[1]
               
        # Several different output formats are possible:
        # FLEPO reports everything in the output file
        # non-FLEPO reports everything in the log file
        # Who knows where the aux file fits in?
        
        if self.outfile.testp_within_next(self.re_opt_cycle_begin, 7):
            self.parse_opt_flepo()
        elif self.outfile.testp_within_next(ContainsText('CYCLE'), 5):
            self.outfile.discard_until_match(ContainsText('CYCLE'))
            self.parse_opt_abbrev()
        else:
            # TODO this should be a warning
            raise ParseError('incomplete/invalid geometry optimization output')
        
        
    def parse_opt_flepo(self):
        '''Parse an optimization with intermediate results in FLEPO format'''
        
        '''\
 AT THE BEGINNING OF CYCLE    1  THE FUNCTION VALUE IS    169.213512
  THE CURRENT POINT IS ...
  GRADIENT NORM =   262.5702
  ANGLE COSINE =    0.3110
  THE CURRENT POINT IS ...


    I           1          2          3          4          5          6
  XPARAM(I)    0.0100     0.0100    -0.0100     1.5171     0.0100    -0.0100
  GRAD  (I)   -2.6168    -0.2773    15.2831    -7.5808     5.4951     4.1743
  PVECT (I)    0.008330   0.000183  -0.013098   0.040203  -0.003942  -0.022138


    I           7          8          9         10         11         12
  XPARAM(I)    2.2082     1.2072    -0.0100     3.5906     1.2072    -0.0100
  GRAD  (I)    4.3132    -4.4019    -1.1382   -11.7968    -5.2627    -0.9246
  PVECT (I)   -0.009302   0.007961   0.008486   0.062562   0.008574   0.008385
'''
        
        # above: beginning of cycle, below: end of cycle
        
        '''\
    I         103        104        105        106        107        108
  XPARAM(I)   -0.4983    -0.7413     1.1490     0.2465    -1.1342     1.9429
  GRAD  (I)  143.1710   -50.1257   119.0323   -96.0928    46.1118  -100.2140
  PVECT (I)   -0.019338   0.013849  -0.015626   0.015420  -0.013403   0.013869


    I         109        110        111
  XPARAM(I)    0.9712    -1.5072     2.7170
  GRAD  (I)  -36.4285    30.8778   -59.8688
  PVECT (I)    0.044993  -0.073390   0.061366
  -ALPHA.P.G =         26.345300


           NUMBER OF COUNTS =     3         COS    =     0.3110
  ABSOLUTE  CHANGE IN X     =     0.115902  ALPHA  =     0.3593
  PREDICTED CHANGE IN F     =   -26.35      ACTUAL =   -4.502    
  GRADIENT NORM             =    231.7    


 CYCLE:     1 TIME:   0.211 TIME LEFT:  2.00D  GRAD.:   231.723 HEAT:  164.7111
'''
               
        # The initial point of the optimization is reflected in the output file
        # but not the aux file. The optimization initial point appears to be 
        # adjusted from the initial coordinates.
        # We have to get the very first set of the coordinates (and associated
        # energy and gradient norm) from the output file, and then we can switch
        # to the aux file for geometries, energies, and gradient norms, but the
        # gradient itself appears to have to come from the output file.
        
        initial_point = True
        while self.outfile.testp_within_next(self.re_opt_cycle_begin, 11):
            self.geoms.append(self.parse_flepo_structure(initial_point=initial_point))
            initial_point = False        
        
        
    def parse_flepo_structure(self, initial_point):
        '''Parse a structure from output and (if applicable) aux file.
        The first point is not stored in the aux file, so initial_point==True
        prevents attempting to read the geometry from the aux file.
        In all cases, the gradient is read from the output file'''
        
        
        
        all_coords = []
        all_grads = []
        natoms = len(self.initial_geometry.atoms)
        
        # Process coordinates in output file if we are on the first point
        # (for which the aux file has no structure) or if there is no aux file
        process_outfile_coords = True
        if initial_point is True:
            process_outfile_coords = True
        else:
            if self.auxnfile:
                process_outfile_coords = False
            else:
                process_outfile_coords = True

        self.outfile.discard_until_match(self.re_opt_cycle_begin)
        energy = ffloat(self.outfile.presult[2])
        self.outfile.discard_until_match(self.re_opt_cycle_gradnorm)
        grad_norm = ffloat(self.outfile.presult[1])
        
        
        while True:    
            self.outfile.discard_until_match(self.re_flepo_coord_line)
            
            # Read geometry if it's not in an aux file    
            if process_outfile_coords:
                fields = self.outfile.line.split()[1:]
                all_coords.extend(float(field) for field in fields)
            
            # Read gradient
            self.outfile.read_and_assertp(self.re_flepo_grad_line)
            fields = self.outfile.line.split()[2:]
            all_grads.extend(float(field) for field in fields)
            
            # Read and discard PVECT
            self.outfile.nextline()
            
            nfields = len(all_grads)
            if nfields == natoms*3:
                break

        grad = np.array(all_grads).reshape(natoms, 3)        
        if process_outfile_coords:
            coords = np.array(all_coords).reshape(natoms, 3)
            geom = Geometry(atoms=self.initial_geometry.atoms,
                            coords=coords)
            ps = PropertySet()
            geom.property_sets.append(ps)
            ps['energy'] = energy
            ps['gradnorm'] = grad_norm
        else:
            # Get coordinates from auxnfile
            geom = self.parse_aux_geom(self.auxnfile, geom_type='opt')
        
        geom.property_sets[0]['gradient'] = grad
        return geom
        
    def parse_opt_abbrev(self):
        '''Parse an optimization in standard output format, with or without
        detailed log information.'''
        
        '''          Geometry optimization using L-BFGS
 CYCLE:     1 TIME:  13.977 TIME LEFT:  2.00D  GRAD.:   849.028 HEAT: -5754.576
 CYCLE:     2 TIME:   5.523 TIME LEFT:  2.00D  GRAD.:   829.477 HEAT: -5761.702
 CYCLE:     3 TIME:   7.398 TIME LEFT:  2.00D  GRAD.:   625.118 HEAT: -6042.151
 CYCLE:     4 TIME:   7.332 TIME LEFT:  2.00D  GRAD.:   830.612 HEAT: -5761.701
 CYCLE:     5 TIME:   7.359 TIME LEFT:  2.00D  GRAD.:   625.857 HEAT: -6042.301'''
        
        
        energies = []
        grads = []
        geoms = []
                
        while self.outfile.read_until_match(whitespace_only):
            fields = self.outfile.line.split()
            if fields[0] == 'CYCLE:':
                grads.append(float(fields[8]))
                energies.append(float(fields[10]))
                if self.auxnfile:
                    geoms.append(self.parse_aux_geom(self.auxnfile, geom_type='opt'))
            elif fields[0] == 'RESTART':
                grads.append(float(fields[7]))
                energies.append(float(fields[9]))
                if self.auxnfile:
                    geoms.append(self.parse_aux_geom(self.auxnfile, geom_type='opt'))
            else:
                break            
        
        self.geoms = geoms      
                            
    def parse_aux_geom(self, parser, geom_type):
        
        # Should know an atom count by now
        natoms = len(self.initial_geometry.atoms)
        
        # Properties/results for this structure
        properties = PropertySet()
        if geom_type == 'opt':
            # Optimization point
            parser.discard_until_match(self.re_energy_updated)
            properties['energy'] = ffloat(parser.presult[2])
        
        
            parser.read_and_assertp(self.re_grad_updated)
            properties['grad_norm'] = ffloat(parser.presult[2])
        
            parser.read_and_assertp(self.re_geom_updated)
        elif geom_type == 'final':
            # HEAT_OF_FORMATION:KCAL/MOL=+0.15193345362744D+03
            parser.discard_until_match(RegexpMatch(r'HEAT_OF_FORMATION:([^=]+)=(\S+)'))
            properties['energy'] = ffloat(parser.presult[2])
            parser.discard_until_match(RegexpMatch(r'GRADIENT_NORM:([^=]+)=(\S+)'))
            properties['grad_norm'] = ffloat(parser.presult[2]) 
            
            parser.discard_until_match(ContainsText('ATOM_X_OPT:'))
                

        assert parser.buffer_is_empty(), 'cannot raw read with non-empty buffer'
        coords = read_positions(parser.textfile, natoms)
        # Parser has had read natoms lines without knowing about it
        parser.fast_forward_linenum(natoms)
        
        return Geometry(self.initial_geometry.atoms, coords, [properties])        
        
    def parse_scf_result_outfile(self, _):
        
        energy = ffloat(self.outfile.presult[1])
        
        self.outfile.discard_until_match(ContainsText('COMPUTATION'))
        
        self.final_geometry = self.parse_geometry()
        try:
            self.final_geometry.property_sets[0]['energy'] = energy
        except IndexError:
            self.final_geometry.property_sets.append(PropertySet({'energy': energy}))

    def parse_scf_result_auxfile(self, _):
        self.final_geometry = self.parse_aux_geom(self.auxfile, geom_type='final')    
    
    def parse_time(self, _):
        self.provenance['jobtime'] = float(self.outfile.presult[1])
        
        