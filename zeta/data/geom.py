'''
Created on Oct 31, 2021

@author: mzwier
'''

import numpy as np

import zeta
from .elements import element_label_to_number

def normalize_atoms(atoms):
    return np.array([element_label_to_number(atom) for atom in atoms], 
                    dtype=zeta.data.atomicnumber_dtype)

class Geometry:
    def __init__(self, atoms, coords=None, property_sets = None):
        self.atoms = atoms # atomic numbers
        self.coords = coords # coordinates of atoms
        self.property_sets = property_sets or []
        
    def __repr__(self):
        return '<Geometry at {:#x}, {:d} atoms, first_coords={}>'\
                .format(id(self),
                        len(self.atoms) if self.atoms is not None else 0,
                        self.coords[0] if self.coords is not None else None)
                
    def has_same_atoms(self, geometry_or_atoms):
        myatoms = np.asarray(self.atoms)
        theiratoms = np.asarray(geometry_or_atoms)
        return (myatoms == theiratoms).all()
