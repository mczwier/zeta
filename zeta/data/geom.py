'''
Created on Oct 31, 2021

@author: mzwier
'''

from collections import deque
import numpy as np
import pandas as pd

from .elements import element_label_to_number

def normalize_atoms(atoms):
    return np.array([element_label_to_number(atom) for atom in atoms])

class Geometry:
    def __init__(self, atoms, coords=None, property_sets = None):
        self.atoms = atoms # atomic numbers
        self.coords = coords # coordinates of atoms
        self.property_sets = property_sets or []
        
    def __repr__(self):
        return '<Geometry at {:#x}, {:d} atoms, first_coords={}>'.format(id(self),
                                                                  len(self.atoms) if self.atoms is not None else 0,
                                                                  self.coords[0] if self.coords is not None else None)
    def has_same_atoms(self, geometry_or_atoms):
        myatoms = np.asarray(self.atoms)
        theiratoms = np.asarray(geometry_or_atoms)
        return (myatoms == theiratoms).all()

class GeometrySequence:
    @staticmethod
    def from_coalesced(atom_data, coords, properties=None, properties_structure='dict_of_lists'):
        '''Create a GeometrySequence directly from pre-compiled sequences of 
        coordinates and corresponding properties. The properties can be a list of
        dicts (properties_structure='list_of_dicts', default) or a dict of lists
        (prroperties_structure='dict_of_lists'). If properties are a dict of lists,
        then missing entries must be explicitly included.'''
        
        gseq = GeometrySequence()
        gseq.atom_data = atom_data
        gseq.coords = coords
        
        if properties_structure == 'list_of_dicts':
            gseq.properties = pd.DataFrame(properties or [])
        elif properties_structure == 'dict_of_lists':
            gseq.properties = pd.DataFrame.from_dict(properties)
        else:
            raise ValueError('invalid property structuring {!r}'.format(properties_structure))
        
        gseq.geometries = [Geometry(atom_data, coords[i], gseq.properties.iloc[i]) 
                           for i in range(len(atom_data))] 
        gseq.finalized = True
        
        return gseq       
    
    def __init__(self, geometries=None):
        self.geometries = geometries or deque()
        self.finalized = False
        
        self.atom_data = None
        self.properties = None
        self.coords = None
                        
    def finalize(self):
        '''Finalize incrementally accumulated coordinates into a (possibly)
        more efficient structure. The ability to prepend or append geometries
        (or do anything other than read, in general) is not guaranteed after
        finalization.'''
        
        self.atoms = self.geometries[0].atoms
        self.coords = [geom.coords for geom in self.geometries]
        
        # Deduplicate coordinates, if we can
        try:
            self.coords = np.array(self.coords)
        except:
            pass
        else:
            for i, geom in enumerate(self.geometries):
                geom.coords = self.coords[i]
        
        # Not sue it's worth deduplicating properties        
        self.properties = pd.DataFrame(geom.properties for geom in self.geometries)
        self.finalized = True

